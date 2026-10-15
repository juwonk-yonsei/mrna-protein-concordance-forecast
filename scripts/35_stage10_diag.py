"""Rebuild RPPA and ProCan F1s features with the freeze-commit code.

Reads the project data. Writes only under /tmp/b10_freeze_out.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, "/tmp/b10_freeze")
sys.path.insert(0, "/HDD8T3/kangjw/Cisplatin/2026_NewProj/proj_B_adc_protein")

import importlib

import numpy as np
import pandas as pd
from lightgbm import Booster

REAL = Path("/HDD8T3/kangjw/Cisplatin/2026_NewProj/proj_B_adc_protein")
OUT = Path("/tmp/b10_freeze_out")
OUT.mkdir(parents=True, exist_ok=True)

stage3 = importlib.import_module("scripts.08_stage3")
stage3.ROOT = REAL
stage3.RES = REAL / "results"
stage3.RAW = REAL / "data" / "raw"
import scripts.b2_lib as b2
b2.ROOT = REAL
b2.RES = stage3.RES
b2.RAW = stage3.RAW

genes, p_genes, adc = stage3.load_lists()
mapped = pd.read_csv(stage3.RES / "tcga_rna_id_map_b3.tsv", sep="\t", dtype=str)
rho = stage3.defined_rho(p_genes)
pset = set(p_genes)
adc_loco = pd.read_csv(stage3.RES / "adc_loco.tsv", sep="\t")
for rec in adc_loco.itertuples(index=False):
    if rec.gene_symbol in pset or rec.gene_symbol not in set(adc):
        continue
    if int(rec.n) >= stage3.MIN_N and np.isfinite(rec.s0):
        rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
_stats, mean_pct, sd_pct, _ = stage3.cptac_percentile_tables(genes)
rates = pd.read_csv(stage3.RES / "cptac_tumor_protein_rate.tsv", sep="\t")
rate = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
rho = stage3.forecast_rho(rho, genes)
gmean = stage3.gmean_all(rho, p_genes)

samples, _path = stage3.disease_samples()
print("reading TCGA RNA", flush=True)
mat = stage3.load_tcga_matrix(mapped, samples)
antibodies = pd.read_csv(stage3.RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
model = pd.read_parquet(stage3.RES / "ar_final_models.parquet", columns=["gene_symbol"])
kept = antibodies[antibodies["gene_symbol"].isin(set(model["gene_symbol"]))].copy()
ab_for_gene = {r.gene_symbol: r.xena_antibody for r in kept.itertuples(index=False)}
rows = []
for disease, ids in samples.items():
    mu, sd = stage3.gene_mean_sd(mat, ids)
    mu_pct = stage3.percentile_map(mu)
    sd_map = stage3.percentile_map(sd)
    for gene in antibodies["gene_symbol"]:
        ext_mu_pct = mu_pct.get(gene, (np.nan, np.nan))[1]
        ext_sd_pct = sd_map.get(gene, (np.nan, np.nan))[1]
        rec = stage3.f1s_feature_row(gene, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, ext_sd_pct)
        antibody = ab_for_gene.get(gene)
        if antibody is None:
            antibody = antibodies.loc[antibodies.gene_symbol == gene, "xena_antibody"].iloc[0]
        rec.update({"xena_antibody": antibody, "gene_symbol": gene, "cancer": disease})
        rows.append(rec)
    print("rppa", disease, flush=True)
rppa = pd.DataFrame(rows)
rppa.to_csv(OUT / "rppa_features_freeze.tsv", sep="\t", index=False)
del mat

print("reading Sanger RNA", flush=True)
lin_table = pd.read_csv(stage3.RES / "procan_sanger_by_lineage.tsv", sep="\t")
lineages = lin_table.loc[lin_table["n_ge_30"].astype(str).isin(["True", "true", "1"]), "Tissue_type"].tolist()
mapping = pd.read_csv(stage3.RAW / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t", dtype=str)
sample_ids = set(pd.read_csv(stage3.RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
mapping = mapping[mapping["Project_Identifier"].isin(sample_ids)]
cands = mapping[mapping["Tissue_type"].isin(lineages)]["SIDM"].tolist()
gindex = {g: i for i, g in enumerate(genes)}
mindex = {m: i for i, m in enumerate(cands)}
arr = np.full((len(cands), len(genes)), np.nan, dtype=np.float32)
filled = np.zeros((len(cands), len(genes)), dtype=bool)
sanger_ids = set()
path = stage3.RAW / "cmp" / "rnaseq_merged" / "rnaseq_merged_20260323.csv"
import csv
with open(path, newline="") as handle:
    reader = csv.reader(handle)
    header = next(reader)
    i_model = header.index("model_id")
    i_gene = header.index("gene_symbol")
    i_tpm = header.index("rsem_tpm")
    i_source = header.index("data_source")
    for parts in reader:
        if len(parts) <= i_source or parts[i_source] != "Sanger":
            continue
        model_id = parts[i_model]
        sanger_ids.add(model_id)
        if model_id not in mindex or parts[i_gene] not in gindex:
            continue
        if filled[mindex[model_id], gindex[parts[i_gene]]]:
            continue
        try:
            tpm = float(parts[i_tpm])
        except ValueError:
            continue
        if not np.isfinite(tpm):
            continue
        arr[mindex[model_id], gindex[parts[i_gene]]] = np.float32(np.log2(tpm + 1.0))
        filled[mindex[model_id], gindex[parts[i_gene]]] = True
use_rows = [i for i, sidm in enumerate(cands) if sidm in sanger_ids]
use_ids = [cands[i] for i in use_rows]
arr = arr[use_rows]
sidm_lineage = dict(zip(mapping["SIDM"], mapping["Tissue_type"]))
proteins = pd.read_csv(stage3.RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
proteins = proteins[(proteins["symbol"] != "NA") & (proteins["symbol"] != "") & proteins["symbol"].isin(set(genes))]
pro_rows = []
for lin in lineages:
    idx = [i for i, sidm in enumerate(use_ids) if sidm_lineage.get(sidm) == lin]
    block = arr[idx]
    with np.errstate(all="ignore"):
        mu = np.nanmean(block, axis=0)
        sd = np.nanstd(block, axis=0, ddof=stage3.DDOF)
    mu_d = {genes[i]: float(mu[i]) if np.isfinite(mu[i]) else np.nan for i in range(len(genes))}
    sd_d = {genes[i]: float(sd[i]) if np.isfinite(sd[i]) else np.nan for i in range(len(genes))}
    mu_pct = stage3.percentile_map(mu_d)
    sd_map = stage3.percentile_map(sd_d)
    for row in proteins.itertuples(index=False):
        ext_mu_pct = mu_pct.get(row.symbol, (np.nan, np.nan))[1]
        ext_sd_pct = sd_map.get(row.symbol, (np.nan, np.nan))[1]
        rec = stage3.f1s_feature_row(row.symbol, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, ext_sd_pct)
        rec.update({"protein_id": row.protein_id, "gene_symbol": row.symbol, "lineage": lin})
        pro_rows.append(rec)
    print("procan", lin, len(idx), flush=True)
pro = pd.DataFrame(pro_rows)
pro.to_csv(OUT / "procan_features_freeze.tsv", sep="\t", index=False)

meta = json_load = __import__("json").loads((REAL / "release" / "models" / "F1s_lightgbm_features.json").read_text())
booster = Booster(model_file=str(REAL / "release" / "models" / "F1s_lightgbm.txt"))
med = np.asarray(meta["median_fill"], dtype=np.float64)
names = meta["features"]

def predict(frame):
    x = frame[names].to_numpy(dtype=np.float64)
    filled_x = np.array(x, dtype=np.float64, copy=True)
    inds = np.where(~np.isfinite(filled_x))
    if inds[0].size:
        filled_x[inds] = np.take(med, inds[1])
    return booster.predict(filled_x)

rppa["pred"] = predict(rppa)
pro["pred"] = predict(pro)
frozen_r = pd.read_csv(REAL / "results" / "frozen_forecasts_rppa_F1s.tsv", sep="\t")
frozen_p = pd.read_csv(REAL / "results" / "frozen_forecasts_procan_F1s.tsv", sep="\t")
cache_r = pd.read_csv(REAL / "results" / "stage9" / "cache_rppa_features.tsv", sep="\t")
cache_p = pd.read_csv(REAL / "results" / "stage9" / "cache_procan_features.tsv", sep="\t")

def pred_gap(new, old, keys):
    m = new.merge(old, on=keys, suffixes=("_new", "_old"))
    a = m["pred"].to_numpy(dtype=np.float64) if "pred" in m.columns else m["f1s_lgbm_new"].to_numpy(dtype=np.float64)
    # after merge, pred may be unsuffixed if only on the left
    if "pred" in m.columns:
        a = m["pred"].to_numpy(dtype=np.float64)
    else:
        a = m["pred_new"].to_numpy(dtype=np.float64)
    b = m["f1s_lgbm"].to_numpy(dtype=np.float64)
    both = np.isfinite(a) & np.isfinite(b)
    return float(np.max(np.abs(a[both] - b[both]))), int(len(m)), int(both.sum())

def feat_gap(new, old, keys):
    m = new.merge(old, on=keys, suffixes=("_new", "_old"))
    rows = []
    changed = np.zeros(len(m), dtype=bool)
    for name in names:
        a = m[f"{name}_new"].to_numpy(dtype=np.float64)
        b = m[f"{name}_old"].to_numpy(dtype=np.float64)
        both = np.isfinite(a) & np.isfinite(b)
        one = np.isfinite(a) != np.isfinite(b)
        diff = np.zeros(len(m), dtype=np.float64)
        diff[both] = np.abs(a[both] - b[both])
        changed |= (diff > 1e-8) | one
        rows.append({"feature": name, "max_abs": float(diff.max()) if diff.size else np.nan, "n_differ_1e-8": int(((diff > 1e-8) | one).sum())})
    return pd.DataFrame(rows), int(changed.sum()), int(len(m))

rg, rn, rf = pred_gap(rppa, frozen_r, ["xena_antibody", "cancer"])
pg, pn, pf = pred_gap(pro, frozen_p, ["protein_id", "lineage"])
rf_tab, rc, rm = feat_gap(rppa, cache_r, ["xena_antibody", "gene_symbol", "cancer"])
pf_tab, pc, pm = feat_gap(pro, cache_p, ["protein_id", "gene_symbol", "lineage"])
summary = pd.DataFrame([
    {"item": "rppa_pred_max_abs", "value": rg, "n_merged": rn, "n_finite": rf},
    {"item": "procan_pred_max_abs", "value": pg, "n_merged": pn, "n_finite": pf},
    {"item": "rppa_feature_rows_differ", "value": rc, "n_merged": rm, "n_finite": rm},
    {"item": "procan_feature_rows_differ", "value": pc, "n_merged": pm, "n_finite": pm},
])
summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)
rf_tab.to_csv(OUT / "rppa_feature_gap.tsv", sep="\t", index=False)
pf_tab.to_csv(OUT / "procan_feature_gap.tsv", sep="\t", index=False)
print(summary.to_string(index=False), flush=True)
print(rf_tab.to_string(index=False), flush=True)
print(pf_tab.to_string(index=False), flush=True)
print("diag done", flush=True)
