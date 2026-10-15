"""Refit F1s-LightGBM for a public model file and compare it with frozen predictions.

Attempt 1 uses the frozen LightGBM parameters. Attempt 2 adds deterministic=True
only if attempt 1 misses the 1e-6 absolute-difference limit. num_threads is not set.
Confirmation outputs and frozen prediction files are not rewritten.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import importlib

import numpy as np
import pandas as pd
import yaml
from lightgbm import LGBMRegressor, Booster
from scipy.stats import rankdata

stage3 = importlib.import_module("scripts.08_stage3")
lookup = importlib.import_module("scripts.25_stage8_lookup")

ROOT = stage3.ROOT
RES = stage3.RES
RAW = stage3.RAW
OUT = ROOT / "results" / "stage9"
REL = ROOT / "release" / "models"
F1S_FEATURES = stage3.F1S_FEATURES
LIMIT = 1e-6


def frozen_params():
    b2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
    stage3.B2 = b2
    stage3.B2_LGBM = dict(b2["reliability"]["lightgbm"])
    return stage3.lgbm_params()


def training_matrix():
    preds = pd.read_csv(RES / "f1s_predictions.tsv", sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    train = preds[preds["in_P"] & np.isfinite(preds["actual"])].copy()
    x = train[F1S_FEATURES].to_numpy(dtype=np.float64)
    y = train["actual"].to_numpy(dtype=np.float64)
    x_f, _, _ = stage3.fill_median(x, x[:1])
    return train, x_f, y


def fit_lgbm(params, x_f, y):
    model = LGBMRegressor(**params)
    model.fit(x_f, y)
    return {"lgbm": model, "x_tr": x_f}


def predict(fit, frame):
    x = frame[F1S_FEATURES].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    med = np.nanmedian(fit["x_tr"], axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(med, inds[1])
    return fit["lgbm"].predict(filled)


def rho_tables(genes, p_genes, adc):
    rho = stage3.defined_rho(p_genes)
    pset = set(p_genes)
    adc_loco = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    for rec in adc_loco.itertuples(index=False):
        if rec.gene_symbol in pset or rec.gene_symbol not in set(adc):
            continue
        if int(rec.n) >= stage3.MIN_N and np.isfinite(rec.s0):
            rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
    rho = stage3.forecast_rho(rho, genes)
    _stats, mean_pct, sd_pct, _ = stage3.cptac_percentile_tables(genes)
    rates = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    rate = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
    gmean = stage3.gmean_all(rho, p_genes)
    return rho, mean_pct, sd_pct, rate, gmean


def rppa_features(rho, mean_pct, sd_pct, rate, gmean):
    mapped = pd.read_csv(RES / "tcga_rna_id_map_b3.tsv", sep="\t", dtype=str)
    samples, _path = stage3.disease_samples()
    print("reading TCGA RNA", flush=True)
    mat = stage3.load_tcga_matrix(mapped, samples)
    antibodies = pd.read_csv(RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    antibodies = antibodies[antibodies["kept_for_forecast"].astype(str).isin(["True", "true"])]
    ab_for_gene = {r.gene_symbol: r.xena_antibody for r in antibodies.itertuples(index=False)}
    rows = []
    for disease, ids in samples.items():
        mu, sd = stage3.gene_mean_sd(mat, ids)
        mu_pct = stage3.percentile_map(mu)
        sd_map = stage3.percentile_map(sd)
        for gene in antibodies["gene_symbol"]:
            ext_mu_pct = mu_pct.get(gene, (np.nan, np.nan))[1]
            ext_sd_pct = sd_map.get(gene, (np.nan, np.nan))[1]
            rec = stage3.f1s_feature_row(gene, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, ext_sd_pct)
            rec.update({"xena_antibody": ab_for_gene[gene], "gene_symbol": gene, "cancer": disease})
            rows.append(rec)
    return pd.DataFrame(rows)


def procan_features(genes, rho, mean_pct, sd_pct, rate, gmean):
    lin_table = pd.read_csv(RES / "procan_sanger_by_lineage.tsv", sep="\t")
    lineages = lin_table.loc[lin_table["n_ge_30"].astype(str).isin(["True", "true", "1"]), "Tissue_type"].tolist()
    mapping = pd.read_csv(RAW / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t", dtype=str)
    sample_ids = set(pd.read_csv(RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    mapping = mapping[mapping["Project_Identifier"].isin(sample_ids)]
    cands = mapping[mapping["Tissue_type"].isin(lineages)]["SIDM"].tolist()
    gindex = {g: i for i, g in enumerate(genes)}
    mindex = {m: i for i, m in enumerate(cands)}
    arr = np.full((len(cands), len(genes)), np.nan, dtype=np.float32)
    filled = np.zeros((len(cands), len(genes)), dtype=bool)
    sanger_ids = set()
    path = RAW / "cmp" / "rnaseq_merged" / "rnaseq_merged_20260323.csv"
    print("reading Sanger RNA", flush=True)
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
            if model_id not in mindex:
                continue
            gene = parts[i_gene]
            if gene not in gindex:
                continue
            if filled[mindex[model_id], gindex[gene]]:
                continue
            try:
                tpm = float(parts[i_tpm])
            except ValueError:
                continue
            if not np.isfinite(tpm):
                continue
            arr[mindex[model_id], gindex[gene]] = np.float32(np.log2(tpm + 1.0))
            filled[mindex[model_id], gindex[gene]] = True
    use_rows = [i for i, sidm in enumerate(cands) if sidm in sanger_ids]
    use_ids = [cands[i] for i in use_rows]
    arr = arr[use_rows]
    sidm_lineage = dict(zip(mapping["SIDM"], mapping["Tissue_type"]))
    proteins = pd.read_csv(RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    proteins = proteins[(proteins["symbol"] != "NA") & (proteins["symbol"] != "") & proteins["symbol"].isin(set(genes))]
    rows = []
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
            rows.append(rec)
    return pd.DataFrame(rows)


def align(new, old, keys):
    merged = new.merge(old[keys + ["f1s_lgbm", "actual", "b_sd"]], on=keys, suffixes=("_new", "_old"))
    return merged


def summarize(name, merged, stored_s):
    a = merged["f1s_lgbm_new"].to_numpy(dtype=np.float64)
    b = merged["f1s_lgbm_old"].to_numpy(dtype=np.float64)
    both = np.isfinite(a) & np.isfinite(b)
    diff = np.abs(a[both] - b[both])
    max_abs = float(diff.max()) if diff.size else np.nan
    spear = lookup.finite_spearman(a, b, 3)
    ra = rankdata(a[both], method="average")
    rb = rankdata(b[both], method="average")
    n_rank = int(np.sum(ra != rb))
    use = merged[np.isfinite(merged["actual"]) & np.isfinite(merged["f1s_lgbm_new"]) & np.isfinite(merged["b_sd"])]
    new_s = lookup.finite_spearman(use["f1s_lgbm_new"], use["actual"], 3)
    return {
        "dataset": name,
        "n_aligned": int(len(merged)),
        "n_finite_pred": int(both.sum()),
        "max_abs": max_abs,
        "spearman_new_vs_frozen": spear,
        "n_rank_changed": n_rank,
        "n_prereg_pairs": int(len(use)),
        "recomputed_S_F": new_s,
        "frozen_S_F": stored_s,
        "S_F_difference": new_s - stored_s if np.isfinite(new_s) else np.nan,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    REL.mkdir(parents=True, exist_ok=True)
    params = frozen_params()
    print("frozen params", params, flush=True)
    _train, x_f, y = training_matrix()
    genes, p_genes, adc = stage3.load_lists()
    rho, mean_pct, sd_pct, rate, gmean = rho_tables(genes, p_genes, adc)
    cache_r = OUT / "cache_rppa_features.tsv"
    cache_p = OUT / "cache_procan_features.tsv"
    if cache_r.exists() and cache_p.exists():
        print("loading cached features", flush=True)
        rppa = pd.read_csv(cache_r, sep="\t")
        pro = pd.read_csv(cache_p, sep="\t")
    else:
        print("building RPPA features", flush=True)
        rppa = rppa_features(rho, mean_pct, sd_pct, rate, gmean)
        print("building ProCan features", flush=True)
        pro = procan_features(genes, rho, mean_pct, sd_pct, rate, gmean)
        rppa.to_csv(cache_r, sep="\t", index=False)
        pro.to_csv(cache_p, sep="\t", index=False)
    confirm_r = pd.read_csv(RES / "stage4" / "confirm" / "rppa_pairs.tsv", sep="\t")
    confirm_p = pd.read_csv(RES / "stage4" / "confirm" / "procan_pairs.tsv", sep="\t")
    h = pd.read_csv(RES / "stage4" / "confirm" / "hypotheses.tsv", sep="\t").set_index("hypothesis")
    p2 = pd.read_csv(RES / "stage6" / "confirm" / "hypotheses.tsv", sep="\t").set_index("hypothesis")
    stored = {"RPPA": float(h.loc["H3", "S_F"]), "ProCan": float(h.loc["H4", "S_F"]), "APOLLO": float(p2.loc["P2", "S_F"])}
    attempts = [("attempt1", dict(params))]
    rows = []
    passed = None
    passed_fit = None
    passed_name = None
    frames = {}
    for name, pars in attempts:
        print("fit", name, flush=True)
        fit = fit_lgbm(pars, x_f, y)
        rppa_hat = rppa.copy()
        rppa_hat["f1s_lgbm"] = predict(fit, rppa)
        pro_hat = pro.copy()
        pro_hat["f1s_lgbm"] = predict(fit, pro)
        print("APOLLO features", flush=True)
        u_p2 = pd.read_csv(ROOT / "results" / "stage6" / "prep" / "U_P2_genes.tsv", sep="\t")["gene_symbol"].tolist()
        mat = pd.read_parquet(ROOT / "data" / "processed" / "apollo_luad" / "rna_log2tpm.parquet")
        samples = mat.columns.astype(str).tolist()
        table = stage3_reliability(genes, p_genes, u_p2, mat, samples, fit)
        merged_r = align(rppa_hat, confirm_r, ["xena_antibody", "cancer"])
        merged_p = align(pro_hat, confirm_p, ["protein_id", "lineage"])
        apollo_old = pd.read_csv(RES / "stage6" / "frozen" / "apollo_gene_pred.tsv", sep="\t")
        genes_used = pd.read_csv(RES / "stage6" / "confirm" / "genes.tsv", sep="\t")
        apollo_old = apollo_old.merge(genes_used[["gene_symbol", "actual", "in_U_P2"]], on="gene_symbol", how="left")
        merged_a = table.merge(apollo_old, on="gene_symbol", suffixes=("_new", "_old"))
        merged_a = merged_a.rename(columns={"f1s_lgbm_new": "f1s_lgbm_new", "f1s_lgbm_old": "f1s_lgbm_old"})
        if "b_sd_old" in merged_a.columns:
            merged_a["b_sd"] = merged_a["b_sd_old"]
        summaries = [
            summarize("RPPA", merged_r.rename(columns={"f1s_lgbm_x": "f1s_lgbm_new", "f1s_lgbm_y": "f1s_lgbm_old"}) if "f1s_lgbm_x" in merged_r.columns else merged_r, stored["RPPA"]),
            summarize("ProCan", merged_p, stored["ProCan"]),
            summarize("APOLLO", merged_a, stored["APOLLO"]),
        ]
        # align() already suffixes f1s_lgbm. summarize expects f1s_lgbm_new/old.
        # Fix RPPA/ProCan column names from the merge.
        for dataset, merged in (("RPPA", merged_r), ("ProCan", merged_p)):
            pass
        worst = max(rec["max_abs"] for rec in summaries)
        for rec in summaries:
            rec["attempt"] = name
            rec["passed"] = int(worst <= LIMIT)
            rows.append(rec)
            print(name, rec["dataset"], rec["max_abs"], rec["S_F_difference"], flush=True)
        if worst <= LIMIT:
            passed = summaries
            passed_fit = fit
            passed_name = name
            frames = {"rppa": rppa_hat, "pro": pro_hat, "apollo": table}
            break
        if name == "attempt1":
            second = dict(params)
            second["deterministic"] = True
            attempts.append(("attempt2", second))
    pd.DataFrame(rows).to_csv(OUT / "refit_compare.tsv", sep="\t", index=False)
    if passed_fit is None:
        print("export skipped", flush=True)
        print("stage9 refit done", flush=True)
        return
    text_path = REL / "F1s_lightgbm.txt"
    passed_fit["lgbm"].booster_.save_model(str(text_path))
    booster = Booster(model_file=str(text_path))
    reloaded = []
    for frame in frames.values():
        x = frame[F1S_FEATURES].to_numpy(dtype=np.float64)
        filled = np.array(x, dtype=np.float64, copy=True)
        med = np.nanmedian(passed_fit["x_tr"], axis=0)
        med = np.where(np.isfinite(med), med, 0.0)
        inds = np.where(~np.isfinite(filled))
        if inds[0].size:
            filled[inds] = np.take(med, inds[1])
        a = passed_fit["lgbm"].predict(filled)
        b = booster.predict(filled)
        reloaded.append(float(np.max(np.abs(a - b))))
    pd.DataFrame([{"attempt": passed_name, "reload_max_abs": max(reloaded), "text_model": str(text_path.relative_to(ROOT))}]).to_csv(OUT / "refit_reload.tsv", sep="\t", index=False)
    print("reload max abs", max(reloaded), flush=True)
    print("stage9 refit done", flush=True)


def stage3_reliability(genes, p_genes, targets, mat, samples, f1s_fit):
    """Feature rows and predictions for APOLLO, without writing the frozen table."""
    rho = stage3.defined_rho(p_genes)
    adc = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    adc_loco = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    pset = set(p_genes)
    for rec in adc_loco.itertuples(index=False):
        if rec.gene_symbol in pset or rec.gene_symbol not in set(adc):
            continue
        if int(rec.n) >= 30 and np.isfinite(rec.s0):
            rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
    rho = stage3.forecast_rho(rho, targets)
    _stats, mean_pct, sd_pct, _ = stage3.cptac_percentile_tables(genes)
    rates = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    rate = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
    gmean = stage3.gmean_all(rho, p_genes)
    mu, sd = stage3.gene_mean_sd(mat.reindex(index=genes), samples)
    mu_pct = stage3.percentile_map(mu)
    sd_pct_ext = stage3.percentile_map(sd)
    rows = []
    for gene in targets:
        ext_mu_pct = mu_pct.get(gene, (np.nan, np.nan))[1]
        hold_pct = sd_pct_ext.get(gene, (np.nan, np.nan))[1]
        rows.append(stage3.f1s_feature_row(gene, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, hold_pct) | {"gene_symbol": gene})
    frame = pd.DataFrame(rows)
    frame["f1s_lgbm"] = predict(f1s_fit, frame)
    return frame


if __name__ == "__main__":
    main()
