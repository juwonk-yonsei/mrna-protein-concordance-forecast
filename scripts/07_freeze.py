"""Section 4. Freeze RNA-proxy forecasts. Does not read locked values."""
import csv
import gzip
import hashlib
import os
import sys
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import rankdata
from sklearn.linear_model import Ridge

from scripts.b2_lib import (
    CANCERS,
    DDOF,
    FEATURE_NAMES,
    MIN_N,
    PLAN,
    RAW,
    RES,
    ROOT,
    conformal_half_width,
    spearman,
    zscore_cols,
)

PROC = ROOT / "data" / "processed"


def fill_median(train_x, apply_x):
    med = np.nanmedian(train_x, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)

    def _fill(x):
        out = np.array(x, dtype=np.float64, copy=True)
        inds = np.where(~np.isfinite(out))
        if inds[0].size:
            out[inds] = np.take(med, inds[1])
        return out

    return _fill(train_x), _fill(apply_x), med


def standardize(train_x, apply_x):
    mu = train_x.mean(axis=0)
    sd = train_x.std(axis=0, ddof=0)
    bad = (sd == 0) | ~np.isfinite(sd)
    sd = np.where(bad, 1.0, sd)

    def _z(x):
        z = (x - mu) / sd
        z[:, bad] = 0.0
        return z

    return _z(train_x), _z(apply_x)


def lgbm_params():
    cfg = dict(PLAN["reliability"]["lightgbm"])
    cfg["verbosity"] = -1
    return cfg


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_lists():
    genes = [g for g in (RES / "G_RNA_genes.txt").read_text().split() if g]
    p_genes = [g for g in (RES / "P_genes.txt").read_text().split() if g]
    adc = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t", dtype=str)["gene_symbol"].tolist()
    return genes, p_genes, adc


def load_defined_rho(p_genes, adc):
    s0 = pd.read_csv(RES / "s0_pairs.tsv", sep="\t")
    use = s0[(s0["n"] >= MIN_N) & np.isfinite(s0["s0"])]
    rho = {(r.cancer, r.gene_symbol): float(r.s0) for r in use.itertuples(index=False)}
    adc_path = RES / "adc_loco.tsv"
    if adc_path.exists():
        extra = pd.read_csv(adc_path, sep="\t")
        pset = set(p_genes)
        for r in extra.itertuples(index=False):
            if r.gene_symbol in pset:
                continue
            if r.gene_symbol not in set(adc):
                continue
            if int(r.n) >= MIN_N and np.isfinite(r.s0):
                rho[(r.heldout_cancer, r.gene_symbol)] = float(r.s0)
    return rho


def extra_cptac_rho(genes, rho):
    """S0 for forecast genes outside the stage-1 P table. Protein values are CPTAC, not locked."""
    need = [g for g in genes if g not in {k[1] for k in rho}]
    if not need:
        return rho, 0
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "tumor"]
    added = 0
    for cancer in CANCERS:
        prot = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__protein.parquet")
        rna = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__rna.parquet")
        present = [g for g in need if g in prot.index and g in rna.index]
        if not present:
            continue
        sids = paired.loc[paired["cancer"] == cancer, "sample_id"].tolist()
        sids = [s for s in sids if s in rna.columns and s in prot.columns]
        x = zscore_cols(rna.reindex(index=present, columns=sids).to_numpy(dtype=np.float64).T)
        y = zscore_cols(prot.reindex(index=present, columns=sids).to_numpy(dtype=np.float64).T)
        for j, gene in enumerate(present):
            n = int(np.isfinite(y[:, j]).sum())
            if n < MIN_N:
                continue
            val = spearman(x[:, j], y[:, j])
            if np.isfinite(val):
                rho[(cancer, gene)] = float(val)
                added += 1
    return rho, added


def load_rates(extra_genes):
    rates = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    out = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
    missing = [g for g in extra_genes if not any((c, g) in out for c in CANCERS)]
    if not missing:
        return out
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "tumor"]
    for cancer in CANCERS:
        prot = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__protein.parquet")
        sids = [s for s in paired.loc[paired["cancer"] == cancer, "sample_id"] if s in prot.columns]
        present = [g for g in missing if g in prot.index]
        if not present or not sids:
            for gene in missing:
                out.setdefault((cancer, gene), 0.0)
            continue
        mat = prot.reindex(index=present, columns=sids).to_numpy(dtype=np.float64)
        finite = np.isfinite(mat).sum(axis=1) / float(len(sids))
        for gene, rate in zip(present, finite):
            out[(cancer, gene)] = float(rate)
        for gene in missing:
            if gene not in prot.index:
                out[(cancer, gene)] = 0.0
    return out


def cptac_rna_tables(genes):
    stats = pd.read_csv(RES / "cptac_rna_gene_stats.tsv", sep="\t")
    stats = stats[stats["gene_symbol"].isin(set(genes))]
    rna_mean = {(r.cancer, r.gene_symbol): float(r.rna_mean) for r in stats.itertuples(index=False)}
    pct = {}
    for cancer, sub in stats.groupby("cancer"):
        sd = sub["rna_sd"].to_numpy(dtype=np.float64)
        names = sub["gene_symbol"].tolist()
        finite = np.isfinite(sd)
        rank = np.full(len(sd), np.nan)
        if finite.any():
            rank[finite] = 100.0 * rankdata(sd[finite], method="average") / float(finite.sum())
        for gene, s, r in zip(names, sd, rank):
            pct[(cancer, gene)] = (float(s) if np.isfinite(s) else np.nan, float(r) if np.isfinite(r) else np.nan)
    return rna_mean, pct


def feature_row(gene, source, rho, rna_mean, rate, global_mean, ext_mean, ext_sd, ext_pct):
    vals = [rho[(c, gene)] for c in source if (c, gene) in rho]
    if vals:
        rho_mean = float(np.mean(vals))
        rho_sd = float(np.std(vals, ddof=DDOF)) if len(vals) >= 2 else 0.0
        rho_min = float(np.min(vals))
        rho_max = float(np.max(vals))
        n_def = len(vals)
    else:
        rho_mean = global_mean
        rho_sd = 0.0
        rho_min = global_mean
        rho_max = global_mean
        n_def = 0
    means = [rna_mean.get((c, gene), np.nan) for c in source]
    means = [v for v in means if np.isfinite(v)]
    rates = [rate.get((c, gene), 0.0) for c in source]
    return {
        "rho_mean": rho_mean,
        "rho_sd": rho_sd,
        "rho_min": rho_min,
        "rho_max": rho_max,
        "n_defined_cancers": n_def,
        "holdout_rna_mean": ext_mean,
        "holdout_rna_sd": ext_sd,
        "holdout_rna_sd_percentile": ext_pct,
        "train_rna_mean": float(np.mean(means)) if means else np.nan,
        "train_protein_rate_mean": float(np.mean(rates)) if rates else np.nan,
    }


def global_mean(rho, p_genes):
    pset = set(p_genes)
    vals = [v for (c, g), v in rho.items() if c in CANCERS and g in pset]
    if not vals:
        raise RuntimeError("no defined CPTAC rho for the global mean")
    return float(np.mean(vals))


def fit_frozen(train):
    train = train[np.isfinite(train["actual"])].copy()
    alphas = [float(a) for a in PLAN["reliability"]["ridge_alphas"]]
    cancers = list(dict.fromkeys(train["query_cancer"]))
    x_all = train[FEATURE_NAMES].to_numpy(dtype=np.float64)
    y_all = train["actual"].to_numpy(dtype=np.float64)
    groups = [train["query_cancer"].to_numpy() == c for c in cancers]
    maes = []
    for alpha in alphas:
        per = []
        for mask in groups:
            tr = ~mask & np.isfinite(y_all)
            te = mask & np.isfinite(y_all)
            if tr.sum() < 1 or te.sum() < 1:
                continue
            xtr, xte, _ = fill_median(x_all[tr], x_all[te])
            ztr, zte = standardize(xtr, xte)
            model = Ridge(alpha=alpha, fit_intercept=True)
            model.fit(ztr, y_all[tr])
            per.append(float(np.mean(np.abs(model.predict(zte) - y_all[te]))))
        maes.append(float(np.mean(per)) if per else np.nan)
    maes = np.asarray(maes, dtype=np.float64)
    if not np.isfinite(maes).any():
        chosen = float(PLAN["reliability"]["ridge_fallback_alpha"])
    else:
        best = np.nanmin(maes)
        chosen = alphas[0]
        for alpha, val in zip(alphas, maes):
            if np.isfinite(val) and val == best:
                chosen = alpha
                break
    # Section 3.2 holdout predictions are the pooled out-of-fold scores.
    widths = {}
    for name, col in (("f0", "f0"), ("ridge", "f1_ridge"), ("lgbm", "f1_lgbm")):
        resid = np.abs(train[col].to_numpy(dtype=np.float64) - y_all)
        widths[(name, 80)] = conformal_half_width(resid, 80)
        widths[(name, 90)] = conformal_half_width(resid, 90)
    x_tr, _, med = fill_median(x_all, x_all[:1])
    z_tr, _ = standardize(x_tr, x_tr[:1])
    ridge = Ridge(alpha=chosen, fit_intercept=True)
    ridge.fit(z_tr, y_all)
    lgbm = LGBMRegressor(**lgbm_params())
    lgbm.fit(x_tr, y_all)
    return {
        "alpha": chosen,
        "maes": maes,
        "alphas": alphas,
        "widths": widths,
        "ridge": ridge,
        "lgbm": lgbm,
        "median": med,
        "mu": x_tr.mean(axis=0),
        "sd": None,
        "x_tr": x_tr,
        "n_train": int(len(train)),
    }


def apply_models(frame, fit):
    x = frame[FEATURE_NAMES].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(fit["median"], inds[1])
    mu = fit["x_tr"].mean(axis=0)
    sd = fit["x_tr"].std(axis=0, ddof=0)
    bad = (sd == 0) | ~np.isfinite(sd)
    sd = np.where(bad, 1.0, sd)
    z = (filled - mu) / sd
    z[:, bad] = 0.0
    out = frame.copy()
    out["f0"] = out["rho_mean"]
    out["f1_ridge"] = fit["ridge"].predict(z)
    out["f1_lgbm"] = fit["lgbm"].predict(filled)
    for name, col in (("f0", "f0"), ("ridge", "f1_ridge"), ("lgbm", "f1_lgbm")):
        for level, tag in ((80, "80"), (90, "90")):
            half = fit["widths"][(name, level)]
            pred = out[col].to_numpy(dtype=np.float64)
            out[f"{col}_low_{tag}"] = pred - half
            out[f"{col}_high_{tag}"] = pred + half
            out[f"{col}_half_{tag}"] = half
    out["ridge_alpha"] = fit["alpha"]
    return out


def percentile_map(sd_by_gene):
    genes = list(sd_by_gene)
    sd = np.asarray([sd_by_gene[g] for g in genes], dtype=np.float64)
    finite = np.isfinite(sd)
    pct = np.full(len(sd), np.nan)
    if finite.any():
        pct[finite] = 100.0 * rankdata(sd[finite], method="average") / float(finite.sum())
    return {g: (float(s) if np.isfinite(s) else np.nan, float(p) if np.isfinite(p) else np.nan) for g, s, p in zip(genes, sd, pct)}


def entrez_map(g_rna):
    hgnc = pd.read_csv(
        RAW / "hgnc" / "hgnc_complete_set.txt",
        sep="\t",
        dtype=str,
        usecols=["symbol", "locus_group", "entrez_id"],
    )
    pc = hgnc[hgnc["locus_group"] == "protein-coding gene"].copy()
    pc = pc[pc["entrez_id"].notna() & (pc["entrez_id"] != "")]
    n_entrez_multi_symbol = int((pc.groupby("entrez_id")["symbol"].nunique() > 1).sum())
    pc = pc.groupby("entrez_id").filter(lambda d: d["symbol"].nunique() == 1)
    pc = pc.drop_duplicates(["entrez_id", "symbol"])
    n_symbol_multi_entrez = int((pc.groupby("symbol")["entrez_id"].nunique() > 1).sum())
    pc = pc.groupby("symbol").filter(lambda d: d["entrez_id"].nunique() == 1)
    pc = pc[pc["symbol"].isin(set(g_rna))]
    mapping = dict(zip(pc["entrez_id"], pc["symbol"]))
    summary = {
        "n_entrez_mapping_to_multiple_symbols_dropped": n_entrez_multi_symbol,
        "n_symbols_mapping_to_multiple_entrez_dropped": n_symbol_multi_entrez,
        "n_g_rna_symbols_with_one_entrez": len(mapping),
    }
    return mapping, summary


def tcga_rna_features(g_rna, diseases):
    mapping, summary = entrez_map(g_rna)
    rna_path = RAW / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    ph = pd.read_csv(RAW / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", dtype=str)
    with gzip.open(rna_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    from collections import Counter

    counts = Counter(header[1:])
    rna_unique = {s for s, n in counts.items() if n == 1}
    dup_ids = set(pd.read_csv(RES / "rppa_duplicate_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    overlap = pd.read_csv(RES / "tcga_rppa_rna_overlap_by_disease.tsv", sep="\t")
    keep_dis = set(overlap.loc[overlap["n_ge_30"].astype(str).isin(["True", "true", "1"]), "_primary_disease"])
    if set(diseases) != keep_dis:
        raise RuntimeError("disease list does not match the section-2 overlap table")
    # Header only. The value rows of the locked RPPA matrix are not read.
    rppa_path = RAW / "tcga" / "TCGA-RPPA-pancan-clean.xena.gz"
    with gzip.open(rppa_path, "rt", encoding="utf-8", errors="replace") as handle:
        rppa_header = handle.readline().rstrip("\n").split("\t")
    rppa_ids = {s for s in rppa_header[1:] if s not in dup_ids}
    both = rppa_ids & rna_unique
    samples_by = {}
    wanted = set()
    for disease in diseases:
        ids = [s for s in ph.loc[ph["_primary_disease"] == disease, "sample"] if s in both]
        if len(ids) != int(overlap.loc[overlap["_primary_disease"] == disease, "n_samples"].iloc[0]):
            raise RuntimeError(f"sample count mismatch for {disease}: {len(ids)}")
        samples_by[disease] = ids
        wanted.update(ids)
    usecols = ["sample"] + sorted(wanted)
    df = pd.read_csv(rna_path, sep="\t", usecols=usecols, dtype={"sample": str}, low_memory=False)
    gset = set(g_rna)
    # Row ids are gene symbols, except a leading block of numeric Entrez ids.
    symbol = []
    n_exact = 0
    n_entrez = 0
    n_unmatched = 0
    for raw_id in df["sample"].astype(str):
        if raw_id in gset:
            symbol.append(raw_id)
            n_exact += 1
        elif raw_id.isdigit() and raw_id in mapping and mapping[raw_id] in gset:
            symbol.append(mapping[raw_id])
            n_entrez += 1
        else:
            symbol.append(None)
            n_unmatched += 1
    df["symbol"] = symbol
    mapped = df[df["symbol"].notna()].copy()
    n_rows = int(df.shape[0])
    n_mapped_rows = int(mapped.shape[0])
    dup_symbols = set(mapped.loc[mapped["symbol"].duplicated(keep=False), "symbol"])
    n_ambiguous_symbol_rows = int(mapped["symbol"].isin(dup_symbols).sum())
    if dup_symbols:
        mapped = mapped[~mapped["symbol"].isin(dup_symbols)]
    mapped = mapped.set_index("symbol")
    mapped = mapped.drop(columns=["sample"])
    features = {}
    for disease, ids in samples_by.items():
        sub = mapped.reindex(columns=ids)
        mat = sub.to_numpy(dtype=np.float64)
        with np.errstate(all="ignore"):
            mu = np.nanmean(mat, axis=1)
            sd = np.nanstd(mat, axis=1, ddof=DDOF)
        sd_map = {g: float(s) if np.isfinite(s) else np.nan for g, s in zip(sub.index, sd)}
        pct = percentile_map(sd_map)
        features[disease] = {
            g: (float(m) if np.isfinite(m) else np.nan, pct[g][0], pct[g][1])
            for g, m in zip(sub.index, mu)
        }
    summary.update(
        {
            "n_rna_rows": n_rows,
            "n_exact_symbol_in_g_rna": n_exact,
            "n_numeric_entrez_mapped_to_g_rna": n_entrez,
            "n_unmatched_row_ids": n_unmatched,
            "n_rna_rows_mapped_to_g_rna": n_mapped_rows,
            "n_duplicate_symbol_rows_dropped": n_ambiguous_symbol_rows,
            "duplicate_symbols_dropped": ";".join(sorted(dup_symbols)) if dup_symbols else "NA",
            "n_samples_used": len(wanted),
            "n_diseases": len(diseases),
        }
    )
    pd.DataFrame([{"item": k, "value": v} for k, v in summary.items()]).to_csv(
        RES / "tcga_entrez_map_summary.tsv", sep="\t", index=False
    )
    return features, samples_by


def procan_rna_features(g_rna, lineages):
    gset = set(g_rna)
    mapping = pd.read_csv(RAW / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t", dtype=str)
    sample_ids = set(pd.read_csv(RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    mapping = mapping[mapping["Project_Identifier"].isin(sample_ids)]
    sidm_lineage = dict(zip(mapping["SIDM"], mapping["Tissue_type"]))
    wanted = {sidm for sidm, lin in sidm_lineage.items() if lin in set(lineages)}
    path = RAW / "cmp" / "rnaseq_merged" / "rnaseq_merged_20260323.csv"
    stored = {}
    n_dup = 0
    n_sanger_rows = 0
    n_kept_rows = 0
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
            n_sanger_rows += 1
            model = parts[i_model]
            if model not in wanted:
                continue
            gene = parts[i_gene]
            if gene not in gset:
                continue
            key = (model, gene)
            if key in stored:
                n_dup += 1
                continue
            try:
                tpm = float(parts[i_tpm])
            except ValueError:
                continue
            if not np.isfinite(tpm):
                continue
            stored[key] = np.log2(tpm + 1.0)
            n_kept_rows += 1
    by_lineage = {lin: {} for lin in lineages}
    for (model, gene), val in stored.items():
        lin = sidm_lineage[model]
        bucket = by_lineage[lin].setdefault(gene, [])
        bucket.append(val)
    features = {}
    for lin, genes in by_lineage.items():
        sd_map = {}
        mu_map = {}
        n_map = {}
        for gene, vals in genes.items():
            arr = np.asarray(vals, dtype=np.float64)
            mu_map[gene] = float(arr.mean())
            sd_map[gene] = float(arr.std(ddof=DDOF))
            n_map[gene] = int(arr.size)
        pct = percentile_map(sd_map)
        features[lin] = {
            g: (mu_map[g], pct[g][0], pct[g][1], n_map[g]) for g in mu_map
        }
    meta = {
        "n_sanger_rows": n_sanger_rows,
        "n_kept_model_gene": n_kept_rows,
        "n_duplicate_model_gene_skipped": n_dup,
        "n_models_wanted": len(wanted),
    }
    pd.DataFrame([{"item": k, "value": v} for k, v in meta.items()]).to_csv(
        RES / "procan_rna_agg_summary.tsv", sep="\t", index=False
    )
    return features


def build_rows(records, fit, rho, rna_mean, rate, p_genes, gmean):
    rows = []
    for rec in records:
        feat = feature_row(
            rec["gene_symbol"],
            CANCERS,
            rho,
            rna_mean,
            rate,
            gmean,
            rec["ext_mean"],
            rec["ext_sd"],
            rec["ext_pct"],
        )
        feat.update(rec)
        rows.append(feat)
    frame = pd.DataFrame(rows)
    return apply_models(frame, fit)


def prediction_columns(prefix_rows):
    cols = []
    for name in ("f0", "f1_ridge", "f1_lgbm"):
        cols.extend([name, f"{name}_low_80", f"{name}_high_80", f"{name}_low_90", f"{name}_high_90", f"{name}_half_80", f"{name}_half_90"])
    return ["ridge_alpha"] + cols


def main():
    genes, p_genes, adc = load_lists()
    preds = pd.read_csv(RES / "reliability_predictions.tsv", sep="\t")
    train = preds[preds["in_P"].astype(bool) & np.isfinite(preds["actual"])].copy()
    if train["query_cancer"].nunique() != 10:
        raise RuntimeError("freeze training rows are not one block per cancer")
    fit = fit_frozen(train)
    rho = load_defined_rho(p_genes, adc)
    antibodies = pd.read_csv(RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    antibodies = antibodies[antibodies["kept_for_forecast"].astype(str).isin(["True", "true"])]
    proteins = pd.read_csv(RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    proteins = proteins[(proteins["symbol"] != "NA") & (proteins["symbol"] != "")]
    proteins = proteins[proteins["symbol"].isin(set(genes))]
    forecast_genes = sorted(set(antibodies["gene_symbol"]) | set(proteins["symbol"]))
    rho, n_added = extra_cptac_rho(forecast_genes, rho)
    print(f"extra cptac rho pairs {n_added}", flush=True)
    rna_mean, _ = cptac_rna_tables(genes)
    rate = load_rates([g for g in forecast_genes if g not in set(p_genes)])
    gmean = global_mean(rho, p_genes)

    overlap = pd.read_csv(RES / "tcga_rppa_rna_overlap_by_disease.tsv", sep="\t")
    diseases = overlap.loc[overlap["n_ge_30"], "_primary_disease"].tolist()
    print(f"tcga rna diseases {len(diseases)}", flush=True)
    tcga_feat, samples_by = tcga_rna_features(genes, diseases)
    rppa_records = []
    for disease in diseases:
        block = tcga_feat[disease]
        n = len(samples_by[disease])
        for row in antibodies.itertuples(index=False):
            ext = block.get(row.gene_symbol, (np.nan, np.nan, np.nan))
            rppa_records.append(
                {
                    "xena_antibody": row.xena_antibody,
                    "gene_symbol": row.gene_symbol,
                    "cancer": disease,
                    "n": n,
                    "ext_mean": ext[0],
                    "ext_sd": ext[1],
                    "ext_pct": ext[2],
                    "is_adc": row.gene_symbol in set(adc),
                }
            )
    rppa = build_rows(rppa_records, fit, rho, rna_mean, rate, p_genes, gmean)
    rppa_path = RES / "frozen_forecasts_rppa.tsv"
    rppa_cols = ["xena_antibody", "gene_symbol", "cancer", "n", "is_adc"] + prediction_columns(rppa)
    # rho_mean was overwritten into f0; keep the feature columns out of the frozen file.
    rppa[rppa_cols].to_csv(rppa_path, sep="\t", index=False)

    lin_table = pd.read_csv(RES / "procan_sanger_by_lineage.tsv", sep="\t")
    lineages = lin_table.loc[lin_table["n_ge_30"], "Tissue_type"].tolist()
    n_by_lin = dict(zip(lin_table["Tissue_type"], lin_table["n_cell_lines"].astype(int)))
    print(f"procan lineages {len(lineages)}", flush=True)
    procan_feat = procan_rna_features(genes, lineages)
    procan_records = []
    for lin in lineages:
        block = procan_feat[lin]
        for row in proteins.itertuples(index=False):
            ext = block.get(row.symbol)
            if ext is None:
                mean, sd, pct, n_gene = np.nan, np.nan, np.nan, 0
            else:
                mean, sd, pct, n_gene = ext
            procan_records.append(
                {
                    "protein_id": row.protein_id,
                    "gene_symbol": row.symbol,
                    "lineage": lin,
                    "n": int(n_by_lin[lin]),
                    "n_rna_gene": int(n_gene),
                    "ext_mean": mean,
                    "ext_sd": sd,
                    "ext_pct": pct,
                    "is_adc": row.symbol in set(adc),
                }
            )
    procan = build_rows(procan_records, fit, rho, rna_mean, rate, p_genes, gmean)
    procan_path = RES / "frozen_forecasts_procan.tsv"
    procan_cols = ["protein_id", "gene_symbol", "lineage", "n", "n_rna_gene", "is_adc"] + prediction_columns(procan)
    procan[procan_cols].to_csv(procan_path, sep="\t", index=False)

    hash_rows = []
    for key, path in (("rppa", rppa_path), ("procan", procan_path)):
        hash_rows.append(
            {
                "key": key,
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256_file(path),
                "size": path.stat().st_size,
            }
        )
    hash_path = ROOT / "config" / "frozen_hashes.tsv"
    pd.DataFrame(hash_rows).to_csv(hash_path, sep="\t", index=False)
    pd.DataFrame(
        [
            {"item": "ridge_alpha", "value": fit["alpha"]},
            {"item": "n_train_pairs", "value": fit["n_train"]},
            {"item": "q80_f0", "value": fit["widths"][("f0", 80)]},
            {"item": "q90_f0", "value": fit["widths"][("f0", 90)]},
            {"item": "q80_ridge", "value": fit["widths"][("ridge", 80)]},
            {"item": "q90_ridge", "value": fit["widths"][("ridge", 90)]},
            {"item": "q80_lgbm", "value": fit["widths"][("lgbm", 80)]},
            {"item": "q90_lgbm", "value": fit["widths"][("lgbm", 90)]},
            {"item": "n_rppa_rows", "value": len(rppa)},
            {"item": "n_procan_rows", "value": len(procan)},
            {"item": "n_extra_cptac_rho", "value": n_added},
        ]
    ).to_csv(RES / "frozen_model_fit.tsv", sep="\t", index=False)
    print("freeze done", fit["alpha"], len(rppa), len(procan), flush=True)


if __name__ == "__main__":
    main()
