"""Stage 3 development checks and freezes. Does not read locked values."""
import gzip
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
from sklearn.metrics import roc_auc_score

from scripts.b2_lib import (
    ALPHAS,
    CANCERS,
    DDOF,
    FEATURE_NAMES,
    K_EXTRA,
    LAMBDAS,
    MIN_N,
    PROC,
    RAW,
    RES,
    ROOT,
    conformal_half_width,
    fit_slope,
    residuals,
    select_residual_features,
    spearman,
    zscore_cols,
)
from scripts.b3_lib import (
    F1S_FEATURES,
    PLAN,
    append_frozen_hashes,
    centered_spearman,
    fit_ar_coefficients,
    map_row_ids,
    percentile_map,
    predict_ar_matrix,
    quantiles,
    self_check_map,
    stage2_map_row_ids,
    within_gene_summary,
)
from scripts.common import sha256_file

warnings_filter = True


def quiet_warnings():
    import warnings

    warnings.filterwarnings("ignore", message="Mean of empty slice")
    warnings.filterwarnings("ignore", message="Degrees of freedom")


def load_lists():
    genes = [g for g in (RES / "G_RNA_genes.txt").read_text().split() if g]
    p_genes = [g for g in (RES / "P_genes.txt").read_text().split() if g]
    adc = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t", dtype=str)["gene_symbol"].tolist()
    return genes, p_genes, adc


def tcga_row_ids():
    path = RAW / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    ids = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        handle.readline()
        for line in handle:
            ids.append(line.split("\t", 1)[0])
    return ids


def write_id_map(genes):
    ids = tcga_row_ids()
    mapped = map_row_ids(ids, genes)
    mapped.to_csv(RES / "tcga_rna_id_map_b3.tsv", sep="\t", index=False)
    counts = mapped["rule"].value_counts()
    summary = [{"item": f"n_{rule}", "value": int(counts.get(rule, 0))} for rule in ["exact", "entrez", "prev_symbol", "unmatched"]]
    summary.append({"item": "n_rows", "value": int(len(mapped))})
    summary.append({"item": "n_symbols", "value": int(mapped.loc[mapped.rule != "unmatched", "gene_symbol"].nunique())})
    pd.DataFrame(summary).to_csv(RES / "tcga_rna_id_summary_b3.tsv", sep="\t", index=False)
    antibodies = pd.read_csv(RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    by_symbol = {}
    for row in mapped.itertuples(index=False):
        if row.rule == "unmatched":
            continue
        by_symbol.setdefault(row.gene_symbol, set()).add(row.rule)
    order = ["exact", "entrez", "prev_symbol"]
    routes = []
    for gene in antibodies["gene_symbol"]:
        rules = by_symbol.get(gene, set())
        route = next((name for name in order if name in rules), "none")
        routes.append({"gene_symbol": gene, "route": route})
    route_df = pd.DataFrame(routes)
    route_df.to_csv(RES / "rppa134_id_route.tsv", sep="\t", index=False)
    stage2 = stage2_map_row_ids(ids, genes)
    stage2_symbols = set(stage2.values())
    missing = sorted(set(antibodies["gene_symbol"]) - stage2_symbols)
    pd.DataFrame({"gene_symbol": missing}).to_csv(RES / "stage2_rppa_genes_without_rna_row.tsv", sep="\t", index=False)
    print("id map", {r["item"]: r["value"] for r in summary}, "stage2 missing genes", len(missing), flush=True)
    return mapped


def defined_rho(p_genes):
    s0 = pd.read_csv(RES / "s0_pairs.tsv", sep="\t")
    use = s0[(s0["n"] >= MIN_N) & np.isfinite(s0["s0"])]
    rho = {(r.cancer, r.gene_symbol): float(r.s0) for r in use.itertuples(index=False)}
    return rho


def cptac_percentile_tables(genes):
    stats = pd.read_csv(RES / "cptac_rna_gene_stats.tsv", sep="\t")
    stats = stats[stats["gene_symbol"].isin(set(genes))]
    mean_pct = {}
    sd_pct = {}
    mean_val = {}
    for cancer, sub in stats.groupby("cancer"):
        genes_c = sub["gene_symbol"].tolist()
        mu = dict(zip(genes_c, sub["rna_mean"]))
        sd = dict(zip(genes_c, sub["rna_sd"]))
        mu_map = percentile_map(mu)
        sd_map = percentile_map(sd)
        for gene in genes_c:
            mean_pct[(cancer, gene)] = mu_map[gene][1]
            sd_pct[(cancer, gene)] = sd_map[gene][1]
            mean_val[(cancer, gene)] = mu_map[gene][0]
    return stats, mean_pct, sd_pct, mean_val


def f1s_feature_row(gene, source, rho, mean_pct, sd_pct, rate, global_mean, ext_mean_pct, ext_sd_pct):
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
    pcts = [mean_pct.get((c, gene), np.nan) for c in source]
    pcts = [v for v in pcts if np.isfinite(v)]
    rates = [rate.get((c, gene), 0.0) for c in source]
    return {
        "rho_mean": rho_mean,
        "rho_sd": rho_sd,
        "rho_min": rho_min,
        "rho_max": rho_max,
        "n_defined_cancers": n_def,
        "holdout_rna_mean_percentile": ext_mean_pct,
        "holdout_rna_sd_percentile": ext_sd_pct,
        "train_rna_mean_percentile_mean": float(np.mean(pcts)) if pcts else np.nan,
        "train_protein_rate_mean": float(np.mean(rates)) if rates else np.nan,
    }


def global_mean(source, rho, p_genes):
    pset = set(p_genes)
    vals = [rho[(c, g)] for c in source for g in p_genes if (c, g) in rho and g in pset]
    return float(np.mean(vals)) if vals else np.nan


def design_f1s(genes, source, query, rho, mean_pct, sd_pct, rate, p_genes):
    gmean = global_mean(source, rho, p_genes)
    rows = []
    for gene in genes:
        rec = f1s_feature_row(
            gene,
            source,
            rho,
            mean_pct,
            sd_pct,
            rate,
            gmean,
            mean_pct.get((query, gene), np.nan),
            sd_pct.get((query, gene), np.nan),
        )
        rec["gene_symbol"] = gene
        rec["query_cancer"] = query
        rec["actual"] = rho.get((query, gene), np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


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
    sd = train_x.std(axis=0, ddof=0)
    bad = (sd == 0) | ~np.isfinite(sd)
    sd = np.where(bad, 1.0, sd)
    mu = train_x.mean(axis=0)

    def _z(x):
        z = (x - mu) / sd
        z[:, bad] = 0.0
        return z

    return _z(train_x), _z(apply_x)


def lgbm_params():
    cfg = dict(B2_LGBM)
    cfg["verbosity"] = -1
    return cfg


def ridge_alpha_and_oof(frame, features):
    alphas = [float(a) for a in B2["reliability"]["ridge_alphas"]]
    cancers = list(dict.fromkeys(frame["query_cancer"]))
    x_all = frame[features].to_numpy(dtype=np.float64)
    y_all = frame["actual"].to_numpy(dtype=np.float64)
    groups = [frame["query_cancer"].to_numpy() == c for c in cancers]
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
        chosen = float(B2["reliability"]["ridge_fallback_alpha"])
    else:
        best = np.nanmin(maes)
        chosen = next(a for a, v in zip(alphas, maes) if np.isfinite(v) and v == best)
    oof = np.full(len(frame), np.nan)
    for mask in groups:
        tr = ~mask & np.isfinite(y_all)
        te = mask & np.isfinite(y_all)
        if tr.sum() < 1 or te.sum() < 1:
            continue
        xtr, xte, _ = fill_median(x_all[tr], x_all[te])
        ztr, zte = standardize(xtr, xte)
        model = Ridge(alpha=chosen, fit_intercept=True)
        model.fit(ztr, y_all[tr])
        oof[np.flatnonzero(te)] = model.predict(zte)
    return chosen, oof


def lgbm_oof(frame, features):
    cancers = list(dict.fromkeys(frame["query_cancer"]))
    x_all = frame[features].to_numpy(dtype=np.float64)
    y_all = frame["actual"].to_numpy(dtype=np.float64)
    oof = np.full(len(frame), np.nan)
    for cancer in cancers:
        mask = frame["query_cancer"].to_numpy() == cancer
        tr = ~mask & np.isfinite(y_all)
        te = mask & np.isfinite(y_all)
        if tr.sum() < 1 or te.sum() < 1:
            continue
        xtr, xte, _ = fill_median(x_all[tr], x_all[te])
        model = LGBMRegressor(**lgbm_params())
        model.fit(xtr, y_all[tr])
        oof[np.flatnonzero(te)] = model.predict(xte)
    return oof


def fit_apply(train, apply, features):
    train = train[np.isfinite(train["actual"])].copy()
    alpha, oof_ridge = ridge_alpha_and_oof(train, features)
    oof_lgb = lgbm_oof(train, features)
    train = train.copy()
    train["oof_ridge"] = oof_ridge
    train["oof_lgbm"] = oof_lgb
    train["oof_f0"] = train["rho_mean"]
    x_tr = train[features].to_numpy(dtype=np.float64)
    y_tr = train["actual"].to_numpy(dtype=np.float64)
    x_ap = apply[features].to_numpy(dtype=np.float64)
    x_tr_f, x_ap_f, _ = fill_median(x_tr, x_ap)
    z_tr, z_ap = standardize(x_tr_f, x_ap_f)
    ridge = Ridge(alpha=alpha, fit_intercept=True)
    ridge.fit(z_tr, y_tr)
    lgbm = LGBMRegressor(**lgbm_params())
    lgbm.fit(x_tr_f, y_tr)
    pred = apply.copy()
    pred["f0"] = pred["rho_mean"]
    pred["f1s_ridge"] = ridge.predict(z_ap)
    pred["f1s_lgbm"] = lgbm.predict(x_ap_f)
    pred["ridge_alpha"] = alpha
    widths = {}
    for name, col in (("f0", "oof_f0"), ("ridge", "oof_ridge"), ("lgbm", "oof_lgbm")):
        resid = np.abs(train[col].to_numpy(dtype=np.float64) - y_tr)
        widths[(name, 80)] = conformal_half_width(resid, 80)
        widths[(name, 90)] = conformal_half_width(resid, 90)
        pred[f"q80_{name}"] = widths[(name, 80)]
        pred[f"q90_{name}"] = widths[(name, 90)]
    model = {"alpha": alpha, "ridge": ridge, "lgbm": lgbm, "x_tr": x_tr_f, "widths": widths, "features": features}
    return pred, model


def pearson(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 3:
        return np.nan
    x = np.asarray(a, dtype=np.float64)[m]
    y = np.asarray(b, dtype=np.float64)[m]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def auroc(actual, pred):
    label = np.asarray(actual) >= 0.5
    if label.size < 1 or np.unique(label).size < 2:
        return np.nan
    return float(roc_auc_score(label, pred))


def metric_block(actual, pred, q80, q90):
    actual = np.asarray(actual, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    q80 = np.asarray(q80, dtype=np.float64)
    q90 = np.asarray(q90, dtype=np.float64)
    m = np.isfinite(actual) & np.isfinite(pred)
    a = actual[m]
    p = pred[m]
    if q80.size == 1:
        q80 = np.full(m.sum(), q80.ravel()[0])
        q90 = np.full(m.sum(), q90.ravel()[0])
    else:
        q80 = q80[m]
        q90 = q90[m]
    return {
        "n": int(m.sum()),
        "spearman": spearman(p, a) if m.sum() else np.nan,
        "pearson": pearson(p, a) if m.sum() else np.nan,
        "mae": float(np.mean(np.abs(p - a))) if m.sum() else np.nan,
        "coverage_80": float(np.mean(np.abs(p - a) <= q80)) if m.sum() else np.nan,
        "width_80": float(np.mean(2 * q80)) if m.sum() else np.nan,
        "coverage_90": float(np.mean(np.abs(p - a) <= q90)) if m.sum() else np.nan,
        "width_90": float(np.mean(2 * q90)) if m.sum() else np.nan,
        "auroc": auroc(a, p) if m.sum() else np.nan,
        "mean_actual_minus_pred": float(np.mean(a - p)) if m.sum() else np.nan,
    }


def section_31(adc):
    preds = pd.read_csv(RES / "reliability_predictions.tsv", sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    min_c = int(PLAN["within_gene"]["min_defined_cancers"])
    min_center = int(PLAN["within_gene"]["gene_centered_min_cancers"])
    models = [("f0", "f0", "q80_f0", "q90_f0"), ("f1_ridge", "f1_ridge", "q80_ridge", "q90_ridge"), ("f1_lgbm", "f1_lgbm", "q80_lgbm", "q90_lgbm")]
    rows = []
    bias = []
    for universe, frame in (
        ("P", preds[preds["in_P"] & np.isfinite(preds["actual"])]),
        ("ADC", preds[preds["gene_symbol"].isin(adc) & np.isfinite(preds["actual"])]),
    ):
        for name, col, q80, q90 in models:
            rec = within_gene_summary(frame, col, min_c)
            rec.update({"universe": universe, "model": name, "metric": "within_gene"})
            if name == "f0":
                rec["note"] = PLAN["within_gene"]["f0_note"]
            else:
                rec["note"] = ""
            rows.append(rec)
            cen = centered_spearman(frame, col, min_center)
            cen.update({"universe": universe, "model": name, "metric": "gene_centered"})
            rows.append(cen)
            for hold, sub in frame.groupby("heldout_cancer"):
                block = metric_block(sub["actual"], sub[col], sub[q80], sub[q90])
                block.update({"universe": universe, "model": name, "heldout_cancer": hold})
                bias.append(block)
    pd.DataFrame(rows).to_csv(RES / "within_gene_stage2.tsv", sep="\t", index=False)
    pd.DataFrame(bias).to_csv(RES / "cancer_bias_stage2.tsv", sep="\t", index=False)
    print("section 3.1 written", flush=True)


def run_f1s(p_genes, adc, mean_pct, sd_pct, rate, rho):
    pred_rows = []
    for hold in CANCERS:
        source = [c for c in CANCERS if c != hold]
        pieces = []
        for inner in source:
            base = [c for c in source if c != inner]
            part = design_f1s(p_genes, base, inner, rho, mean_pct, sd_pct, rate, p_genes)
            part = part[np.isfinite(part["actual"])]
            pieces.append(part)
        train = pd.concat(pieces, ignore_index=True)
        apply_genes = list(dict.fromkeys(list(p_genes) + [g for g in adc if g not in set(p_genes)]))
        apply = design_f1s(apply_genes, source, hold, rho, mean_pct, sd_pct, rate, p_genes)
        pred, _ = fit_apply(train, apply, F1S_FEATURES)
        pred["heldout_cancer"] = hold
        pred["in_P"] = pred["gene_symbol"].isin(set(p_genes))
        pred_rows.append(pred)
        print(f"F1s {hold}", flush=True)
    preds = pd.concat(pred_rows, ignore_index=True)
    preds.to_csv(RES / "f1s_predictions.tsv", sep="\t", index=False)
    summarize_f1s(preds, adc)
    return preds


def summarize_f1s(preds, adc):
    models = [("f0", "f0", "q80_f0", "q90_f0"), ("f1s_ridge", "f1s_ridge", "q80_ridge", "q90_ridge"), ("f1s_lgbm", "f1s_lgbm", "q80_lgbm", "q90_lgbm")]
    base = preds[preds["in_P"] & np.isfinite(preds["actual"])]
    by = []
    overall = []
    for hold, sub in base.groupby("heldout_cancer"):
        for name, col, q80, q90 in models:
            rec = metric_block(sub["actual"], sub[col], sub[q80], sub[q90])
            rec.update({"heldout_cancer": hold, "model": name})
            by.append(rec)
    for name, col, q80, q90 in models:
        rec = metric_block(base["actual"], base[col], base[q80], base[q90])
        rec.update({"model": name})
        overall.append(rec)
    pd.DataFrame(by).to_csv(RES / "f1s_by_heldout.tsv", sep="\t", index=False)
    pd.DataFrame(overall).to_csv(RES / "f1s_overall.tsv", sep="\t", index=False)
    adc_base = preds[preds["gene_symbol"].isin(adc) & np.isfinite(preds["actual"])]
    adc_rows = []
    for name, col, q80, q90 in models:
        rec = metric_block(adc_base["actual"], adc_base[col], adc_base[q80], adc_base[q90])
        rec.update({"model": name, "heldout_cancer": "ALL"})
        adc_rows.append(rec)
        for hold, sub in adc_base.groupby("heldout_cancer"):
            rec = metric_block(sub["actual"], sub[col], sub[q80], sub[q90])
            rec.update({"model": name, "heldout_cancer": hold})
            adc_rows.append(rec)
    pd.DataFrame(adc_rows).to_csv(RES / "f1s_adc.tsv", sep="\t", index=False)
    min_c = int(PLAN["within_gene"]["min_defined_cancers"])
    min_center = int(PLAN["within_gene"]["gene_centered_min_cancers"])
    rows = []
    for universe, frame in (("P", base), ("ADC", adc_base)):
        for name, col, _, _ in models:
            rec = within_gene_summary(frame, col, min_c)
            rec.update({"universe": universe, "model": name, "metric": "within_gene", "note": PLAN["within_gene"]["f0_note"] if name == "f0" else ""})
            rows.append(rec)
            cen = centered_spearman(frame, col, min_center)
            cen.update({"universe": universe, "model": name, "metric": "gene_centered"})
            rows.append(cen)
    pd.DataFrame(rows).to_csv(RES / "within_gene_f1s.tsv", sep="\t", index=False)


def final_f1s_model(preds):
    train = preds[preds["in_P"].astype(bool) & np.isfinite(preds["actual"])].copy()
    alpha, _ = ridge_alpha_and_oof(train, F1S_FEATURES)
    x = train[F1S_FEATURES].to_numpy(dtype=np.float64)
    y = train["actual"].to_numpy(dtype=np.float64)
    x_f, _, _ = fill_median(x, x[:1])
    z, _ = standardize(x_f, x_f[:1])
    ridge = Ridge(alpha=alpha, fit_intercept=True)
    ridge.fit(z, y)
    lgbm = LGBMRegressor(**lgbm_params())
    lgbm.fit(x_f, y)
    widths = {}
    for name, col in (("f0", "f0"), ("ridge", "f1s_ridge"), ("lgbm", "f1s_lgbm")):
        resid = np.abs(train[col].to_numpy(dtype=np.float64) - y)
        widths[(name, 80)] = conformal_half_width(resid, 80)
        widths[(name, 90)] = conformal_half_width(resid, 90)
    pd.DataFrame(
        [{"item": "ridge_alpha", "value": alpha}, {"item": "n_train", "value": len(train)}]
        + [{"item": f"q{level}_{name}", "value": widths[(name, level)]} for name, level in widths]
    ).to_csv(RES / "f1s_model_fit.tsv", sep="\t", index=False)
    return {"alpha": alpha, "ridge": ridge, "lgbm": lgbm, "x_tr": x_f, "widths": widths}


def final_f1_model():
    preds = pd.read_csv(RES / "reliability_predictions.tsv", sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    train = preds[preds["in_P"] & np.isfinite(preds["actual"])].copy()
    alpha, _ = ridge_alpha_and_oof(train, FEATURE_NAMES)
    if float(alpha) != 100.0:
        raise RuntimeError(f"refit F1 alpha {alpha} does not match the frozen alpha 100")
    x = train[FEATURE_NAMES].to_numpy(dtype=np.float64)
    y = train["actual"].to_numpy(dtype=np.float64)
    x_f, _, _ = fill_median(x, x[:1])
    z, _ = standardize(x_f, x_f[:1])
    ridge = Ridge(alpha=alpha, fit_intercept=True)
    ridge.fit(z, y)
    lgbm = LGBMRegressor(**lgbm_params())
    lgbm.fit(x_f, y)
    fit = pd.read_csv(RES / "frozen_model_fit.tsv", sep="\t")
    stored = dict(zip(fit["item"], fit["value"]))
    widths = {}
    for name, col, key in (("f0", "f0", "f0"), ("ridge", "f1_ridge", "ridge"), ("lgbm", "f1_lgbm", "lgbm")):
        resid = np.abs(train[col].to_numpy(dtype=np.float64) - y)
        for level in (80, 90):
            widths[(name, level)] = conformal_half_width(resid, level)
            if abs(widths[(name, level)] - float(stored[f"q{level}_{key}"])) > 1e-8:
                raise RuntimeError(f"F1 width mismatch {name} {level}")
    return {"alpha": alpha, "ridge": ridge, "lgbm": lgbm, "x_tr": x_f, "widths": widths, "train": train}


def apply_saved(frame, features, fit, pred_names):
    x = frame[features].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    med = np.nanmedian(fit["x_tr"], axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(med, inds[1])
    mu = fit["x_tr"].mean(axis=0)
    sd = fit["x_tr"].std(axis=0, ddof=0)
    bad = (sd == 0) | ~np.isfinite(sd)
    sd = np.where(bad, 1.0, sd)
    z = (filled - mu) / sd
    z[:, bad] = 0.0
    out = frame.copy()
    out[pred_names[0]] = fit["ridge"].predict(z)
    out[pred_names[1]] = fit["lgbm"].predict(filled)
    return out


def ar_gene_list(p_genes, genes):
    gset = set(genes)
    pset = set(p_genes)
    rate = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    ok = rate[rate["n_measured"] >= int(PLAN["ar_set"]["min_measured"])]
    n_c = ok.groupby("gene_symbol").size()
    extra = [g for g, n in n_c.items() if n >= int(PLAN["ar_set"]["min_cancers"]) and g not in pset and g in gset]
    chosen = list(p_genes) + sorted(extra)
    pd.Series(chosen, name="gene_symbol").to_csv(RES / "ar_set_genes.txt", index=False, header=False)
    pd.DataFrame(
        [
            {"item": "n_P", "value": len(p_genes)},
            {"item": "n_extra", "value": len(extra)},
            {"item": "n_ar_set", "value": len(chosen)},
            {"item": "n_extra_not_in_G_RNA", "value": int(sum(g not in gset for g in n_c.index if g not in pset and n_c[g] >= 3))},
        ]
    ).to_csv(RES / "ar_set_summary.tsv", sep="\t", index=False)
    print("AR set", len(chosen), flush=True)
    return chosen


def load_tumor(genes, targets):
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "tumor"]
    data = {}
    all_ids = []
    for cancer in CANCERS:
        rna = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__rna.parquet")
        prot = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__protein.parquet")
        sids = [s for s in paired.loc[paired["cancer"] == cancer, "sample_id"] if s in rna.columns and s in prot.columns]
        x = zscore_cols(rna.reindex(index=genes, columns=sids).to_numpy(dtype=np.float64).T)
        y = zscore_cols(prot.reindex(index=targets, columns=sids).to_numpy(dtype=np.float64).T)
        data[cancer] = {"sample_id": sids, "x": x, "y": y}
        all_ids.extend(sids)
        print(f"loaded {cancer} n={len(sids)}", flush=True)
    if len(all_ids) != len(set(all_ids)):
        raise RuntimeError("tumor sample ids overlap across cancers")
    return data


def stack(data, cancers):
    x = np.vstack([data[c]["x"] for c in cancers])
    y = np.vstack([data[c]["y"] for c in cancers])
    labels = np.concatenate([np.full(data[c]["x"].shape[0], i, dtype=np.int32) for i, c in enumerate(cancers)])
    return x, y, labels


def fit_final_ar(genes, ar_genes):
    import importlib.util
    import multiprocessing as mp

    spec = importlib.util.spec_from_file_location("b2develop", ROOT / "scripts" / "06_develop.py")
    develop = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(develop)
    data = load_tumor(genes, ar_genes)
    gene_index = {g: i for i, g in enumerate(genes)}
    self_index = np.array([gene_index[g] for g in ar_genes], dtype=np.int32)
    x, y, labels = stack(data, CANCERS)
    ctx = mp.get_context("spawn")
    with ctx.Pool(os.cpu_count() or 1) as pool:
        alpha, lam, n_fb = develop.inner_scores(pool, x, y, labels, self_index, os.cpu_count() or 1)
    b_sel = fit_slope(x[:, self_index], y)
    e = residuals(x[:, self_index], y, b_sel)
    feats = select_residual_features(x, e, self_index)
    b, coefs = fit_ar_coefficients(x, y, self_index, feats, alpha, lam)
    rows = []
    for j, gene in enumerate(ar_genes):
        extras = [genes[i] for i in feats[:, j]]
        rows.append(
            {
                "gene_symbol": gene,
                "b": float(b[j]),
                "lambda": float(lam[j]),
                "alpha": float(alpha[j]),
                "n_fallback_in_model": int(n_fb),
                "extra_genes": extras,
                "coefficients": coefs[:, j].tolist(),
            }
        )
    # fallback count is the number of proteins whose inner folds were all invalid.
    # inner_scores returns that count, not a per-protein flag. Recompute from the
    # same rule by asking develop only for the count already returned.
    model = pd.DataFrame(rows)
    model.attrs["n_fallback"] = int(n_fb)
    model.to_parquet(RES / "ar_final_models.parquet", index=False)
    lam_counts = model["lambda"].value_counts().rename_axis("lambda").reset_index(name="n")
    alpha_counts = model["alpha"].value_counts().rename_axis("alpha").reset_index(name="n")
    lam_counts.to_csv(RES / "ar_final_lambda_counts.tsv", sep="\t", index=False)
    alpha_counts.to_csv(RES / "ar_final_alpha_counts.tsv", sep="\t", index=False)
    pd.DataFrame([{"item": "n_genes", "value": len(model)}, {"item": "n_fallback", "value": int(n_fb)}, {"item": "n_nonfinite_coef", "value": int((~np.isfinite(coefs)).any(axis=0).sum())}]).to_csv(
        RES / "ar_final_summary.tsv", sep="\t", index=False
    )
    print("AR final", len(model), "fallback", n_fb, flush=True)
    del data, x, y
    return model, feats, gene_index


def disease_samples():
    from collections import Counter

    overlap = pd.read_csv(RES / "tcga_rppa_rna_overlap_by_disease.tsv", sep="\t")
    keep = overlap.loc[overlap["n_ge_30"].astype(str).isin(["True", "true", "1"])]
    ph = pd.read_csv(RAW / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", dtype=str, usecols=["sample", "_primary_disease"])
    rna_path = RAW / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    with gzip.open(rna_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    rna_unique = {s for s, n in Counter(header[1:]).items() if n == 1}
    dup = set(pd.read_csv(RES / "rppa_duplicate_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    with gzip.open(RAW / "tcga" / "TCGA-RPPA-pancan-clean.xena.gz", "rt", encoding="utf-8", errors="replace") as handle:
        rppa_header = handle.readline().rstrip("\n").split("\t")
    both = {s for s in rppa_header[1:] if s not in dup} & rna_unique
    samples = {}
    for disease, n_expected in zip(keep["_primary_disease"], keep["n_samples"].astype(int)):
        ids = [s for s in ph.loc[ph["_primary_disease"] == disease, "sample"] if s in both]
        if len(ids) != int(n_expected):
            raise RuntimeError(f"sample count mismatch {disease}")
        samples[disease] = ids
    return samples, rna_path


def load_tcga_matrix(mapped, samples):
    symbol = dict(zip(mapped["row_id"], mapped["gene_symbol"]))
    rule = dict(zip(mapped["row_id"], mapped["rule"]))
    wanted = sorted({s for ids in samples.values() for s in ids})
    usecols = ["sample"] + wanted
    df = pd.read_csv(
        RAW / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
        sep="\t",
        usecols=usecols,
        dtype={"sample": str},
        low_memory=False,
    )
    df["gene_symbol"] = df["sample"].map(lambda rid: symbol.get(rid) if rule.get(rid) != "unmatched" else None)
    df = df[df["gene_symbol"].notna()]
    if df["gene_symbol"].duplicated().any():
        raise RuntimeError("TCGA map assigned one symbol more than once")
    mat = df.drop(columns=["sample"]).set_index("gene_symbol")
    return mat


def gene_mean_sd(mat, ids):
    sub = mat.reindex(columns=ids)
    arr = sub.to_numpy(dtype=np.float64)
    with np.errstate(all="ignore"):
        mu = np.nanmean(arr, axis=1)
        sd = np.nanstd(arr, axis=1, ddof=DDOF)
    return dict(zip(sub.index, mu)), dict(zip(sub.index, sd))


def z_for_genes(mat, ids, genes):
    """samples x genes. Absent genes are 0. Observed gaps stay NA."""
    n = len(ids)
    out = np.zeros((n, len(genes)), dtype=np.float64)
    present = []
    cols = []
    for j, gene in enumerate(genes):
        if gene not in mat.index:
            continue
        present.append(mat.loc[gene, ids].to_numpy(dtype=np.float64))
        cols.append(j)
    if present:
        block = np.column_stack(present)
        out[:, cols] = zscore_cols(block)
    return out


def external_f1s_rows(records, fit):
    frame = pd.DataFrame(records)
    scored = apply_saved(frame, F1S_FEATURES, fit, ("f1s_ridge", "f1s_lgbm"))
    scored["f0"] = frame["rho_mean"].to_numpy()
    for name, col in (("f0", "f0"), ("ridge", "f1s_ridge"), ("lgbm", "f1s_lgbm")):
        for level in (80, 90):
            half = fit["widths"][(name, level)]
            scored[f"{col}_low_{level}"] = scored[col] - half
            scored[f"{col}_high_{level}"] = scored[col] + half
            scored[f"{col}_half_{level}"] = half
    scored["ridge_alpha"] = fit["alpha"]
    return scored


def f1s_prediction_columns():
    cols = []
    for name in ("f0", "f1s_ridge", "f1s_lgbm"):
        cols.extend(
            [
                name,
                f"{name}_low_80",
                f"{name}_high_80",
                f"{name}_low_90",
                f"{name}_high_90",
                f"{name}_half_80",
                f"{name}_half_90",
            ]
        )
    return cols


def fit_and_freeze_ar(genes, p_genes, adc):
    ar_genes = ar_gene_list(p_genes, genes)
    model, _, gene_index = fit_final_ar(genes, ar_genes)
    return model


def models_from_parquet(path, genes):
    model = pd.read_parquet(path)
    gene_index = {g: i for i, g in enumerate(genes)}
    extras = []
    coefs = []
    self_index = []
    for row in model.itertuples(index=False):
        extra = list(row.extra_genes)
        extras.append([gene_index[g] for g in extra])
        coefs.append(list(row.coefficients))
        self_index.append(gene_index[row.gene_symbol])
    return model, np.array(self_index), np.array(extras), np.array(coefs, dtype=np.float64)


def write_scale_row(rows, context, bundle, stat, universe, values):
    rec = quantiles(values)
    rec.update({"context": context, "bundle": bundle, "stat": stat, "universe": universe})
    rows.append(rec)


def cptac_scale_rows(stats, p_genes):
    rows = []
    pset = set(p_genes)
    for cancer, sub in stats.groupby("cancer"):
        write_scale_row(rows, cancer, "CPTAC", "rna_mean", "G_RNA", sub["rna_mean"])
        write_scale_row(rows, cancer, "CPTAC", "rna_sd", "G_RNA", sub["rna_sd"])
        psub = sub[sub["gene_symbol"].isin(pset)]
        write_scale_row(rows, cancer, "CPTAC", "rna_mean", "P", psub["rna_mean"])
        write_scale_row(rows, cancer, "CPTAC", "rna_sd", "P", psub["rna_sd"])
    return rows


def training_pair_scale():
    preds = pd.read_csv(RES / "reliability_predictions.tsv", sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    base = preds[preds["in_P"] & np.isfinite(preds["actual"])]
    rows = []
    write_scale_row(rows, "CPTAC_LOCO_pairs", "training_F1_features", "holdout_rna_mean", "P", base["holdout_rna_mean"])
    write_scale_row(rows, "CPTAC_LOCO_pairs", "training_F1_features", "holdout_rna_sd", "P", base["holdout_rna_sd"])
    return rows


B2 = None
B2_LGBM = None


def main():
    global B2, B2_LGBM
    import yaml

    quiet_warnings()
    self_check_map()
    B2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
    B2_LGBM = dict(B2["reliability"]["lightgbm"])
    genes, p_genes, adc = load_lists()
    mapped = write_id_map(genes)
    section_31(adc)
    rho = defined_rho(p_genes)
    # ADC actuals used only at apply time.
    adc_loco = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    pset = set(p_genes)
    for r in adc_loco.itertuples(index=False):
        if r.gene_symbol in pset or r.gene_symbol not in set(adc):
            continue
        if int(r.n) >= MIN_N and np.isfinite(r.s0):
            rho[(r.heldout_cancer, r.gene_symbol)] = float(r.s0)
    stats, mean_pct, sd_pct, _ = cptac_percentile_tables(genes)
    rates = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    rate = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
    print("F1s start", flush=True)
    f1s_preds = run_f1s(p_genes, adc, mean_pct, sd_pct, rate, rho)
    f1s_fit = final_f1s_model(f1s_preds)
    f1_fit = final_f1_model()
    print("AR start", flush=True)
    model = fit_and_freeze_ar(genes, p_genes, adc)
    # The rest of the external application is in the second half, after the
    # CPTAC matrices are released.
    del model
    freeze_external(genes, p_genes, adc, mapped, mean_pct, sd_pct, rate, rho, f1s_fit, f1_fit, stats)
    print("stage3 develop done", flush=True)


def gmean_all(rho, p_genes):
    vals = [v for (c, g), v in rho.items() if c in CANCERS and g in set(p_genes)]
    return float(np.mean(vals))


def forecast_rho(rho, genes):
    """Same CPTAC S0 values stage 2 used for forecast genes outside P."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("b2freeze", ROOT / "scripts" / "07_freeze.py")
    freeze = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(freeze)
    antibodies = pd.read_csv(RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    antibodies = antibodies[antibodies["kept_for_forecast"].astype(str).isin(["True", "true"])]
    proteins = pd.read_csv(RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    proteins = proteins[(proteins["symbol"] != "NA") & (proteins["symbol"] != "") & proteins["symbol"].isin(set(genes))]
    forecast = sorted(set(antibodies["gene_symbol"]) | set(proteins["symbol"]))
    rho, n_added = freeze.extra_cptac_rho(forecast, rho)
    print(f"extra cptac rho pairs {n_added}", flush=True)
    return rho


def assert_same_f0(new, old_path, keys):
    old = pd.read_csv(old_path, sep="\t")
    merged = new.merge(old[keys + ["f0"]], on=keys, suffixes=("_new", "_old"))
    if len(merged) != len(new) or len(merged) != len(old):
        raise RuntimeError(f"F0 row alignment failed for {old_path}")
    diff = (merged["f0_new"] - merged["f0_old"]).abs()
    if float(diff.max()) > 1e-6:
        raise RuntimeError(f"F0 does not match the stage-2 file {old_path}")


def freeze_external(genes, p_genes, adc, mapped, mean_pct, sd_pct, rate, rho, f1s_fit, f1_fit, stats):
    samples, _ = disease_samples()
    print("reading TCGA RNA", flush=True)
    mat = load_tcga_matrix(mapped, samples)
    # Stage-2 fill check: a gene with no stage-2 symbol, or a non-finite mean.
    stage2_missing = set(pd.read_csv(RES / "stage2_rppa_genes_without_rna_row.tsv", sep="\t", dtype=str)["gene_symbol"]) if (RES / "stage2_rppa_genes_without_rna_row.tsv").exists() else set()
    antibodies = pd.read_csv(RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    fill_rows = []
    for disease, ids in samples.items():
        mu, sd = gene_mean_sd(mat, ids)
        for gene in antibodies["gene_symbol"]:
            filled = gene in stage2_missing or not np.isfinite(mu.get(gene, np.nan)) or not np.isfinite(sd.get(gene, np.nan))
            # Stage-2 map did not include previous symbols, so a gene reached only
            # by prev_symbol was missing then even if it is in `mat` now.
            if gene in stage2_missing:
                filled = True
            fill_rows.append({"gene_symbol": gene, "cancer": disease, "stage2_rna_feature_median_filled": bool(filled)})
    fill_df = pd.DataFrame(fill_rows)
    fill_df.to_csv(RES / "stage2_f1_rppa_fill.tsv", sep="\t", index=False)
    print("stage2 filled rows", int(fill_df["stage2_rna_feature_median_filled"].sum()), flush=True)
    scale_rows = cptac_scale_rows(stats, p_genes) + training_pair_scale()
    pset = set(p_genes)
    adc_set = set(adc)
    rho = forecast_rho(rho, genes)
    gmean = gmean_all(rho, p_genes)
    model, self_index, feat_index, coefs = models_from_parquet(RES / "ar_final_models.parquet", genes)
    ar_symbols = set(model["gene_symbol"])
    kept = antibodies[antibodies["gene_symbol"].isin(ar_symbols)].copy()
    missing_ab = antibodies.loc[~antibodies["gene_symbol"].isin(ar_symbols), "xena_antibody"].tolist()
    pd.DataFrame({"xena_antibody": missing_ab}).to_csv(RES / "rppa_antibodies_not_in_ar_set.tsv", sep="\t", index=False)
    available = set(mat.index)
    extra_missing = []
    rppa_pred_rows = []
    f1s_records = []
    b = model["b"].to_numpy(dtype=np.float64)
    lams = model["lambda"].to_numpy(dtype=np.float64)
    # Restrict AR application to the kept antibodies, in model-row order.
    keep_pos = [i for i, g in enumerate(model["gene_symbol"]) if g in set(kept["gene_symbol"])]
    ab_for_gene = {r.gene_symbol: r.xena_antibody for r in kept.itertuples(index=False)}
    for disease, ids in samples.items():
        mu, sd = gene_mean_sd(mat, ids)
        mu_pct = percentile_map(mu)
        sd_pct = percentile_map(sd)
        for gene in antibodies["gene_symbol"]:
            ext_mu, ext_mu_pct = mu_pct.get(gene, (np.nan, np.nan))
            ext_sd, ext_sd_pct = sd_pct.get(gene, (np.nan, np.nan))
            rec = f1s_feature_row(gene, CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, ext_sd_pct)
            rec.update({"xena_antibody": ab_for_gene.get(gene, antibodies.loc[antibodies.gene_symbol == gene, "xena_antibody"].iloc[0]), "gene_symbol": gene, "cancer": disease, "n": len(ids), "is_adc": gene in adc_set})
            f1s_records.append(rec)
        write_scale_row(scale_rows, disease, "TCGA", "rna_mean", "G_RNA", list(mu.values()))
        write_scale_row(scale_rows, disease, "TCGA", "rna_sd", "G_RNA", list(sd.values()))
        write_scale_row(scale_rows, disease, "TCGA", "rna_mean", "P", [mu.get(g, np.nan) for g in p_genes])
        write_scale_row(scale_rows, disease, "TCGA", "rna_sd", "P", [sd.get(g, np.nan) for g in p_genes])
        if not keep_pos:
            continue
        needed = []
        pos_of = {}
        for j in keep_pos:
            for g_i in [self_index[j], *feat_index[j]]:
                g = genes[g_i]
                if g not in pos_of:
                    pos_of[g] = len(needed)
                    needed.append(g)
        z_need = z_for_genes(mat, ids, needed)
        # Rebuild a full G_RNA-shaped z only for the columns the predictor indexes.
        z_full = np.zeros((len(ids), len(genes)), dtype=np.float64)
        for g, col in pos_of.items():
            z_full[:, gene_index_of(genes, g)] = z_need[:, col]
        z_full = np.nan_to_num(z_full, nan=0.0)
        pred = predict_ar_matrix(
            z_full,
            self_index[keep_pos],
            feat_index[keep_pos].T,
            b[keep_pos],
            lams[keep_pos],
            np.asarray(coefs[keep_pos], dtype=np.float64).T,
        )
        for t, j in enumerate(keep_pos):
            gene = model["gene_symbol"].iloc[j]
            n_miss = int(sum(genes[g_i] not in available for g_i in feat_index[j]))
            extra_missing.append({"cohort": "TCGA", "context": "all_diseases", "gene_symbol": gene, "n_extra_missing": n_miss})
            for s_i, sample in enumerate(ids):
                rppa_pred_rows.append(
                    {
                        "xena_antibody": ab_for_gene[gene],
                        "gene_symbol": gene,
                        "cancer": disease,
                        "sample_id": sample,
                        "rna_z": z_full[s_i, self_index[j]],
                        "ar_pred": pred[s_i, t],
                    }
                )
        print(f"TCGA AR {disease}", flush=True)
    # The extra-missing loop above repeats the same gene once per disease. Keep one.
    extra_df = pd.DataFrame(extra_missing).drop_duplicates(["cohort", "gene_symbol"])
    pd.DataFrame(rppa_pred_rows).to_parquet(RES / "frozen_ar_tcga.parquet", index=False)
    f1s_rppa = external_f1s_rows(f1s_records, f1s_fit)
    rppa_cols = ["xena_antibody", "gene_symbol", "cancer", "n", "is_adc", "ridge_alpha"] + f1s_prediction_columns()
    assert_same_f0(f1s_rppa, RES / "frozen_forecasts_rppa.tsv", ["xena_antibody", "gene_symbol", "cancer"])
    f1s_rppa[rppa_cols].to_csv(RES / "frozen_forecasts_rppa_F1s.tsv", sep="\t", index=False)
    del mat
    procan_scale, procan_extra, procan_ar, procan_f1s = freeze_procan(genes, p_genes, adc, model, self_index, feat_index, coefs, b, lams, mean_pct, sd_pct, rate, rho, gmean, f1s_fit)
    scale_rows.extend(procan_scale)
    extra_df = pd.concat([extra_df, procan_extra], ignore_index=True)
    extra_df.to_csv(RES / "ar_extra_missing.tsv", sep="\t", index=False)
    pd.DataFrame(scale_rows).to_csv(RES / "rna_scale_quantiles.tsv", sep="\t", index=False)
    ccle_check(genes, p_genes, adc, model, self_index, feat_index, coefs, b, lams, mean_pct, sd_pct, rate, rho, gmean, f1_fit, f1s_fit, scale_rows)
    write_hashes()


def gene_index_of(genes, gene):
    return genes.index(gene)


def freeze_procan(genes, p_genes, adc, model, self_index, feat_index, coefs, b, lams, mean_pct, sd_pct, rate, rho, gmean, f1s_fit):
    import csv

    lin_table = pd.read_csv(RES / "procan_sanger_by_lineage.tsv", sep="\t")
    lineages = lin_table.loc[lin_table["n_ge_30"].astype(str).isin(["True", "true", "1"]), "Tissue_type"].tolist()
    n_by = dict(zip(lin_table["Tissue_type"], lin_table["n_cell_lines"].astype(int)))
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
            r = mindex[model_id]
            c = gindex[gene]
            if filled[r, c]:
                continue
            try:
                tpm = float(parts[i_tpm])
            except ValueError:
                continue
            if not np.isfinite(tpm):
                continue
            arr[r, c] = np.float32(np.log2(tpm + 1.0))
            filled[r, c] = True
    use_rows = [i for i, sidm in enumerate(cands) if sidm in sanger_ids]
    use_ids = [cands[i] for i in use_rows]
    arr = arr[use_rows]
    sidm_lineage = dict(zip(mapping["SIDM"], mapping["Tissue_type"]))
    proteins = pd.read_csv(RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    proteins = proteins[(proteins["symbol"] != "NA") & (proteins["symbol"] != "") & proteins["symbol"].isin(set(genes))]
    ar_symbols = set(model["gene_symbol"])
    proteins_ar = proteins[proteins["symbol"].isin(ar_symbols)]
    keep_pos = [i for i, g in enumerate(model["gene_symbol"]) if g in set(proteins_ar["symbol"])]
    pset = set(p_genes)
    adc_set = set(adc)
    scale_rows = []
    extra_rows = []
    ar_rows = []
    f1s_records = []
    by_lineage_index = {}
    for lin in lineages:
        idx = [i for i, sidm in enumerate(use_ids) if sidm_lineage[sidm] == lin]
        by_lineage_index[lin] = idx
        block = arr[idx]
        with np.errstate(all="ignore"):
            mu = np.nanmean(block, axis=0)
            sd = np.nanstd(block, axis=0, ddof=DDOF)
        mu_d = {genes[i]: float(mu[i]) if np.isfinite(mu[i]) else np.nan for i in range(len(genes))}
        sd_d = {genes[i]: float(sd[i]) if np.isfinite(sd[i]) else np.nan for i in range(len(genes))}
        # A gene with no observation has a nan mean. Percentiles use finite values.
        mu_pct = percentile_map(mu_d)
        sd_pct = percentile_map(sd_d)
        write_scale_row(scale_rows, lin, "Sanger", "rna_mean", "G_RNA", list(mu_d.values()))
        write_scale_row(scale_rows, lin, "Sanger", "rna_sd", "G_RNA", list(sd_d.values()))
        write_scale_row(scale_rows, lin, "Sanger", "rna_mean", "P", [mu_d.get(g, np.nan) for g in p_genes])
        write_scale_row(scale_rows, lin, "Sanger", "rna_sd", "P", [sd_d.get(g, np.nan) for g in p_genes])
        available = {g for g, v in mu_d.items() if np.isfinite(v)}
        for j in keep_pos:
            n_miss = int(sum(genes[g_i] not in available for g_i in feat_index[j]))
            extra_rows.append({"cohort": "ProCan", "context": lin, "gene_symbol": model["gene_symbol"].iloc[j], "n_extra_missing": n_miss})
        for row in proteins.itertuples(index=False):
            ext_mu_pct = mu_pct.get(row.symbol, (np.nan, np.nan))[1]
            ext_sd_pct = sd_pct.get(row.symbol, (np.nan, np.nan))[1]
            rec = f1s_feature_row(row.symbol, CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, ext_sd_pct)
            rec.update(
                {
                    "protein_id": row.protein_id,
                    "gene_symbol": row.symbol,
                    "lineage": lin,
                    "n": int(n_by[lin]),
                    "n_rna_gene": int(np.isfinite(block[:, gindex[row.symbol]]).sum()) if row.symbol in gindex else 0,
                    "is_adc": row.symbol in adc_set,
                }
            )
            f1s_records.append(rec)
        # AR predictions for proteins in the AR set.
        if not idx or not keep_pos:
            continue
        raw = np.array(block, dtype=np.float64)
        # Absent genes (all nan) become 0 before z-scoring. Sample gaps stay nan.
        absent = ~np.isfinite(raw).any(axis=0)
        raw[:, absent] = 0.0
        z = zscore_cols(raw)
        z = np.nan_to_num(z, nan=0.0)
        pred = predict_ar_matrix(
            z,
            self_index[keep_pos],
            feat_index[keep_pos].T,
            b[keep_pos],
            lams[keep_pos],
            np.asarray(coefs[keep_pos], dtype=np.float64).T,
        )
        id_to_protein = {}
        for row in proteins_ar.itertuples(index=False):
            id_to_protein.setdefault(row.symbol, []).append(row.protein_id)
        for t, j in enumerate(keep_pos):
            gene = model["gene_symbol"].iloc[j]
            for protein_id in id_to_protein.get(gene, []):
                for local_i, row_i in enumerate(idx):
                    ar_rows.append(
                        {
                            "protein_id": protein_id,
                            "gene_symbol": gene,
                            "lineage": lin,
                            "model_id": use_ids[row_i],
                            "rna_z": z[local_i, self_index[j]],
                            "ar_pred": pred[local_i, t],
                        }
                    )
        print(f"ProCan AR {lin} models={len(idx)}", flush=True)
    pd.DataFrame(ar_rows).to_parquet(RES / "frozen_ar_procan.parquet", index=False)
    f1s = external_f1s_rows(f1s_records, f1s_fit)
    cols = ["protein_id", "gene_symbol", "lineage", "n", "n_rna_gene", "is_adc", "ridge_alpha"] + f1s_prediction_columns()
    assert_same_f0(f1s, RES / "frozen_forecasts_procan.tsv", ["protein_id", "gene_symbol", "lineage"])
    f1s[cols].to_csv(RES / "frozen_forecasts_procan_F1s.tsv", sep="\t", index=False)
    return scale_rows, pd.DataFrame(extra_rows), ar_rows, f1s


def ccle_check(genes, p_genes, adc, model, self_index, feat_index, coefs, b, lams, mean_pct, sd_pct, rate, rho, gmean, f1_fit, f1s_fit, scale_rows):
    meta = pd.read_csv(PROC / "ccle" / "models.tsv", sep="\t", dtype=str)
    rna = pd.read_parquet(PROC / "ccle" / "rna.parquet")
    prot = pd.read_parquet(PROC / "ccle" / "protein.parquet")
    meta = meta[meta["ModelID"].isin(set(rna.columns) & set(prot.columns))].copy()
    gindex = {g: i for i, g in enumerate(genes)}
    pset = set(p_genes)
    delta_rows = []
    rel_rows = []
    stash = []

    def append_reliability(frame, lineage):
        for universe, mask in (
            ("P", frame["gene_symbol"].isin(pset)),
            ("ADC", frame["gene_symbol"].isin(adc)),
            ("ALL_EVAL", np.ones(len(frame), dtype=bool)),
        ):
            for model_name in ("f0", "f1_ridge", "f1_lgbm", "f1s_ridge", "f1s_lgbm"):
                if model_name == "f0":
                    q80, q90 = f1_fit["widths"][("f0", 80)], f1_fit["widths"][("f0", 90)]
                elif model_name.startswith("f1s"):
                    key = "ridge" if "ridge" in model_name else "lgbm"
                    q80, q90 = f1s_fit["widths"][(key, 80)], f1s_fit["widths"][(key, 90)]
                else:
                    key = "ridge" if "ridge" in model_name else "lgbm"
                    q80, q90 = f1_fit["widths"][(key, 80)], f1_fit["widths"][(key, 90)]
                block = metric_block(frame.loc[mask, "actual"], frame.loc[mask, model_name], np.array([q80]), np.array([q90]))
                block.update({"lineage": lineage, "universe": universe, "model": model_name})
                rel_rows.append(block)

    for lin in PLAN["ccle"]["lineages"]:
        ids = meta.loc[meta["OncotreeLineage"] == lin, "ModelID"].tolist()
        raw = np.zeros((len(ids), len(genes)), dtype=np.float64)
        present = set(rna.index)
        for j, gene in enumerate(genes):
            if gene in present:
                raw[:, j] = pd.to_numeric(rna.loc[gene, ids], errors="coerce").to_numpy(dtype=np.float64)
            else:
                raw[:, j] = np.nan
        # Reliability features use raw values. Absent genes stay missing.
        with np.errstate(all="ignore"):
            mu = np.nanmean(raw, axis=0)
            sd = np.nanstd(raw, axis=0, ddof=DDOF)
        mu_d = {genes[i]: float(mu[i]) if np.isfinite(mu[i]) else np.nan for i in range(len(genes))}
        sd_d = {genes[i]: float(sd[i]) if np.isfinite(sd[i]) else np.nan for i in range(len(genes))}
        write_scale_row(scale_rows, lin, "CCLE", "rna_mean", "G_RNA", list(mu_d.values()))
        write_scale_row(scale_rows, lin, "CCLE", "rna_sd", "G_RNA", list(sd_d.values()))
        write_scale_row(scale_rows, lin, "CCLE", "rna_mean", "P", [mu_d.get(g, np.nan) for g in p_genes])
        write_scale_row(scale_rows, lin, "CCLE", "rna_sd", "P", [sd_d.get(g, np.nan) for g in p_genes])
        mu_pct = percentile_map({g: v for g, v in mu_d.items() if g in present or np.isfinite(v)})
        sd_pct_c = percentile_map({g: v for g, v in sd_d.items() if g in present or np.isfinite(v)})
        # AR z: absent gene is 0, then z-score.
        raw_z = raw.copy()
        raw_z[:, ~np.isfinite(raw_z).any(axis=0)] = 0.0
        z = zscore_cols(raw_z)
        z_pred = np.nan_to_num(z, nan=0.0)
        pred = predict_ar_matrix(z_pred, self_index, feat_index.T, b, lams, np.asarray(coefs, dtype=np.float64).T)
        # Protein z and S0 for genes in the AR set that are measured.
        for j, gene in enumerate(model["gene_symbol"]):
            if gene not in prot.index:
                continue
            yy = pd.to_numeric(prot.loc[gene, ids], errors="coerce").to_numpy(dtype=np.float64)
            n = int(np.isfinite(yy).sum())
            if n < MIN_N:
                continue
            yz = zscore_cols(yy.reshape(-1, 1)).ravel()
            s0 = spearman(z[:, self_index[j]], yz)
            ar = spearman(pred[:, j], yz)
            delta_rows.append({"lineage": lin, "gene_symbol": gene, "in_P": gene in pset, "n": n, "s0": s0, "ar": ar, "delta_ar": ar - s0 if np.isfinite(ar) and np.isfinite(s0) else np.nan})
        # Reliability on P genes with defined S0, plus ADC rows stored in the same file.
        eval_genes = list(dict.fromkeys(list(p_genes) + list(adc)))
        f1_recs = []
        f1s_recs = []
        actuals = []
        # Reuse CPTAC feature builders. External RNA stats replace the holdout fields.
        # F1 raw mean/sd come from mu_d/sd_d. Percentile of sd uses the CCLE percentile.
        cptac_sd_pct = {}
        stats = pd.read_csv(RES / "cptac_rna_gene_stats.tsv", sep="\t")
        for cancer, sub in stats.groupby("cancer"):
            mapped = percentile_map(dict(zip(sub["gene_symbol"], sub["rna_sd"])))
            means = dict(zip(sub["gene_symbol"], sub["rna_mean"]))
            for g in eval_genes:
                cptac_sd_pct[(cancer, g)] = mapped.get(g, (np.nan, np.nan))
                # mean stored below
            cptac_sd_pct[("mean", cancer)] = means
        for gene in eval_genes:
            if gene not in prot.index:
                actual = np.nan
                n = 0
            else:
                yy = pd.to_numeric(prot.loc[gene, ids], errors="coerce").to_numpy(dtype=np.float64)
                n = int(np.isfinite(yy).sum())
                yz = zscore_cols(yy.reshape(-1, 1)).ravel()
                actual = spearman(z[:, gindex[gene]], yz) if n >= MIN_N and gene in gindex else np.nan
            actuals.append(actual)
            # F1 feature row, stage-2 names.
            vals = [rho[(c, gene)] for c in CANCERS if (c, gene) in rho]
            if vals:
                rho_mean, rho_sd, rho_min, rho_max, n_def = float(np.mean(vals)), float(np.std(vals, ddof=0)) if len(vals) >= 2 else 0.0, float(np.min(vals)), float(np.max(vals)), len(vals)
            else:
                rho_mean = gmean
                rho_sd, rho_min, rho_max, n_def = 0.0, gmean, gmean, 0
            means = [cptac_sd_pct[("mean", c)].get(gene, np.nan) for c in CANCERS]
            means = [v for v in means if np.isfinite(v)]
            rates = [rate.get((c, gene), 0.0) for c in CANCERS]
            hold_sd, hold_pct = sd_pct_c.get(gene, (np.nan, np.nan))
            f1_recs.append(
                {
                    "rho_mean": rho_mean,
                    "rho_sd": rho_sd,
                    "rho_min": rho_min,
                    "rho_max": rho_max,
                    "n_defined_cancers": n_def,
                    "holdout_rna_mean": mu_d.get(gene, np.nan),
                    "holdout_rna_sd": hold_sd,
                    "holdout_rna_sd_percentile": hold_pct,
                    "train_rna_mean": float(np.mean(means)) if means else np.nan,
                    "train_protein_rate_mean": float(np.mean(rates)) if rates else np.nan,
                    "gene_symbol": gene,
                    "actual": actual,
                }
            )
            f1s_recs.append(f1s_feature_row(gene, CANCERS, rho, mean_pct, sd_pct, rate, gmean, mu_pct.get(gene, (np.nan, np.nan))[1], sd_pct_c.get(gene, (np.nan, np.nan))[1]) | {"gene_symbol": gene, "actual": actual, "query_cancer": lin})
        f1_frame = pd.DataFrame(f1_recs)
        f1_hat = apply_saved(f1_frame, FEATURE_NAMES, f1_fit, ("f1_ridge", "f1_lgbm"))
        f1_hat["f0"] = f1_frame["rho_mean"]
        f1s_frame = pd.DataFrame(f1s_recs)
        f1s_hat = apply_saved(f1s_frame, F1S_FEATURES, f1s_fit, ("f1s_ridge", "f1s_lgbm"))
        piece = pd.DataFrame(
            {
                "gene_symbol": f1_frame["gene_symbol"].to_numpy(),
                "actual": f1_frame["actual"].to_numpy(),
                "f0": f1_hat["f0"].to_numpy(),
                "f1_ridge": f1_hat["f1_ridge"].to_numpy(),
                "f1_lgbm": f1_hat["f1_lgbm"].to_numpy(),
                "f1s_ridge": f1s_hat["f1s_ridge"].to_numpy(),
                "f1s_lgbm": f1s_hat["f1s_lgbm"].to_numpy(),
            }
        )
        stash.append(piece)
        append_reliability(piece, lin)
        print(f"CCLE {lin} n={len(ids)}", flush=True)
    if stash:
        append_reliability(pd.concat(stash, ignore_index=True), "pooled")
    deltas = pd.DataFrame(delta_rows)
    deltas.to_csv(RES / "ccle_ar_pairs.tsv", sep="\t", index=False)
    summary = []
    for lin, sub in deltas.groupby("lineage"):
        for universe, part in (("P", sub[sub["in_P"]]), ("AR_set", sub)):
            d = part["delta_ar"].dropna()
            summary.append(
                {
                    "lineage": lin,
                    "universe": universe,
                    "n": int(d.size),
                    "median": float(d.median()) if len(d) else np.nan,
                    "q1": float(d.quantile(0.25)) if len(d) else np.nan,
                    "q3": float(d.quantile(0.75)) if len(d) else np.nan,
                    "frac_gt_0.05": float((d > 0.05).mean()) if len(d) else np.nan,
                    "frac_lt_m0.05": float((d < -0.05).mean()) if len(d) else np.nan,
                }
            )
    pd.DataFrame(summary).to_csv(RES / "ccle_ar_delta.tsv", sep="\t", index=False)
    pd.DataFrame(rel_rows).to_csv(RES / "ccle_reliability.tsv", sep="\t", index=False)
    pd.DataFrame(scale_rows).to_csv(RES / "rna_scale_quantiles.tsv", sep="\t", index=False)
    write_overlap(meta)
    print("CCLE done", flush=True)


def write_overlap(meta):
    mapping = pd.read_csv(RAW / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t", dtype=str)
    sample_ids = set(pd.read_csv(RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    mapping = mapping[mapping["Project_Identifier"].isin(sample_ids)]
    cmp = pd.read_csv(RAW / "cmp" / "model_list_20260921.csv", usecols=["model_id", "BROAD_ID", "tissue"], dtype=str)
    joined = mapping.merge(cmp, left_on="SIDM", right_on="model_id", how="left")
    rows = []
    for lin in PLAN["ccle"]["lineages"]:
        ccle_ids = set(meta.loc[meta["OncotreeLineage"] == lin, "ModelID"])
        procan = joined[joined["Tissue_type"] == lin]
        both = procan[procan["BROAD_ID"].isin(ccle_ids)]
        rows.append({"lineage": lin, "n_ccle": len(ccle_ids), "n_procan": int(procan.shape[0]), "n_id_overlap_same_lineage": int(both["SIDM"].nunique())})
    pd.DataFrame(rows).to_csv(RES / "ccle_procan_overlap.tsv", sep="\t", index=False)


def write_hashes():
    rows = []
    for key, rel in (
        ("ar_models", "results/ar_final_models.parquet"),
        ("ar_tcga", "results/frozen_ar_tcga.parquet"),
        ("ar_procan", "results/frozen_ar_procan.parquet"),
        ("f1s_rppa", "results/frozen_forecasts_rppa_F1s.tsv"),
        ("f1s_procan", "results/frozen_forecasts_procan_F1s.tsv"),
    ):
        path = ROOT / rel
        rows.append({"key": key, "path": rel, "sha256": sha256_file(path), "size": path.stat().st_size})
    append_frozen_hashes(rows)
    print("hashes", rows, flush=True)


if __name__ == "__main__":
    main()
