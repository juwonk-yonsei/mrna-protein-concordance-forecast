"""Freeze R1/R2 weights and APOLLO predictions. Does not read APOLLO protein values.

R1/R2 are refit with the stage-5 functions. RPPA and ProCan deltas must match
the stage-5 table within 1e-9 before any APOLLO prediction is written.
"""
import importlib.util
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.b2_lib import FEATURE_NAMES, RES, ROOT
from scripts.b3_lib import F1S_FEATURES, percentile_map, predict_ar_matrix
from scripts.common import sha256_file

OUT = RES / "stage6" / "frozen"
PREP = RES / "stage6" / "prep"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def unit_means(rows):
    frame = pd.DataFrame(rows)
    out = []
    for platform, sub in frame.groupby("platform", sort=True):
        for model, col in (("R1", "delta_r1"), ("R2", "delta_r2")):
            for unit, u in sub.groupby("unit", sort=True):
                d = u[col].to_numpy(dtype=np.float64)
                d = d[np.isfinite(d)]
                if d.size:
                    out.append({"platform": platform, "unit": unit, "model": model, "mean_delta": float(d.mean())})
    return pd.DataFrame(out)


def compare_stage5(rows):
    new = unit_means(rows)
    old = pd.read_csv(RES / "stage5" / "posthoc" / "posthoc_23_units.tsv", sep="\t")
    old = old[old["model"].isin(["R1", "R2"])][["platform", "unit", "model", "mean_delta"]]
    merged = new.merge(old, on=["platform", "unit", "model"], how="outer", suffixes=("_new", "_old"))
    both = merged[np.isfinite(merged["mean_delta_new"]) & np.isfinite(merged["mean_delta_old"])]
    gap = np.abs(both["mean_delta_new"] - both["mean_delta_old"])
    rec = {
        "n_new": int(len(new)),
        "n_old": int(len(old)),
        "n_both": int(len(both)),
        "n_only_new": int(merged["mean_delta_old"].isna().sum()),
        "n_only_old": int(merged["mean_delta_new"].isna().sum()),
        "max_abs": float(gap.max()) if len(gap) else np.nan,
    }
    print("stage5 match", rec, flush=True)
    if rec["n_only_new"] or rec["n_only_old"] or not np.isfinite(rec["max_abs"]) or rec["max_abs"] > 1e-9:
        raise RuntimeError(f"R1/R2 predictions do not match stage 5: {rec}")
    return rec


def fit_r1r2(r1):
    r1.self_check()
    plan = r1.load_plan()
    cancers = list(plan["cancers"])
    alphas = [float(a) for a in plan["r1"]["alphas"]]
    genes = r1.gene_list("G_RNA_genes.txt")
    p_genes = r1.gene_list("P_genes.txt")
    gindex = {g: i for i, g in enumerate(genes)}
    self_index = np.array([gindex[g] for g in p_genes], dtype=np.int32)
    data = r1.load_tumor(genes, p_genes, cancers)
    x, y, labels = r1.stack(data, cancers)
    del data
    print("stacked", x.shape, y.shape, flush=True)
    feats = r1.select_features(x, y, self_index)
    p = y.shape[1]
    xf = np.empty((p, x.shape[0], feats.shape[0]), dtype=np.float64)
    for j in range(p):
        xf[j] = x[:, feats[:, j]]
    fitted = r1.Parallel(n_jobs=-1, backend="loky")(r1.delayed(r1.fit_one)(xf[j], y[:, j], labels, alphas) for j in range(p))
    chosen = np.array([rec[0] for rec in fitted], dtype=np.float64)
    weights = np.vstack(r1.Parallel(n_jobs=-1, backend="loky")(r1.delayed(r1.final_weights)(xf[j], y[:, j], chosen[j]) for j in range(p))).T
    del xf
    alpha_r2, _means, r2_fallback = r1.r2_alpha(x, y, labels, alphas)
    y_fill = np.nan_to_num(y, nan=0.0)
    print("r1 fallback", int(sum(rec[1] for rec in fitted)), "r2", alpha_r2, r2_fallback, flush=True)
    return genes, p_genes, feats, weights, chosen, x, y_fill, float(alpha_r2)


def fit_reliability(stage3):
    import yaml
    from lightgbm import LGBMRegressor
    from sklearn.linear_model import Ridge
    from scripts.b2_lib import conformal_half_width

    b2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
    stage3.B2 = b2
    stage3.B2_LGBM = dict(b2["reliability"]["lightgbm"])
    preds = pd.read_csv(RES / "f1s_predictions.tsv", sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    train = preds[preds["in_P"] & np.isfinite(preds["actual"])].copy()
    stored = dict(zip(pd.read_csv(RES / "f1s_model_fit.tsv", sep="\t")["item"], pd.read_csv(RES / "f1s_model_fit.tsv", sep="\t")["value"]))
    alpha = 100.0
    x = train[F1S_FEATURES].to_numpy(dtype=np.float64)
    y = train["actual"].to_numpy(dtype=np.float64)
    x_f, _, _ = stage3.fill_median(x, x[:1])
    z, _ = stage3.standardize(x_f, x_f[:1])
    ridge = Ridge(alpha=alpha, fit_intercept=True)
    ridge.fit(z, y)
    lgbm = LGBMRegressor(**stage3.lgbm_params())
    lgbm.fit(x_f, y)
    widths = {}
    for name, col, key in (("f0", "f0", "f0"), ("ridge", "f1s_ridge", "ridge"), ("lgbm", "f1s_lgbm", "lgbm")):
        resid = np.abs(train[col].to_numpy(dtype=np.float64) - y)
        for level in (80, 90):
            half = conformal_half_width(resid, level)
            if abs(half - float(stored[f"q{level}_{key}"])) > 1e-8:
                raise RuntimeError(f"F1s width changed {key} {level}")
            widths[(name, level)] = half
    f1s_fit = {"alpha": alpha, "ridge": ridge, "lgbm": lgbm, "x_tr": x_f, "widths": widths}
    f1_fit = stage3.final_f1_model()
    return f1s_fit, f1_fit


def apollo_z(stage3, genes):
    mat = pd.read_parquet(ROOT / "data" / "processed" / "apollo_luad" / "rna_log2tpm.parquet")
    samples = mat.columns.astype(str).tolist()
    z = stage3.z_for_genes(mat, samples, genes)
    z = np.nan_to_num(z, nan=0.0)
    return mat, samples, z


def reliability_table(stage3, genes, p_genes, targets, mat, samples, f1s_fit, f1_fit):
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
    stats = pd.read_csv(RES / "cptac_rna_gene_stats.tsv", sep="\t")
    cptac_mean = {}
    for cancer, sub in stats.groupby("cancer"):
        cptac_mean[cancer] = dict(zip(sub["gene_symbol"], sub["rna_mean"]))
    mu, sd = stage3.gene_mean_sd(mat.reindex(index=genes), samples)
    mu_pct = percentile_map(mu)
    sd_pct_ext = percentile_map(sd)
    f1_rows = []
    f1s_rows = []
    for gene in targets:
        vals = [rho[(c, gene)] for c in stage3.CANCERS if (c, gene) in rho]
        if vals:
            rho_mean = float(np.mean(vals))
            rho_sd = float(np.std(vals, ddof=0)) if len(vals) >= 2 else 0.0
            rho_min, rho_max, n_def = float(np.min(vals)), float(np.max(vals)), len(vals)
        else:
            rho_mean = gmean
            rho_sd, rho_min, rho_max, n_def = 0.0, gmean, gmean, 0
        means = [cptac_mean[c].get(gene, np.nan) for c in stage3.CANCERS]
        means = [v for v in means if np.isfinite(v)]
        rates_g = [rate.get((c, gene), 0.0) for c in stage3.CANCERS]
        hold_sd, hold_pct = sd_pct_ext.get(gene, (np.nan, np.nan))
        f1_rows.append(
            {
                "rho_mean": rho_mean,
                "rho_sd": rho_sd,
                "rho_min": rho_min,
                "rho_max": rho_max,
                "n_defined_cancers": n_def,
                "holdout_rna_mean": mu.get(gene, np.nan),
                "holdout_rna_sd": hold_sd,
                "holdout_rna_sd_percentile": hold_pct,
                "train_rna_mean": float(np.mean(means)) if means else np.nan,
                "train_protein_rate_mean": float(np.mean(rates_g)) if rates_g else np.nan,
                "gene_symbol": gene,
                "b_sd": hold_pct,
            }
        )
        ext_mu_pct = mu_pct.get(gene, (np.nan, np.nan))[1]
        f1s_rows.append(stage3.f1s_feature_row(gene, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, hold_pct) | {"gene_symbol": gene})
    f1_hat = stage3.apply_saved(pd.DataFrame(f1_rows), FEATURE_NAMES, f1_fit, ("f1_ridge", "f1_lgbm"))
    f1s_hat = stage3.apply_saved(pd.DataFrame(f1s_rows), F1S_FEATURES, f1s_fit, ("f1s_ridge", "f1s_lgbm"))
    out = pd.DataFrame({"gene_symbol": targets})
    out["f0"] = f1_hat["rho_mean"].to_numpy()
    out["f1_ridge"] = f1_hat["f1_ridge"].to_numpy()
    out["f1_lgbm"] = f1_hat["f1_lgbm"].to_numpy()
    out["f1s_ridge"] = f1s_hat["f1s_ridge"].to_numpy()
    out["f1s_lgbm"] = f1s_hat["f1s_lgbm"].to_numpy()
    out["b_sd"] = [row["b_sd"] for row in f1_rows]
    width_key = {"f0": "f0", "f1_ridge": "ridge", "f1_lgbm": "lgbm", "f1s_ridge": "ridge", "f1s_lgbm": "lgbm"}
    fit_of = {"f0": f1_fit, "f1_ridge": f1_fit, "f1_lgbm": f1_fit, "f1s_ridge": f1s_fit, "f1s_lgbm": f1s_fit}
    for name in ("f0", "f1_ridge", "f1_lgbm", "f1s_ridge", "f1s_lgbm"):
        for level in (80, 90):
            half = fit_of[name]["widths"][(width_key[name], level)]
            out[f"{name}_half_{level}"] = half
            out[f"{name}_low_{level}"] = out[name] - half
            out[f"{name}_high_{level}"] = out[name] + half
    return out


def similarity(scope, c09, genes, mat):
    profiles = scope.cptac_profiles(genes)
    block = mat.reindex(index=genes).to_numpy(dtype=np.float64)
    with np.errstate(all="ignore"):
        med = np.nanmedian(block, axis=1)
    mapped = percentile_map({g: float(v) if np.isfinite(v) else np.nan for g, v in zip(genes, med)})
    profile = {g: mapped[g][1] for g in genes}
    return scope.max_similarity(c09, profile, profiles, genes)


def main():
    print("thread env", {k: os.environ.get(k, "") for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    r1 = load_module("r1r2", ROOT / "scripts" / "12_stage5_r1r2.py")
    stage3 = load_module("b3stage", ROOT / "scripts" / "08_stage3.py")
    genes, p_genes, feats, weights, chosen, x, y_fill, alpha_r2 = fit_r1r2(r1)
    p_index = {g: i for i, g in enumerate(p_genes)}
    rows = r1.eval_rppa(genes, p_index, feats, weights, x, y_fill, alpha_r2, stage3)
    rows.extend(r1.eval_procan(genes, p_index, feats, weights, x, y_fill, alpha_r2))
    match = compare_stage5(rows)
    pd.DataFrame([match]).to_csv(OUT / "r1r2_stage5_match.tsv", sep="\t", index=False)
    np.savez_compressed(
        OUT / "r1r2_model.npz",
        weights=weights,
        feats=feats,
        alpha_r1=chosen,
        alpha_r2=np.array([alpha_r2]),
        x_tr=x,
        y_fill=y_fill,
        genes_g=np.array(genes),
        genes_p=np.array(p_genes),
    )
    f1s_fit, f1_fit = fit_reliability(stage3)
    mat, samples, z = apollo_z(stage3, genes)
    u_p1 = pd.read_csv(PREP / "U_P1_genes.tsv", sep="\t")["gene_symbol"].tolist()
    u_p2 = pd.read_csv(PREP / "U_P2_genes.tsv", sep="\t")["gene_symbol"].tolist()
    u_q = pd.read_csv(PREP / "U_Q_genes.tsv", sep="\t")["gene_symbol"].tolist()
    model, self_index, feat_index, coefs = stage3.models_from_parquet(RES / "ar_final_models.parquet", genes)
    symbol_at = list(model["gene_symbol"])
    pos = {g: i for i, g in enumerate(symbol_at)}
    keep = [pos[g] for g in u_p1 if g in pos]
    b = model["b"].to_numpy(dtype=np.float64)
    lams = model["lambda"].to_numpy(dtype=np.float64)
    ar = predict_ar_matrix(z, self_index[keep], feat_index[keep].T, b[keep], lams[keep], np.asarray(coefs[keep], dtype=np.float64).T)
    gindex = {g: i for i, g in enumerate(genes)}
    p1_genes = [symbol_at[i] for i in keep]
    rna_z = np.column_stack([z[:, gindex[g]] for g in p1_genes])
    q_index = [p_index[g] for g in u_q]
    r1_pred = np.column_stack([z[:, feats[:, j]] @ weights[:, j] for j in q_index])
    r2_all = r1.dual_predict(x, y_fill, z, alpha_r2)
    r2_pred = r2_all[:, q_index]
    np.savez_compressed(
        OUT / "apollo_sample_pred.npz",
        samples=np.array(samples),
        p1_genes=np.array(p1_genes),
        rna_z=rna_z,
        ar_pred=ar,
        q_genes=np.array(u_q),
        r1_pred=r1_pred,
        r2_pred=r2_pred,
    )
    table = reliability_table(stage3, genes, p_genes, u_p2, mat, samples, f1s_fit, f1_fit)
    table.to_csv(OUT / "apollo_gene_pred.tsv", sep="\t", index=False)
    scope = load_module("scope", ROOT / "scripts" / "12_stage5_scope.py")
    c09 = load_module("c09", ROOT / "scripts" / "09_confirm.py")
    s_value = similarity(scope, c09, genes, mat)
    pd.DataFrame([{"item": "s", "value": s_value}]).to_csv(OUT / "apollo_s.tsv", sep="\t", index=False)
    counts = pd.DataFrame(
        [
            {"item": "U_P1", "value": len(u_p1)},
            {"item": "n_ar_genes", "value": len(p1_genes)},
            {"item": "U_P2", "value": len(u_p2)},
            {"item": "n_gene_pred", "value": int(table["f1s_lgbm"].notna().sum())},
            {"item": "U_Q", "value": len(u_q)},
            {"item": "n_r1", "value": r1_pred.shape[1]},
            {"item": "n_samples", "value": len(samples)},
            {"item": "s", "value": s_value},
            {"item": "U_P2_absent_from_rna", "value": int(sum(g not in mat.index for g in u_p2))},
        ]
    )
    counts.to_csv(OUT / "prediction_counts.tsv", sep="\t", index=False)
    print(counts.to_string(index=False), flush=True)
    hash_rows = []
    for path in sorted(OUT.iterdir()):
        if path.is_file():
            hash_rows.append({"key": path.name, "path": str(path.relative_to(ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    hashes = pd.DataFrame(hash_rows)
    hashes.to_csv(ROOT / "config" / "frozen_hashes_B6.tsv", sep="\t", index=False)
    print("freeze done", flush=True)


if __name__ == "__main__":
    main()
