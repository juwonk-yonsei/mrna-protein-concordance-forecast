"""ADC target reliability atlas for TCGA primary tumors.

F1s is refit with the frozen procedure and the stored CPTAC half-widths.
The fit file from stage 3 is not rewritten.
"""
import importlib.util
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml
from lightgbm import LGBMRegressor
from sklearn.linear_model import Ridge

from scripts.b2_lib import RES, ROOT, conformal_half_width
from scripts.b3_lib import F1S_FEATURES, percentile_map

LABEL = "사후"
OUT = RES / "stage5" / "atlas"
POST = RES / "stage5" / "posthoc"
HPA_NOTE = "HPA는 암종 사이 수준 비교이고, ρ̂는 암종 안 환자 사이 일치도이다."
CPTAC_OF = {
    "breast invasive carcinoma": "BRCA",
    "kidney clear cell carcinoma": "CCRCC",
    "colon adenocarcinoma": "COAD",
    "glioblastoma multiforme": "GBM",
    "head & neck squamous cell carcinoma": "HNSCC",
    "lung squamous cell carcinoma": "LSCC",
    "lung adenocarcinoma": "LUAD",
    "ovarian serous cystadenocarcinoma": "OV",
    "pancreatic adenocarcinoma": "PDAC",
    "uterine corpus endometrioid carcinoma": "UCEC",
}
PROCAN_OF = {
    "lung adenocarcinoma": "Lung",
    "lung squamous cell carcinoma": "Lung",
    "colon adenocarcinoma": "Large Intestine",
    "rectum adenocarcinoma": "Large Intestine",
    "diffuse large B-cell lymphoma": "Haematopoietic and Lymphoid",
    "acute myeloid leukemia": "Haematopoietic and Lymphoid",
}


def load_stage3():
    spec = importlib.util.spec_from_file_location("b3stage", ROOT / "scripts" / "08_stage3.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    b2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
    mod.B2 = b2
    mod.B2_LGBM = dict(b2["reliability"]["lightgbm"])
    return mod


def write(path, frame):
    frame = frame.copy()
    frame.insert(0, "label", LABEL)
    frame.to_csv(path, sep="\t", index=False)
    print("wrote", path, len(frame), flush=True)


def finite_spearman(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < 3:
        return np.nan, n
    x = a[m]
    y = b[m]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, n
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if den == 0:
        return np.nan, n
    return float((rx * ry).sum() / den), n


def fit_f1s(stage3):
    preds = pd.read_csv(RES / "f1s_predictions.tsv", sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    train = preds[preds["in_P"] & np.isfinite(preds["actual"])].copy()
    stored = pd.read_csv(RES / "f1s_model_fit.tsv", sep="\t")
    stored = dict(zip(stored["item"], stored["value"]))
    alpha = 100.0
    if float(stored["ridge_alpha"]) != alpha:
        raise RuntimeError("stored F1s ridge alpha is not 100")
    x = train[F1S_FEATURES].to_numpy(dtype=np.float64)
    y = train["actual"].to_numpy(dtype=np.float64)
    x_f, _, _ = stage3.fill_median(x, x[:1])
    z, _ = stage3.standardize(x_f, x_f[:1])
    ridge = Ridge(alpha=alpha, fit_intercept=True)
    ridge.fit(z, y)
    lgbm = LGBMRegressor(**stage3.lgbm_params())
    lgbm.fit(x_f, y)
    widths = {}
    for name, col in (("f0", "f0"), ("ridge", "f1s_ridge"), ("lgbm", "f1s_lgbm")):
        resid = np.abs(train[col].to_numpy(dtype=np.float64) - y)
        for level, key in ((80, "q80"), (90, "q90")):
            half = conformal_half_width(resid, level)
            stored_half = float(stored[f"{key}_{name if name != 'ridge' else 'ridge'}" if False else f"{key}_{ {'f0':'f0','ridge':'ridge','lgbm':'lgbm'}[name] }"])
            if abs(half - stored_half) > 1e-8:
                raise RuntimeError(f"F1s width changed {name} {level}")
            widths[(name, level)] = half
    return {"alpha": alpha, "ridge": ridge, "lgbm": lgbm, "x_tr": x_f, "widths": widths}


def static_features(stage3, genes):
    rho = stage3.defined_rho(genes)
    p_genes = [g for g in (RES / "P_genes.txt").read_text().split() if g]
    adc = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    adc_loco = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    pset = set(p_genes)
    for rec in adc_loco.itertuples(index=False):
        if rec.gene_symbol in pset or rec.gene_symbol not in set(adc):
            continue
        if int(rec.n) >= 30 and np.isfinite(rec.s0):
            rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
    _stats, mean_pct, sd_pct, _mean = stage3.cptac_percentile_tables(genes)
    rates = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    rate = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
    gmean = stage3.global_mean(stage3.CANCERS, rho, p_genes)
    rows = {}
    for gene in genes:
        rows[gene] = stage3.f1s_feature_row(gene, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean, np.nan, np.nan)
    return rows


def chosen_tau():
    path = POST / "posthoc_22_tau_development.tsv"
    if not path.exists():
        return None
    frame = pd.read_csv(path, sep="\t")
    if str(frame["identified"].iloc[0]).lower() not in {"true", "1"}:
        return None
    value = frame["chosen_tau"].iloc[0]
    return float(value) if np.isfinite(value) else None


def hpa_table(targets, rna_pct, rho_hat):
    raw = pd.read_csv(ROOT / "data" / "raw" / "hpa" / "cancer_data.tsv", sep="\t")
    mapping = pd.read_csv(ROOT / "config" / "hpa_tcga_map.tsv", sep="\t", dtype=str).fillna("")
    mapped = mapping[mapping["mapped"] == "yes"]
    hpa_to_tcga = {str(k).casefold(): v for k, v in zip(mapped["hpa_cancer"], mapped["tcga_primary_disease"])}
    raw = raw[raw["Gene name"].isin(targets) & raw["Cancer"].astype(str).str.casefold().isin(hpa_to_tcga)].copy()
    for col in ("High", "Medium", "Low", "Not detected"):
        raw[col] = pd.to_numeric(raw[col], errors="coerce")
    raw["n_patients"] = raw[["High", "Medium", "Low", "Not detected"]].sum(axis=1)
    raw["high_medium_fraction"] = (raw["High"] + raw["Medium"]) / raw["n_patients"]
    raw.loc[raw["n_patients"] <= 0, "high_medium_fraction"] = np.nan
    raw["tcga_primary_disease"] = raw["Cancer"].astype(str).str.casefold().map(hpa_to_tcga)
    rows = []
    target_rho = []
    for gene in targets:
        sub = raw[raw["Gene name"] == gene]
        joined = []
        for rec in sub.itertuples(index=False):
            key = (rec.tcga_primary_disease, gene)
            if key not in rna_pct or not np.isfinite(rec.high_medium_fraction):
                continue
            joined.append((rna_pct[key], float(rec.high_medium_fraction)))
        n = len(joined)
        rho = np.nan
        if n >= 6:
            rho, _n = finite_spearman([a for a, _b in joined], [b for _a, b in joined])
        rows.append({"gene_symbol": gene, "n_cancers": n, "spearman_hpa_rna": rho if n >= 6 else np.nan})
        if n >= 6 and np.isfinite(rho) and gene in rho_hat and np.isfinite(rho_hat[gene]):
            target_rho.append((rho, rho_hat[gene]))
    detail = pd.DataFrame(rows)
    write(OUT / "hpa_by_target.tsv", detail)
    pooled_rho, pooled_n = (np.nan, 0)
    if len(target_rho) >= 3:
        pooled_rho, pooled_n = finite_spearman([a for a, _b in target_rho], [b for _a, b in target_rho])
    write(
        OUT / "hpa_summary.tsv",
        pd.DataFrame(
            [
                {
                    "n_targets_with_spearman": int(np.isfinite(detail["spearman_hpa_rna"]).sum()),
                    "n_targets": int(len(detail)),
                    "spearman_of_target_concordance_vs_mean_rho_hat": pooled_rho,
                    "n_targets_in_that_spearman": pooled_n,
                    "note": HPA_NOTE,
                    "map_commit": "9e5de1954a5251c9a278df1b34d623bb2c70c974",
                }
            ]
        ),
    )
    raw_out = raw[["Gene name", "Cancer", "tcga_primary_disease", "High", "Medium", "Low", "Not detected", "n_patients", "high_medium_fraction"]]
    write(OUT / "hpa_values.tsv", raw_out)


def figures(atlas, tau, scatter):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    diseases = list(dict.fromkeys(atlas["disease"]))
    targets = list(dict.fromkeys(atlas["gene_symbol"]))
    color = atlas.pivot(index="gene_symbol", columns="disease", values="rho_hat_percentile").reindex(index=targets, columns=diseases)
    size = atlas.pivot(index="gene_symbol", columns="disease", values="rna_median_percentile").reindex(index=targets, columns=diseases)
    scope = atlas.pivot(index="gene_symbol", columns="disease", values="out_of_scope").reindex(index=targets, columns=diseases)
    fig, ax = plt.subplots(figsize=(14, 8))
    xs, ys, cs, ss, edges = [], [], [], [], []
    for i, gene in enumerate(targets):
        for j, disease in enumerate(diseases):
            xs.append(j)
            ys.append(i)
            cs.append(color.loc[gene, disease])
            ss.append(size.loc[gene, disease])
            edges.append("black" if scope.loc[gene, disease] == 1 else "none")
    ss = np.asarray(ss, dtype=np.float64)
    point = 20 + 160 * (np.nan_to_num(ss, nan=0.0) / 100.0)
    ax.scatter(xs, ys, c=cs, s=point, cmap="viridis", edgecolors=edges, linewidths=0.6)
    labels = []
    blood = dict(zip(atlas["disease"], atlas["is_blood"]))
    for disease in diseases:
        labels.append(disease + (" *" if int(blood[disease]) == 1 else ""))
    ax.set_xticks(range(len(diseases)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(targets)))
    ax.set_yticklabels(targets, fontsize=8)
    ax.set_xlabel("TCGA cancer")
    ax.set_ylabel("target")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_adc_heatmap.png", dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for source, sub in scatter.groupby("source", sort=True):
        ax.scatter(sub["rho_hat"], sub["actual"], s=18, label=source)
    ax.set_xlabel("rho hat")
    ax.set_ylabel("actual rho")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_rho_hat_vs_actual.png", dpi=150)
    plt.close(fig)
    print("wrote figures", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("thread env", {k: os.environ.get(k, "") for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, flush=True)
    stage3 = load_stage3()
    genes = [g for g in (RES / "G_RNA_genes.txt").read_text().split() if g]
    targets = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    fit = fit_f1s(stage3)
    print("f1s refit", flush=True)
    features = static_features(stage3, genes)
    pct = pd.read_csv(OUT / "tcga_rna_percentiles.tsv", sep="\t")
    sim = pd.read_csv(POST / "posthoc_22_similarity.tsv", sep="\t")
    s_of = dict(zip(sim["context"], sim["s"]))
    tau = chosen_tau()
    draws = pd.read_csv(POST / "posthoc_24_draws.tsv", sep="\t")
    delta_rows = draws[(draws["platform"] == "rppa") & (draws["predictor"] == "f1s_lgbm") & (draws["K"] == 50)]
    platform_delta = float(delta_rows["delta"].median())
    model = pd.read_parquet(RES / "ar_final_models.parquet", columns=["gene_symbol", "lambda"])
    lam = dict(zip(model["gene_symbol"], model["lambda"]))
    ar_pairs = pd.read_csv(RES / "ar_pairs.tsv", sep="\t")
    loco_delta = ar_pairs.groupby("gene_symbol")["delta_ar"].mean().to_dict()
    rppa = pd.read_csv(RES / "stage4" / "confirm" / "rppa_pairs.tsv", sep="\t")
    pro = pd.read_csv(RES / "stage4" / "confirm" / "procan_pairs.tsv", sep="\t")
    f1s_pred = pd.read_csv(RES / "f1s_predictions.tsv", sep="\t")
    half = fit["widths"][("lgbm", 80)]
    atlas_rows = []
    rna_pct = {}
    mean_hat = {gene: [] for gene in targets}
    for disease, sub in pct.groupby("disease", sort=True):
        records = []
        order = []
        for rec in sub.itertuples(index=False):
            base = dict(features[rec.gene_symbol])
            base["holdout_rna_mean_percentile"] = rec.median_percentile
            base["holdout_rna_sd_percentile"] = rec.sd_percentile
            base["gene_symbol"] = rec.gene_symbol
            records.append(base)
            order.append(rec.gene_symbol)
            rna_pct[(disease, rec.gene_symbol)] = rec.median_percentile
        frame = pd.DataFrame(records)
        scored = stage3.apply_saved(frame, F1S_FEATURES, fit, ("f1s_ridge", "f1s_lgbm"))
        hat = scored["f1s_lgbm"].to_numpy(dtype=np.float64)
        ranked = percentile_map({g: float(v) if np.isfinite(v) else np.nan for g, v in zip(order, hat)})
        lookup = dict(zip(order, hat))
        s = s_of.get(disease, np.nan)
        blood = int(sub["is_blood"].iloc[0])
        n_samples = int(sub["n_samples"].iloc[0])
        for gene in targets:
            if gene not in lookup:
                continue
            value = float(lookup[gene])
            mean_hat[gene].append(value)
            rppa_sub = rppa[(rppa["cancer"] == disease) & (rppa["gene_symbol"] == gene)]
            rppa_actual = float(rppa_sub["actual"].mean()) if len(rppa_sub) and np.isfinite(rppa_sub["actual"]).any() else np.nan
            lin = PROCAN_OF.get(disease)
            pro_sub = pro[(pro["lineage"] == lin) & (pro["gene_symbol"] == gene)] if lin else pro.iloc[0:0]
            pro_actual = float(pro_sub["actual"].mean()) if len(pro_sub) and np.isfinite(pro_sub["actual"]).any() else np.nan
            cptac = CPTAC_OF.get(disease)
            cptac_sub = ar_pairs[(ar_pairs["heldout_cancer"] == cptac) & (ar_pairs["gene_symbol"] == gene)] if cptac else ar_pairs.iloc[0:0]
            cptac_actual = float(cptac_sub["s0"].mean()) if len(cptac_sub) and np.isfinite(cptac_sub["s0"]).any() else np.nan
            gene_pct = pct[(pct["disease"] == disease) & (pct["gene_symbol"] == gene)]
            atlas_rows.append(
                {
                    "gene_symbol": gene,
                    "disease": disease,
                    "n_samples": n_samples,
                    "is_blood": blood,
                    "rna_median_percentile": float(gene_pct["median_percentile"].iloc[0]) if len(gene_pct) else np.nan,
                    "rna_sd_percentile": float(gene_pct["sd_percentile"].iloc[0]) if len(gene_pct) else np.nan,
                    "rho_hat": value,
                    "rho_hat_low_80": value - half,
                    "rho_hat_high_80": value + half,
                    "rho_hat_percentile": ranked[gene][1],
                    "rho_hat_rppa_shift": value + platform_delta,
                    "rppa_delta_k50_median": platform_delta,
                    "in_ar": int(gene in lam),
                    "lambda": lam.get(gene, np.nan),
                    "cptac_loco_mean_delta": loco_delta.get(gene, np.nan),
                    "s": s,
                    "tau": tau if tau is not None else np.nan,
                    "out_of_scope": int(tau is not None and np.isfinite(s) and s < tau),
                    "rppa_actual": rppa_actual,
                    "procan_lineage": lin if lin else "",
                    "procan_actual": pro_actual,
                    "cptac_cancer": cptac if cptac else "",
                    "cptac_loco_actual": cptac_actual,
                }
            )
        print("atlas", disease, flush=True)
    atlas = pd.DataFrame(atlas_rows)
    write(OUT / "adc_atlas.tsv", atlas)
    scatter_rows = []
    loco = f1s_pred[f1s_pred["gene_symbol"].isin(targets)]
    for rec in loco.itertuples(index=False):
        if np.isfinite(rec.f1s_lgbm) and np.isfinite(rec.actual):
            scatter_rows.append({"source": "CPTAC", "gene_symbol": rec.gene_symbol, "context": rec.heldout_cancer, "rho_hat": float(rec.f1s_lgbm), "actual": float(rec.actual)})
    for rec in atlas.itertuples(index=False):
        if np.isfinite(rec.rppa_actual):
            scatter_rows.append({"source": "RPPA", "gene_symbol": rec.gene_symbol, "context": rec.disease, "rho_hat": rec.rho_hat, "actual": rec.rppa_actual})
        if np.isfinite(rec.procan_actual):
            scatter_rows.append({"source": "ProCan", "gene_symbol": rec.gene_symbol, "context": rec.disease, "rho_hat": rec.rho_hat, "actual": rec.procan_actual})
    # ProCan frozen forecasts, separate from the TCGA-row actuals above, are the external predictions.
    pro_f = pd.read_csv(RES / "frozen_forecasts_procan_F1s.tsv", sep="\t")
    pro_obs = pro.merge(pro_f, on=["protein_id", "gene_symbol", "lineage"], suffixes=("", "_f"))
    pro_obs = pro_obs[pro_obs["gene_symbol"].isin(targets)]
    frozen_rows = []
    for rec in pro_obs.itertuples(index=False):
        if np.isfinite(rec.actual) and np.isfinite(rec.f1s_lgbm):
            frozen_rows.append({"source": "ProCan_frozen_forecast", "gene_symbol": rec.gene_symbol, "context": rec.lineage, "rho_hat": float(rec.f1s_lgbm), "actual": float(rec.actual)})
    scatter = pd.DataFrame(scatter_rows + frozen_rows)
    write(OUT / "concordance_pairs.tsv", scatter)
    summary = []
    for source, sub in scatter.groupby("source", sort=True):
        rho, n = finite_spearman(sub["rho_hat"], sub["actual"])
        summary.append({"source": source, "level": "pooled", "gene_symbol": "", "n_pairs": n, "spearman": rho})
        for gene, gsub in sub.groupby("gene_symbol", sort=True):
            rho_g, n_g = finite_spearman(gsub["rho_hat"], gsub["actual"])
            summary.append({"source": source, "level": "gene", "gene_symbol": gene, "n_pairs": n_g, "spearman": rho_g if n_g >= 3 else np.nan})
    write(OUT / "concordance_summary.tsv", pd.DataFrame(summary))
    rho_mean = {gene: float(np.nanmean(vals)) if vals else np.nan for gene, vals in mean_hat.items()}
    hpa_table(targets, rna_pct, rho_mean)
    figures(atlas, tau, scatter[scatter["source"].isin(["CPTAC", "RPPA", "ProCan"])])
    print("atlas done", flush=True)


if __name__ == "__main__":
    main()
