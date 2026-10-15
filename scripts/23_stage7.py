"""Post-hoc F0 comparison, calibration table, and ADC atlas final table.

Does not refit models and does not modify frozen predictions or confirmation outputs.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.common import ROOT

LABEL = "사후"
OUT = ROOT / "results" / "stage7"
SEED = 20261001
DRAWS = 2000
GMEAN_TOL = 1e-8


def average_ranks(x):
    x = np.asarray(x, dtype=np.float64)
    n = int(x.size)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    distinct = np.empty(n, dtype=bool)
    distinct[0] = True
    if n > 1:
        distinct[1:] = xs[1:] != xs[:-1]
    groups = np.cumsum(distinct) - 1
    pos = np.arange(1, n + 1, dtype=np.float64)
    counts = np.bincount(groups)
    ranks = np.empty(n, dtype=np.float64)
    ranks[order] = (np.bincount(groups, weights=pos) / counts)[groups]
    return ranks


def finite_spearman(a, b, min_n=3):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < min_n:
        return np.nan
    x = a[m]
    y = b[m]
    if x.max() == x.min() or y.max() == y.min():
        return np.nan
    rx = average_ranks(x)
    ry = average_ranks(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if den == 0:
        return np.nan
    return float((rx * ry).sum() / den)


def auroc_high(score, positive):
    score = np.asarray(score, dtype=np.float64)
    positive = np.asarray(positive, dtype=bool)
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = average_ranks(score)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def write(name, frame):
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / name, sep="\t", index=False)


def global_mean():
    s0 = pd.read_csv(ROOT / "results" / "s0_pairs.tsv", sep="\t")
    p_genes = set(pd.read_csv(ROOT / "results" / "P_genes.txt", sep="\t", header=None)[0].astype(str))
    use = s0[(s0["n"] >= 30) & np.isfinite(s0["s0"]) & s0["gene_symbol"].isin(p_genes)]
    return float(use["s0"].mean()), int(use.shape[0]), len(p_genes)


def load_frames(p_genes, gmean):
    rppa = pd.read_csv(
        ROOT / "results" / "stage4" / "confirm" / "rppa_pairs.tsv",
        sep="\t",
        usecols=["xena_antibody", "gene_symbol", "cancer", "actual", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"],
    )
    rppa["unit"] = rppa["xena_antibody"].astype(str)
    rppa["context"] = rppa["cancer"].astype(str)
    pro = pd.read_csv(
        ROOT / "results" / "stage4" / "confirm" / "procan_pairs.tsv",
        sep="\t",
        usecols=["protein_id", "gene_symbol", "lineage", "actual", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"],
    )
    pro["unit"] = pro["protein_id"].astype(str)
    pro["context"] = pro["lineage"].astype(str)
    genes = pd.read_csv(ROOT / "results" / "stage6" / "confirm" / "genes.tsv", sep="\t", usecols=["gene_symbol", "in_U_P2", "actual"])
    pred = pd.read_csv(
        ROOT / "results" / "stage6" / "frozen" / "apollo_gene_pred.tsv",
        sep="\t",
        usecols=["gene_symbol", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"],
    )
    apollo = genes.merge(pred, on="gene_symbol", how="left")
    apollo = apollo[apollo["in_U_P2"] == 1].copy()
    apollo["unit"] = apollo["gene_symbol"].astype(str)
    apollo["context"] = "LUAD"
    frames = {"RPPA": rppa, "ProCan": pro, "APOLLO": apollo}
    rows = []
    for name, frame in frames.items():
        gene_f0 = frame.groupby("gene_symbol", sort=True)["f0"].mean()
        for gene, value in gene_f0.items():
            in_p = gene in p_genes
            if in_p:
                fill = "in_P_mean_of_defined_CPTAC_rho"
            elif np.isfinite(value) and abs(float(value) - gmean) <= GMEAN_TOL:
                fill = "outside_P_global_mean"
            else:
                fill = "outside_P_gene_specific"
            rows.append({"label": LABEL, "dataset": name, "gene_symbol": gene, "in_P": int(in_p), "f0": float(value) if np.isfinite(value) else np.nan, "fill": fill})
    return frames, pd.DataFrame(rows)


def draw_indices(groups, rng):
    chosen = rng.integers(0, len(groups), size=len(groups))
    return np.concatenate([groups[i] for i in chosen])


def bootstrap_diff(actual, left, right, groups):
    point_l = finite_spearman(left, actual, 3)
    point_r = finite_spearman(right, actual, 3)
    rng = np.random.default_rng(SEED)
    diffs = np.empty(DRAWS)
    for i in range(DRAWS):
        idx = draw_indices(groups, rng)
        diffs[i] = finite_spearman(left[idx], actual[idx], 3) - finite_spearman(right[idx], actual[idx], 3)
    ok = diffs[np.isfinite(diffs)]
    return {
        "s_model": point_l,
        "s_base": point_r,
        "diff": point_l - point_r if np.isfinite(point_l) and np.isfinite(point_r) else np.nan,
        "ci_low": float(np.quantile(ok, 0.025)) if ok.size else np.nan,
        "ci_high": float(np.quantile(ok, 0.975)) if ok.size else np.nan,
        "n_nonfinite_draws": int(DRAWS - ok.size),
        "n_units": int(len(groups)),
        "n_pairs": int(np.isfinite(actual).sum()),
    }


def bootstrap_auroc(actual, score, groups):
    point = auroc_high(score, actual >= 0.5)
    rng = np.random.default_rng(SEED)
    draws = np.empty(DRAWS)
    for i in range(DRAWS):
        idx = draw_indices(groups, rng)
        draws[i] = auroc_high(score[idx], actual[idx] >= 0.5)
    ok = draws[np.isfinite(draws)]
    return {
        "auroc": point,
        "ci_low": float(np.quantile(ok, 0.025)) if ok.size else np.nan,
        "ci_high": float(np.quantile(ok, 0.975)) if ok.size else np.nan,
        "n_nonfinite_draws": int(DRAWS - ok.size),
        "n_units": int(len(groups)),
        "n_pairs": int(np.isfinite(actual).sum()),
        "n_positive": int(np.sum(np.isfinite(actual) & np.isfinite(score) & (actual >= 0.5))),
    }


def f0_comparisons(frames, p_genes):
    diff_rows = []
    auroc_rows = []
    context_rows = []
    median_rows = []
    predictors = ("f1s_lgbm", "f1s_ridge")
    for dataset, frame in frames.items():
        frame = frame.reset_index(drop=True)
        in_p = frame["gene_symbol"].isin(p_genes).to_numpy()
        strata = {"all": np.ones(len(frame), dtype=bool), "in_P": in_p, "outside_P": ~in_p}
        actual = frame["actual"].to_numpy(dtype=np.float64)
        f0 = frame["f0"].to_numpy(dtype=np.float64)
        for pred in predictors:
            model = frame[pred].to_numpy(dtype=np.float64)
            for stratum, mask0 in strata.items():
                mask = mask0 & np.isfinite(actual) & np.isfinite(f0) & np.isfinite(model)
                sub = frame.loc[mask].reset_index(drop=True)
                groups = [np.asarray(v, dtype=np.int64) for v in sub.groupby("unit", sort=True).indices.values()]
                if not groups:
                    continue
                rec = bootstrap_diff(sub["actual"].to_numpy(dtype=np.float64), sub[pred].to_numpy(dtype=np.float64), sub["f0"].to_numpy(dtype=np.float64), groups)
                rec.update({"label": LABEL, "dataset": dataset, "predictor": pred, "stratum": stratum, "base": "f0"})
                diff_rows.append(rec)
                print(dataset, pred, stratum, rec["diff"], flush=True)
        for pred in ("f0", "f1s_lgbm", "f1s_ridge", "b_sd"):
            model = frame[pred].to_numpy(dtype=np.float64)
            for stratum, mask0 in strata.items():
                mask = mask0 & np.isfinite(actual) & np.isfinite(model)
                sub = frame.loc[mask].reset_index(drop=True)
                groups = [np.asarray(v, dtype=np.int64) for v in sub.groupby("unit", sort=True).indices.values()]
                if not groups:
                    continue
                rec = bootstrap_auroc(sub["actual"].to_numpy(dtype=np.float64), sub[pred].to_numpy(dtype=np.float64), groups)
                rec.update({"label": LABEL, "dataset": dataset, "predictor": pred, "stratum": stratum})
                auroc_rows.append(rec)
        for pred in ("f0", "f1s_lgbm", "f1s_ridge", "b_sd"):
            vals = []
            for context, sub in frame.groupby("context", sort=True):
                rho = finite_spearman(sub[pred], sub["actual"], 3)
                n = int((np.isfinite(sub[pred]) & np.isfinite(sub["actual"])).sum())
                context_rows.append({"label": LABEL, "dataset": dataset, "context": context, "predictor": pred, "n_pairs": n, "spearman": rho})
                if np.isfinite(rho):
                    vals.append(rho)
            median_rows.append(
                {
                    "label": LABEL,
                    "dataset": dataset,
                    "predictor": pred,
                    "n_contexts_defined": len(vals),
                    "median_spearman": float(np.median(vals)) if vals else np.nan,
                }
            )
    write("f0_diff.tsv", pd.DataFrame(diff_rows))
    write("f0_auroc.tsv", pd.DataFrame(auroc_rows))
    write("f0_context.tsv", pd.DataFrame(context_rows))
    write("f0_context_median.tsv", pd.DataFrame(median_rows))


def within_antibody(rppa):
    rows = []
    for antibody, sub in rppa.groupby("xena_antibody", sort=True):
        for pred in ("f1s_lgbm", "f1s_ridge", "b_sd", "f0"):
            m = np.isfinite(sub[pred]) & np.isfinite(sub["actual"])
            n = int(m.sum())
            constant = n > 0 and (sub.loc[m, pred].nunique(dropna=True) < 2)
            rho = np.nan if n < 10 or constant else finite_spearman(sub.loc[m, pred], sub.loc[m, "actual"], 3)
            rows.append(
                {
                    "label": LABEL,
                    "xena_antibody": antibody,
                    "gene_symbol": sub["gene_symbol"].iloc[0],
                    "predictor": pred,
                    "n_cancers": n,
                    "constant_predictor": int(constant),
                    "spearman": rho,
                }
            )
    detail = pd.DataFrame(rows)
    write("f0_within_antibody.tsv", detail)
    stored = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "posthoc_25_diff.tsv", sep="\t")
    summary = []
    for pred, sub in detail.groupby("predictor", sort=True):
        defined = sub[np.isfinite(sub["spearman"])]
        summary.append(
            {
                "label": LABEL,
                "predictor": pred,
                "n_antibodies": int(sub.shape[0]),
                "n_defined": int(defined.shape[0]),
                "median_spearman": float(defined["spearman"].median()) if len(defined) else np.nan,
                "note": "F0 is constant within an antibody, so its within-antibody Spearman is undefined." if pred == "f0" else "",
            }
        )
    for rec in stored.itertuples(index=False):
        summary.append(
            {
                "label": LABEL,
                "predictor": f"{rec.predictor}_minus_b_sd",
                "n_antibodies": int(rec.n_antibodies),
                "n_defined": int(rec.n_antibodies),
                "median_spearman": float(rec.median_diff),
                "ci_low": float(rec.ci_low),
                "ci_high": float(rec.ci_high),
                "bootstrap_draws": int(rec.bootstrap_draws),
                "bootstrap_seed": int(rec.bootstrap_seed),
                "note": "Stored B5 antibody-resample interval of the median within-antibody Spearman difference. Not recomputed.",
            }
        )
    write("f0_within_antibody_summary.tsv", pd.DataFrame(summary))


def calibration():
    b5 = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "posthoc_24_draws.tsv", sep="\t")
    b6 = pd.read_csv(ROOT / "results" / "stage6" / "confirm" / "descriptive_calibration_draws.tsv", sep="\t")
    b6 = b6.copy()
    b6.insert(0, "platform", "APOLLO")
    b6.insert(0, "label", LABEL)
    keep = [
        "platform",
        "predictor",
        "K",
        "cov80_none",
        "cov90_none",
        "cov80_bias",
        "cov90_bias",
        "cov80_both",
        "cov90_both",
        "width80_none",
        "width90_none",
        "width80_bias",
        "width90_bias",
        "width80_both",
        "width90_both",
    ]
    draws = pd.concat([b5[keep], b6[keep]], ignore_index=True)
    methods = {
        "none": ("cov80_none", "cov90_none", "width80_none", "width90_none"),
        "delta": ("cov80_bias", "cov90_bias", "width80_bias", "width90_bias"),
        "delta_quantile": ("cov80_both", "cov90_both", "width80_both", "width90_both"),
    }
    rows = []
    for (platform, predictor, k), sub in draws.groupby(["platform", "predictor", "K"], sort=True):
        for method, cols in methods.items():
            rec = {
                "label": LABEL,
                "platform": platform,
                "predictor": predictor,
                "K": int(k),
                "method": method,
                "n_draws": int(len(sub)),
            }
            names = ("cov80", "cov90", "width80", "width90")
            for name, col in zip(names, cols):
                values = sub[col].to_numpy(dtype=np.float64)
                rec[f"{name}_mean"] = float(np.mean(values))
                rec[f"{name}_p2_5"] = float(np.quantile(values, 0.025))
                rec[f"{name}_p97_5"] = float(np.quantile(values, 0.975))
            rows.append(rec)
    table = pd.DataFrame(rows)
    write("calibration_table.tsv", table)
    primary = table[table["predictor"] == "f1s_lgbm"]
    recs = []
    for (method, k), sub in primary.groupby(["method", "K"], sort=True):
        ok = (sub["cov80_mean"] >= 0.75) & (sub["cov80_mean"] <= 0.85)
        recs.append(
            {
                "label": LABEL,
                "predictor": "f1s_lgbm",
                "method": method,
                "K": int(k),
                "n_platforms_in_0.75_0.85": int(ok.sum()),
                "all_three": int(ok.sum() == 3 and sub["platform"].nunique() == 3),
                "platforms": ",".join(sub.loc[ok, "platform"].astype(str)),
            }
        )
    write("calibration_recommend.tsv", pd.DataFrame(recs))


def atlas(approved):
    src = pd.read_csv(ROOT / "results" / "stage5" / "atlas" / "adc_atlas_v2.tsv", sep="\t")
    notes = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t")
    src = src.merge(notes, on="gene_symbol", how="left")
    src = src.merge(approved, on="gene_symbol", how="left")
    src["label"] = LABEL
    write("adc_atlas_final.tsv", src)
    columns = [
        ("gene_symbol", "HGNC symbol of the ADC target."),
        ("disease", "TCGA primary-disease name."),
        ("n_samples", "Primary-tumor RNA samples used for the RNA percentiles in that disease."),
        ("is_blood", "1 if the disease name contains leukemia, lymphoma, myeloma, haematopoietic, or hematopoietic."),
        ("rna_median_percentile", "Percentile of the target RNA median among G_RNA genes in that disease."),
        ("rna_sd_percentile", "Percentile of the target RNA standard deviation among G_RNA genes in that disease."),
        ("rho_hat", "Frozen F1s-LightGBM prediction of within-context Spearman correlation."),
        ("rho_hat_low_80", "rho_hat minus the stored CPTAC 80% half-width. Not clipped to [-1, 1]."),
        ("rho_hat_high_80", "rho_hat plus the stored CPTAC 80% half-width. Not clipped to [-1, 1]."),
        ("rho_hat_percentile", "Percentile of rho_hat among G_RNA genes in that disease."),
        ("rho_hat_rppa_shift", "rho_hat plus the stage-5 K=50 RPPA delta median. Display only."),
        ("rppa_delta_k50_median", "Stage-5 median of the K=50 RPPA calibration delta."),
        ("in_ar", "1 if the gene is in the anchor-residual set."),
        ("lambda", "Selected anchor weight on the final CPTAC model."),
        ("cptac_loco_mean_delta", "Mean CPTAC leave-one-cancer-out delta for the gene."),
        ("s", "RNA-profile similarity of the disease to the nearest CPTAC tumor profile."),
        ("tau", "Unused threshold. Empty in this atlas."),
        ("out_of_scope", "1 if a threshold tau was set and s < tau. tau is empty, so this is 0."),
        ("rppa_actual", "Observed RPPA Spearman correlation when that disease and gene were measured."),
        ("procan_lineage", "ProCan lineage mapped to the disease, when a map exists."),
        ("procan_actual", "Observed ProCan Spearman correlation for that lineage and gene."),
        ("cptac_cancer", "CPTAC cancer abbreviation mapped to the disease, when the map is one-to-one."),
        ("cptac_loco_actual", "CPTAC leave-one-cancer-out observed Spearman correlation."),
        ("apollo_actual", "Observed APOLLO Spearman correlation. Filled only on lung adenocarcinoma rows."),
        ("apollo_d", "APOLLO rho_AR minus rho_S0. Filled only on lung adenocarcinoma rows."),
        ("apollo_rho_hat", "APOLLO-sample F1s-LightGBM prediction. Filled only on lung adenocarcinoma rows."),
        ("isoform_note", "From config/adc_targets.tsv. CLDN18 states that gene-level assays do not separate CLDN18.2."),
        ("approved_adc", "FDA-label antibody-drug conjugate naming this target. Empty if not confirmed from an FDA label in this step."),
        ("approval_url", "FDA label URL used for approved_adc. Empty if not confirmed."),
    ]
    write("adc_columns.tsv", pd.DataFrame([{"label": LABEL, "column": name, "definition": text} for name, text in columns]))
    conc = pd.read_csv(ROOT / "results" / "stage6" / "posthoc" / "posthoc_04_concordance.tsv", sep="\t")
    conc["label"] = LABEL
    write("adc_concordance.tsv", conc)


def main():
    gmean, n_pairs, n_p = global_mean()
    p_genes = set(pd.read_csv(ROOT / "results" / "P_genes.txt", sep="\t", header=None)[0].astype(str))
    frames, fill = load_frames(p_genes, gmean)
    write("f0_fill.tsv", fill)
    counts = fill.groupby(["dataset", "fill"], sort=True).size().rename("n_genes").reset_index()
    counts.insert(0, "label", LABEL)
    counts["global_mean"] = gmean
    counts["n_s0_pairs_in_P"] = n_pairs
    counts["n_P_genes"] = n_p
    write("f0_fill_counts.tsv", counts)
    print("fill counts written", flush=True)
    f0_comparisons(frames, p_genes)
    within_antibody(frames["RPPA"])
    calibration()
    approved = pd.DataFrame(
        [
            {
                "gene_symbol": "ERBB2",
                "approved_adc": "fam-trastuzumab deruxtecan-nxki (Enhertu), HER2-directed antibody-drug conjugate",
                "approval_url": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/761139s028lbl.pdf",
            },
            {
                "gene_symbol": "TACSTD2",
                "approved_adc": "sacituzumab govitecan-hziy (Trodelvy), Trop-2-directed antibody-drug conjugate",
                "approval_url": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/761115s059lbl.pdf",
            },
            {
                "gene_symbol": "NECTIN4",
                "approved_adc": "enfortumab vedotin-ejfv (Padcev), Nectin-4-directed antibody-drug conjugate",
                "approval_url": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/761137s029lbl.pdf",
            },
            {
                "gene_symbol": "F3",
                "approved_adc": "tisotumab vedotin-tftv (Tivdak), tissue factor-directed antibody-drug conjugate",
                "approval_url": "https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/761208s005lbl.pdf",
            },
        ]
    )
    atlas(approved)
    print("stage7 done", flush=True)


if __name__ == "__main__":
    main()
