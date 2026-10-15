"""Post-hoc tables. Not part of the preregistered confirmation.

Does not read locked values, config/paths.yaml, or thread-cap settings.
"""
import os

os.environ["B_CPU_CAP"] = "none"

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
CONFIRM = ROOT / "results" / "stage4" / "confirm"
OUT = ROOT / "results" / "stage4" / "posthoc"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261001
BOOT_N = 2000
LABEL = "사후"


def write(name, frame):
    out = frame.copy()
    out.insert(0, "label", LABEL)
    out.to_csv(OUT / name, sep="\t", index=False)


def effect(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    rec = {
        "n": int(values.size),
        "n_zero": int(np.sum(values == 0)) if values.size else 0,
        "statistic": np.nan,
        "pvalue": np.nan,
        "mean": np.nan,
        "mean_ci_low": np.nan,
        "mean_ci_high": np.nan,
        "median": np.nan,
        "median_ci_low": np.nan,
        "median_ci_high": np.nan,
        "error": "",
    }
    if values.size == 0:
        rec["error"] = "n=0"
        return rec
    rec["mean"] = float(values.mean())
    rec["median"] = float(np.median(values))
    try:
        result = wilcoxon(values, alternative="greater", zero_method="wilcox")
        rec["statistic"] = float(result.statistic)
        rec["pvalue"] = float(result.pvalue)
    except ValueError as exc:
        rec["error"] = f"ValueError: {exc}"
    rng = np.random.default_rng(SEED)
    drawn = rng.choice(values, size=(BOOT_N, values.size), replace=True)
    means = drawn.mean(axis=1)
    medians = np.median(drawn, axis=1)
    rec["mean_ci_low"] = float(np.quantile(means, 0.025))
    rec["mean_ci_high"] = float(np.quantile(means, 0.975))
    rec["median_ci_low"] = float(np.quantile(medians, 0.025))
    rec["median_ci_high"] = float(np.quantile(medians, 0.975))
    return rec


def main():
    units = pd.read_csv(CONFIRM / "h1_units.tsv", sep="\t")
    models = pd.read_parquet(ROOT / "results" / "ar_final_models.parquet", columns=["gene_symbol", "lambda", "alpha"])
    missing = pd.read_csv(ROOT / "results" / "ar_extra_missing.tsv", sep="\t")
    missing = missing[(missing["cohort"] == "TCGA") & (missing["context"] == "all_diseases")]
    table = units.merge(models, on="gene_symbol", how="left").merge(
        missing[["gene_symbol", "n_extra_missing"]], on="gene_symbol", how="left"
    )
    table = table.rename(columns={"mean_delta": "d_a", "n_contexts": "n_defined_cancers", "n_extra_missing": "n_residual_genes_missing_tcga"})
    table = table.sort_values(["d_a", "xena_antibody"], ascending=[False, True], kind="mergesort")
    write("posthoc_01_rppa_antibodies.tsv", table)
    write("posthoc_01_rppa_top15.tsv", table.head(15).copy())
    write("posthoc_01_rppa_bottom15.tsv", table.tail(15).sort_values(["d_a", "xena_antibody"], ascending=[True, True], kind="mergesort"))

    p_genes = set(g for g in (ROOT / "results" / "P_genes.txt").read_text().split() if g)
    h2 = pd.read_csv(CONFIRM / "h2_units.tsv", sep="\t")
    lam_rows = []
    rppa_lam = table.dropna(subset=["lambda"])
    for lam, sub in rppa_lam.groupby("lambda", sort=True):
        d = sub["d_a"].to_numpy(dtype=np.float64)
        lam_rows.append({"cohort": "rppa", "lambda": float(lam), "n": int(np.isfinite(d).sum()), "median_d": float(np.nanmedian(d))})
    h2j = h2.merge(models, on="gene_symbol", how="left")
    for lam, sub in h2j.dropna(subset=["lambda"]).groupby("lambda", sort=True):
        d = sub["mean_delta"].to_numpy(dtype=np.float64)
        lam_rows.append({"cohort": "procan", "lambda": float(lam), "n": int(np.isfinite(d).sum()), "median_d": float(np.nanmedian(d))})
    write("posthoc_02_lambda_median.tsv", pd.DataFrame(lam_rows))

    procan = pd.read_csv(CONFIRM / "procan_pairs.tsv", sep="\t")
    proteins = set(h2["protein_id"])
    procan = procan[procan["protein_id"].isin(proteins)].copy()
    lin_rows = []
    for lin, sub in procan.groupby("lineage", sort=True):
        d = sub["delta"].to_numpy(dtype=np.float64)
        finite = d[np.isfinite(d)]
        rec = effect(finite)
        rec.update(
            {
                "lineage": lin,
                "q1": float(np.quantile(finite, 0.25)) if finite.size else np.nan,
                "q3": float(np.quantile(finite, 0.75)) if finite.size else np.nan,
                "frac_gt_0.05": float(np.mean(finite > 0.05)) if finite.size else np.nan,
                "frac_lt_m0.05": float(np.mean(finite < -0.05)) if finite.size else np.nan,
            }
        )
        lin_rows.append(rec)
    write("posthoc_03_lineage.tsv", pd.DataFrame(lin_rows))

    err_rows = []
    rppa = pd.read_csv(CONFIRM / "rppa_pairs.tsv", sep="\t")
    for cohort, frame, id_col, context_col in (
        ("rppa", rppa, "xena_antibody", "cancer"),
        ("procan", pd.read_csv(CONFIRM / "procan_pairs.tsv", sep="\t"), "protein_id", "lineage"),
    ):
        use = frame[np.isfinite(frame["actual"]) & np.isfinite(frame["f1s_lgbm"])].copy()
        use["abs_error"] = (use["actual"] - use["f1s_lgbm"]).abs()
        use = use.sort_values(["abs_error", id_col, context_col], ascending=[False, True, True], kind="mergesort").head(20)
        for rank, rec in enumerate(use.itertuples(index=False), start=1):
            err_rows.append(
                {
                    "cohort": cohort,
                    "rank": rank,
                    "id": getattr(rec, id_col),
                    "gene_symbol": rec.gene_symbol,
                    "context": getattr(rec, context_col),
                    "actual": rec.actual,
                    "f1s_lgbm": rec.f1s_lgbm,
                    "abs_error": rec.abs_error,
                }
            )
    write("posthoc_04_abs_error_top20.tsv", pd.DataFrame(err_rows))

    mem_rows = []
    for cohort, frame, value_col in (("rppa", table, "d_a"), ("procan", h2, "mean_delta")):
        for name, mask in (("in_P", frame["gene_symbol"].isin(p_genes)), ("outside_P", ~frame["gene_symbol"].isin(p_genes))):
            d = frame.loc[mask, value_col].to_numpy(dtype=np.float64)
            d = d[np.isfinite(d)]
            mem_rows.append({"cohort": cohort, "group": name, "n": int(d.size), "median_d": float(np.median(d)) if d.size else np.nan})
    write("posthoc_05_p_membership.tsv", pd.DataFrame(mem_rows))
    print("posthoc written", OUT, flush=True)


if __name__ == "__main__":
    main()
