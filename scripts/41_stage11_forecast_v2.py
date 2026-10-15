"""Add observed CPTAC lookup columns to a copy of the all-gene forecast table.

The plan file config/analysis_plan_B11.yaml was committed before this script runs.
The stage-10 table is read and not rewritten. No model is fit.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import importlib

import numpy as np
import pandas as pd
import yaml

from scripts.b2_lib import CANCERS
from scripts.common import ROOT

lookup = importlib.import_module("scripts.25_stage8_lookup")

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B11.yaml"))["forecast_v2"]
OUT = ROOT / "results" / "stage11"


def values_by_gene(rho):
    grouped = {gene: [] for gene in []}
    for (cancer, gene), value in rho.items():
        if cancer not in CANCERS or not np.isfinite(value):
            continue
        grouped.setdefault(gene, []).append(float(value))
    return grouped


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(ROOT / PLAN["source"], sep="\t")
    assignment = pd.read_csv(ROOT / "results" / "stage8" / "lookup_assignment.tsv", sep="\t")
    names = set(pd.read_csv(ROOT / PLAN["cancer_names"], sep="\t")["cancer_type"].astype(str))
    unknown = sorted(set(source["cancer_type"].astype(str)) - names)
    if unknown:
        raise SystemExit(f"cancer names absent from s_by_cancer.tsv: {unknown}")
    eval_genes = sorted(set(assignment["gene_symbol"].astype(str)))
    print("stage-8 rho dictionary", len(eval_genes), flush=True)
    rho, l2, _gmean, n_added, n_g = lookup.lookup_rho(eval_genes)
    grouped = values_by_gene(rho)
    gene_mean = {}
    gene_n = {}
    for gene, vals in grouped.items():
        gene_n[gene] = len(vals)
        gene_mean[gene] = float(np.mean(vals)) if vals else np.nan
    stored = assignment.loc[np.isfinite(assignment["l2"]), ["gene_symbol", "l2"]]
    stored_l2 = {}
    check_vals = []
    for gene, sub in stored.groupby("gene_symbol", sort=False):
        vals = sub["l2"].to_numpy(dtype=np.float64)
        if float(vals.max() - vals.min()) != 0.0:
            raise SystemExit(f"stored l2 is not unique for {gene}")
        stored_l2[gene] = float(vals[0])
        got = gene_mean.get(gene, np.nan)
        if not np.isfinite(got):
            check_vals.append(np.inf)
        else:
            check_vals.append(abs(stored_l2[gene] - got))
    l2_gap = float(np.max(check_vals)) if check_vals else np.nan
    # The stored l2 is the Lpool value. A fresh mean of the same dictionary can
    # differ in the last bit. The column uses the stored value when it exists.
    mapping = PLAN["tcga_to_cptac"]
    cohort = source["cancer_type"].astype(str).map(mapping)
    pooled = source["gene_symbol"].astype(str).map(lambda gene: stored_l2.get(gene, gene_mean.get(gene, np.nan)))
    n_cohorts = source["gene_symbol"].astype(str).map(gene_n).fillna(0).astype(int)
    matched = []
    for gene, code in zip(source["gene_symbol"].astype(str), cohort):
        if not isinstance(code, str) or code == "" or (isinstance(code, float)):
            matched.append(np.nan)
            continue
        matched.append(rho.get((code, gene), np.nan))
    out = source.copy()
    out["cptac_pooled_rho"] = pooled.to_numpy(dtype=np.float64)
    out["cptac_n_cohorts"] = n_cohorts.to_numpy(dtype=np.int64)
    out["cptac_matched_cohort"] = cohort.fillna("").to_numpy()
    out["cptac_matched_rho"] = np.asarray(matched, dtype=np.float64)
    pred_gap = 0.0
    for column in PLAN["prediction_columns"]:
        gap = np.abs(out[column].to_numpy(dtype=np.float64) - source[column].to_numpy(dtype=np.float64))
        pred_gap = max(pred_gap, float(np.nanmax(gap)))
    if pred_gap != 0.0:
        raise SystemExit(f"prediction columns differ: max abs {pred_gap}")
    stored_series = out["gene_symbol"].astype(str).map(stored_l2)
    both = stored_series.notna()
    column_gap = float(np.max(np.abs(out.loc[both, "cptac_pooled_rho"].to_numpy(dtype=np.float64) - stored_series[both].to_numpy(dtype=np.float64))))
    if column_gap != 0.0:
        raise SystemExit(f"stored l2 column check failed: max abs {column_gap}")
    for column in source.columns:
        if column in PLAN["prediction_columns"]:
            continue
        if not np.array_equal(out[column].to_numpy(), source[column].to_numpy()):
            raise SystemExit(f"v1 column changed: {column}")
    path = ROOT / PLAN["output"]
    out.to_csv(path, sep="\t", index=False, compression="gzip", na_rep="")
    genes = out["gene_symbol"].astype(str)
    pooled_missing_genes = int(genes[out["cptac_pooled_rho"].isna()].nunique())
    matched_rho_missing_genes = int(genes[out["cptac_matched_rho"].isna() & out["cptac_matched_cohort"].ne("")].nunique())
    unmapped_cancers = sorted(names - set(mapping))
    check = pd.DataFrame(
        [
            {
                "n_rows": int(len(out)),
                "n_genes": int(genes.nunique()),
                "n_cancers": int(out["cancer_type"].nunique()),
                "size_bytes": int(path.stat().st_size),
                "prediction_max_abs": pred_gap,
                "recomputed_mean_vs_stored_l2_max_abs": l2_gap,
                "column_vs_stored_l2_max_abs": column_gap,
                "n_l2_genes_checked": int(stored["gene_symbol"].nunique()),
                "extra_cptac_rho_pairs": int(n_added),
                "n_p_gene_rho_values": int(n_g),
                "n_genes_pooled_missing": pooled_missing_genes,
                "n_genes_n_cohorts_zero": int(genes[out["cptac_n_cohorts"] == 0].nunique()),
                "n_rows_matched_cohort_empty": int((out["cptac_matched_cohort"] == "").sum()),
                "n_genes_with_any_matched_cohort": int(genes[out["cptac_matched_cohort"] != ""].nunique()),
                "n_genes_matched_rho_missing_on_a_mapped_cancer": matched_rho_missing_genes,
                "n_rows_matched_rho_missing": int(out["cptac_matched_rho"].isna().sum()),
                "n_unmapped_cancers": len(unmapped_cancers),
            }
        ]
    )
    check.to_csv(OUT / "forecast_v2_check.tsv", sep="\t", index=False)
    print(check.to_string(index=False), flush=True)
    print("stage11 forecast v2 done", flush=True)


if __name__ == "__main__":
    main()
