"""All-gene TCGA forecasts from the exported F1s model and the ADC-atlas feature code.

The plan file config/analysis_plan_B10.yaml was committed before this script runs.
The script does not refit and does not rewrite the atlas or frozen predictions.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import importlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml
from lightgbm import Booster

from scripts.common import ROOT

atlas = importlib.import_module("scripts.13_stage5_atlas")
scope = importlib.import_module("scripts.12_stage5_scope")

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B10.yaml"))["forecast_all_genes"]
OUT = ROOT / "results" / "stage10"
FEATURES = json.loads((ROOT / "release" / "models" / "F1s_lightgbm_features.json").read_text())


def predict_filled(booster, frame, med):
    x = frame[FEATURES["features"]].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(np.asarray(med, dtype=np.float64), inds[1])
    return booster.predict(filled)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fit_table = pd.read_csv(ROOT / "results" / "f1s_model_fit.tsv", sep="\t")
    stored = dict(zip(fit_table["item"], fit_table["value"]))
    half80 = float(stored["q80_lgbm"])
    half90 = float(stored["q90_lgbm"])
    p_genes = set(g for g in (ROOT / "results" / "P_genes.txt").read_text().split() if g)
    stage3 = atlas.load_stage3()
    genes = [g for g in (ROOT / "results" / "G_RNA_genes.txt").read_text().split() if g]
    print("static features", len(genes), flush=True)
    features = atlas.static_features(stage3, genes)
    booster = Booster(model_file=str(ROOT / PLAN["model"]))
    med = FEATURES["median_fill"]
    pct = pd.read_csv(ROOT / PLAN["percentile_file"], sep="\t")
    print("percentiles", pct["disease"].nunique(), pct["gene_symbol"].nunique(), len(pct), flush=True)
    rows = []
    for disease, sub in pct.groupby("disease", sort=True):
        records = []
        keep = []
        for rec in sub.itertuples(index=False):
            base = features.get(rec.gene_symbol)
            if base is None:
                continue
            row = dict(base)
            row["holdout_rna_mean_percentile"] = rec.median_percentile
            row["holdout_rna_sd_percentile"] = rec.sd_percentile
            row["gene_symbol"] = rec.gene_symbol
            records.append(row)
            keep.append(rec.gene_symbol)
        frame = pd.DataFrame(records)
        hat = predict_filled(booster, frame, med)
        for gene, value, rec in zip(keep, hat, records):
            value = float(value)
            rows.append(
                {
                    "gene_symbol": gene,
                    "cancer_type": disease,
                    "forecast": value,
                    "lower80": value - half80,
                    "upper80": value + half80,
                    "lower90": value - half90,
                    "upper90": value + half90,
                    "in_P": int(gene in p_genes),
                    "n_defined_cancers": rec["n_defined_cancers"],
                }
            )
        print("forecast", disease, len(keep), flush=True)
    out = pd.DataFrame(rows)
    print("similarity", flush=True)
    c09 = scope.load_09()
    cptac = scope.cptac_profiles(genes)
    s_rows = []
    low = float(PLAN["s_range_low"])
    high = float(PLAN["s_range_high"])
    for disease, sub in pct.groupby("disease", sort=True):
        profile = {r.gene_symbol: r.median_percentile for r in sub.itertuples(index=False)}
        s_value = scope.max_similarity(c09, profile, cptac, genes)
        inside = int(np.isfinite(s_value) and low <= s_value <= high)
        s_rows.append({"cancer_type": disease, "s_context": s_value, "s_in_cptac_range": inside})
        print("s", disease, s_value, inside, flush=True)
    s_frame = pd.DataFrame(s_rows)
    s_frame.to_csv(OUT / "s_by_cancer.tsv", sep="\t", index=False)
    out = out.merge(s_frame, on="cancer_type", how="left")
    path = ROOT / PLAN["output"]
    out.to_csv(path, sep="\t", index=False, compression="gzip")
    atlas_path = ROOT / PLAN["match_file"]
    atlas_df = pd.read_csv(atlas_path, sep="\t")
    merged = out.merge(atlas_df, left_on=["gene_symbol", "cancer_type"], right_on=["gene_symbol", "disease"], how="inner")
    gap = np.abs(merged["forecast"].to_numpy(dtype=np.float64) - merged[PLAN["match_column"]].to_numpy(dtype=np.float64))
    finite = gap[np.isfinite(gap)]
    match = pd.DataFrame(
        [
            {
                "n_atlas": int(len(atlas_df)),
                "n_matched": int(len(merged)),
                "max_abs": float(finite.max()) if finite.size else np.nan,
                "limit": float(PLAN["match_limit"]),
                "passed": int(finite.size and float(finite.max()) <= float(PLAN["match_limit"]) and len(merged) == len(atlas_df)),
                "n_genes": int(out["gene_symbol"].nunique()),
                "n_cancers": int(out["cancer_type"].nunique()),
                "n_rows": int(len(out)),
                "n_s": int(s_frame["s_context"].notna().sum()),
                "n_s_in_range": int(s_frame["s_in_cptac_range"].sum()),
                "size_bytes": int(path.stat().st_size),
            }
        ]
    )
    match.to_csv(OUT / "forecast_match.tsv", sep="\t", index=False)
    print(match.to_string(index=False), flush=True)
    print("stage10 forecast done", flush=True)


if __name__ == "__main__":
    main()
