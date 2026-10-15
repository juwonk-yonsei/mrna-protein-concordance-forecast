"""Post-hoc pooled lookup versus matched lookup on stored stage-8 rows.

The plan file is config/analysis_plan_B10.yaml. It was committed before this script runs.
Lpool is the stored l2 column. This script does not recompute rho and does not refit.
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

from scripts.common import ROOT

lookup = importlib.import_module("scripts.25_stage8_lookup")

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B10.yaml"))
OUT = ROOT / "results" / "stage10"
LABEL = "사후"
SEED = int(PLAN["seed"])
DRAWS = int(PLAN["bootstrap_draws"])
SPEC = PLAN["lpool"]
COVER = PLAN["f0_coverage"]


def write(name, frame):
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / name, sep="\t", index=False)


def contrast(sub, left_name, right_name, left, right):
    point = lookup.finite_spearman(left, sub["actual"], 3) - lookup.finite_spearman(right, sub["actual"], 3)
    groups = lookup.groups_of(sub)
    draws = np.full(DRAWS, np.nan)
    if np.isfinite(point) and groups:
        rng = np.random.default_rng(SEED)
        actual = sub["actual"].to_numpy(dtype=np.float64)
        for i in range(DRAWS):
            idx = lookup.draw_indices(groups, rng)
            draws[i] = lookup.finite_spearman(left[idx], actual[idx], 3) - lookup.finite_spearman(right[idx], actual[idx], 3)
    ok = draws[np.isfinite(draws)]
    return {
        "label": LABEL,
        "dataset": sub["dataset"].iloc[0],
        "stratum": "L1",
        "left": left_name,
        "right": right_name,
        "diff": point,
        "defined": int(np.isfinite(point)),
        "ci_low": float(np.quantile(ok, 0.025)) if ok.size else np.nan,
        "ci_high": float(np.quantile(ok, 0.975)) if ok.size else np.nan,
        "n_nonfinite_draws": int(DRAWS - ok.size) if np.isfinite(point) else DRAWS,
        "n_pairs": int(len(sub)),
        "n_units": int(sub["unit"].nunique()) if len(sub) else 0,
        "n_genes": int(sub["gene_symbol"].nunique()) if len(sub) else 0,
    }


def main():
    pairs = pd.read_csv(ROOT / SPEC["assignment"], sep="\t")
    build = pd.read_csv(ROOT / SPEC["build"], sep="\t")
    gmean = float(build[COVER["global_mean_column"]].iloc[0])
    rows = []
    counts = []
    for dataset, block in pairs.groupby("dataset", sort=True):
        finite_actual = block[np.isfinite(block["actual"])].copy()
        f0_global = np.abs(finite_actual["f0"].to_numpy(dtype=np.float64) - gmean) <= 1e-8
        has_cptac = finite_actual["layer"].isin(["L1", "L2", "L3"]).to_numpy()
        for stratum in COVER["strata"]:
            if stratum == "all":
                mask = np.ones(len(finite_actual), dtype=bool)
            else:
                mask = (finite_actual["layer"] == stratum).to_numpy()
            hit = finite_actual[mask & f0_global & has_cptac]
            counts.append(
                {
                    "label": LABEL,
                    "dataset": dataset,
                    "stratum": stratum,
                    "n_pairs": int(len(hit)),
                    "n_units": int(hit["unit"].nunique()) if len(hit) else 0,
                    "n_genes": int(hit["gene_symbol"].nunique()) if len(hit) else 0,
                    "global_mean": gmean,
                }
            )
        sub = block[block["layer"] == "L1"].copy()
        sub = sub[np.isfinite(sub["actual"]) & np.isfinite(sub["l1"]) & np.isfinite(sub["l2"]) & np.isfinite(sub["f1s_lgbm"])].reset_index(drop=True)
        left = sub["l2"].to_numpy(dtype=np.float64)
        matched = sub["l1"].to_numpy(dtype=np.float64)
        forecast = sub["f1s_lgbm"].to_numpy(dtype=np.float64)
        rows.append(contrast(sub, "lpool", "l1", left, matched))
        rows.append(contrast(sub, "f1s_lgbm", "lpool", forecast, left))
        print(dataset, "L1", rows[-2]["diff"], rows[-1]["diff"], flush=True)
    write("lpool_diff.tsv", pd.DataFrame(rows))
    write("f0_coverage_counts.tsv", pd.DataFrame(counts))
    print("stage10 lpool done", flush=True)


if __name__ == "__main__":
    main()
