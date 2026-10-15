"""Post-hoc F0 minus L* on the stored stage-8 lookup table.

The plan file is config/analysis_plan_B9.yaml. It was committed before this script runs.
This script does not recompute L* and does not refit models.
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

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B9.yaml"))
OUT = ROOT / "results" / "stage9"
LABEL = "사후"
SEED = int(PLAN["seed"])
DRAWS = int(PLAN["bootstrap_draws"])
DIFFS = [tuple(pair) for pair in PLAN["contrasts"]["differences"]]
DIFF_STRATA = list(PLAN["contrasts"]["difference_strata"])
COUNT_STRATA = list(PLAN["contrasts"]["count_only_strata"])


def write(name, frame):
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / name, sep="\t", index=False)


def subset(frame, stratum):
    if stratum == "all":
        return frame
    return frame[frame["layer"] == stratum]


def main():
    pairs = pd.read_csv(ROOT / "results" / "stage8" / "lookup_assignment.tsv", sep="\t")
    rows = []
    counts = []
    for dataset, block in pairs.groupby("dataset", sort=True):
        for stratum in DIFF_STRATA + COUNT_STRATA:
            base = subset(block, stratum)
            finite_actual = base[np.isfinite(base["actual"])].copy()
            sub = finite_actual[np.isfinite(finite_actual["f0"]) & np.isfinite(finite_actual["Lstar"]) & np.isfinite(finite_actual["f1s_lgbm"])].reset_index(drop=True)
            counts.append(
                {
                    "label": LABEL,
                    "dataset": dataset,
                    "stratum": stratum,
                    "n_pairs": int(len(sub)),
                    "n_units": int(sub["unit"].nunique()) if len(sub) else 0,
                    "n_genes": int(sub["gene_symbol"].nunique()) if len(sub) else 0,
                }
            )
            if stratum in COUNT_STRATA:
                continue
            for left, right in DIFFS:
                point = lookup.finite_spearman(sub[left], sub["actual"], 3) - lookup.finite_spearman(sub[right], sub["actual"], 3) if len(sub) else np.nan
                groups = lookup.groups_of(sub) if len(sub) else []
                draws = np.full(DRAWS, np.nan)
                if np.isfinite(point) and groups:
                    rng = np.random.default_rng(SEED)
                    actual = sub["actual"].to_numpy(dtype=np.float64)
                    a = sub[left].to_numpy(dtype=np.float64)
                    b = sub[right].to_numpy(dtype=np.float64)
                    for i in range(DRAWS):
                        idx = lookup.draw_indices(groups, rng)
                        draws[i] = lookup.finite_spearman(a[idx], actual[idx], 3) - lookup.finite_spearman(b[idx], actual[idx], 3)
                ok = draws[np.isfinite(draws)]
                rows.append(
                    {
                        "label": LABEL,
                        "dataset": dataset,
                        "stratum": stratum,
                        "left": left,
                        "right": right,
                        "diff": point,
                        "defined": int(np.isfinite(point)),
                        "ci_low": float(np.quantile(ok, 0.025)) if ok.size else np.nan,
                        "ci_high": float(np.quantile(ok, 0.975)) if ok.size else np.nan,
                        "n_nonfinite_draws": int(DRAWS - ok.size) if np.isfinite(point) else DRAWS,
                        "n_pairs": int(len(sub)),
                        "n_units": int(sub["unit"].nunique()) if len(sub) else 0,
                        "n_genes": int(sub["gene_symbol"].nunique()) if len(sub) else 0,
                    }
                )
                print(dataset, stratum, left, right, point, flush=True)
    write("f0_lstar_diff.tsv", pd.DataFrame(rows))
    write("f0_lstar_counts.tsv", pd.DataFrame(counts))
    print("stage9 f0-lstar done", flush=True)


if __name__ == "__main__":
    main()
