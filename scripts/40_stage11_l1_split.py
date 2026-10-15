"""Post-hoc split of L1 pairs by whether F0 is the global mean.

The plan file config/analysis_plan_B11.yaml was committed before this script runs.
Stored predictions and stored lookup columns are not recomputed. The model is not fit.
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

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B11.yaml"))
SPEC = PLAN["l1_split"]
FEAT = SPEC["feature_default_check"]
OUT = ROOT / "results" / "stage11"
LABEL = PLAN["label"]
SEED = int(PLAN["seed"])
DRAWS = int(PLAN["bootstrap_draws"])
TOL = float(FEAT["tolerance"])


def write(name, frame):
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / name, sep="\t", index=False)


def contrast(sub, layer, left_name, right_name, left, right):
    point = lookup.finite_spearman(left, sub["actual"], 3) - lookup.finite_spearman(right, sub["actual"], 3)
    groups = lookup.groups_of(sub) if len(sub) else []
    draws = np.full(DRAWS, np.nan)
    if np.isfinite(point) and groups:
        rng = np.random.default_rng(SEED)
        actual = sub["actual"].to_numpy(dtype=np.float64)
        for i in range(DRAWS):
            idx = lookup.draw_indices(groups, rng)
            draws[i] = lookup.finite_spearman(left[idx], actual[idx], 3) - lookup.finite_spearman(
                right[idx], actual[idx], 3
            )
    ok = draws[np.isfinite(draws)]
    return {
        "label": LABEL,
        "dataset": sub["dataset"].iloc[0] if len(sub) else "",
        "layer": layer,
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


def is_default(frame, gmean):
    if frame.empty:
        return np.array([], dtype=bool)
    n_def = frame["n_defined_cancers"].to_numpy(dtype=np.float64)
    rho_sd = frame["rho_sd"].to_numpy(dtype=np.float64)
    rho_mean = frame["rho_mean"].to_numpy(dtype=np.float64)
    rho_min = frame["rho_min"].to_numpy(dtype=np.float64)
    rho_max = frame["rho_max"].to_numpy(dtype=np.float64)
    return (
        (n_def == 0)
        & (rho_sd == 0)
        & (np.abs(rho_mean - gmean) <= TOL)
        & (np.abs(rho_min - gmean) <= TOL)
        & (np.abs(rho_max - gmean) <= TOL)
    )


def join_cached(block, cache, left_keys, right_keys):
    first = cache[right_keys + FEAT["columns"]].drop_duplicates(right_keys, keep="first")
    counts = cache.groupby(right_keys, sort=False).size().rename("n_feature_copies").reset_index()
    merged = block.merge(first, left_on=left_keys, right_on=right_keys, how="left")
    merged = merged.merge(counts, on=right_keys, how="left")
    return merged


def apollo_gene_level(genes, gmean_stored):
    stage3 = importlib.import_module("scripts.08_stage3")
    _genes, p_genes, adc = stage3.load_lists()
    rho = stage3.defined_rho(p_genes)
    pset = set(p_genes)
    adc_set = set(adc)
    loco = pd.read_csv(ROOT / "results" / "adc_loco.tsv", sep="\t")
    for rec in loco.itertuples(index=False):
        if rec.gene_symbol in pset or rec.gene_symbol not in adc_set:
            continue
        if int(rec.n) >= stage3.MIN_N and np.isfinite(rec.s0):
            rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
    rho = stage3.forecast_rho(rho, genes)
    gmean = stage3.gmean_all(rho, p_genes)
    rows = []
    for gene in genes:
        vals = [rho[(c, gene)] for c in stage3.CANCERS if (c, gene) in rho]
        if vals:
            rows.append(
                {
                    "gene_symbol": gene,
                    "rho_mean": float(np.mean(vals)),
                    "rho_sd": float(np.std(vals, ddof=stage3.DDOF)) if len(vals) >= 2 else 0.0,
                    "rho_min": float(np.min(vals)),
                    "rho_max": float(np.max(vals)),
                    "n_defined_cancers": len(vals),
                }
            )
        else:
            rows.append(
                {
                    "gene_symbol": gene,
                    "rho_mean": gmean,
                    "rho_sd": 0.0,
                    "rho_min": gmean,
                    "rho_max": gmean,
                    "n_defined_cancers": 0,
                }
            )
    frame = pd.DataFrame(rows)
    frame["gmean_recovered"] = gmean
    frame["gmean_stored"] = gmean_stored
    return frame


def feature_status(dataset, block, gmean, caches):
    if dataset == "RPPA":
        merged = join_cached(block, caches["RPPA"], ["unit", "context"], ["xena_antibody", "cancer"])
    elif dataset == "ProCan":
        merged = join_cached(block, caches["ProCan"], ["unit", "context"], ["protein_id", "lineage"])
    else:
        gene_rows = caches["APOLLO"]
        merged = block.merge(gene_rows, on="gene_symbol", how="left")
        merged["n_feature_copies"] = np.where(merged["n_defined_cancers"].notna(), 1, np.nan)
    present = merged["n_defined_cancers"].notna().to_numpy()
    copies = merged["n_feature_copies"].fillna(0).to_numpy(dtype=np.float64)
    default = np.zeros(len(merged), dtype=bool)
    if present.any():
        default_part = is_default(merged.loc[present], gmean)
        default[present] = default_part.to_numpy() if hasattr(default_part, "to_numpy") else default_part
    mismatch = (~present) | (copies != 1) | (present & ~default)
    return int((~present).sum()), int(mismatch.sum())


def main():
    pairs = pd.read_csv(ROOT / SPEC["assignment"], sep="\t")
    build = pd.read_csv(ROOT / SPEC["build"], sep="\t")
    gmean = float(build[SPEC["global_mean_column"]].iloc[0])
    p_genes = set(g for g in (ROOT / SPEC["p_genes"]).read_text().split() if g)
    caches = {
        "RPPA": pd.read_csv(ROOT / FEAT["rppa"], sep="\t"),
        "ProCan": pd.read_csv(ROOT / FEAT["procan"], sep="\t"),
    }
    apollo_genes = pd.read_csv(ROOT / "results" / "stage6" / "prep" / "U_P2_genes.tsv", sep="\t")["gene_symbol"].tolist()
    print("recovering APOLLO gene-level columns", flush=True)
    caches["APOLLO"] = apollo_gene_level(apollo_genes, gmean)
    rows = []
    counts = []
    finite_cols = SPEC["require_finite"]
    for dataset, block in pairs.groupby("dataset", sort=True):
        sub = block[block["layer"] == "L1"].copy()
        keep = np.ones(len(sub), dtype=bool)
        for column in finite_cols:
            keep &= np.isfinite(sub[column].to_numpy(dtype=np.float64))
        sub = sub.loc[keep].reset_index(drop=True)
        is_b = np.abs(sub["f0"].to_numpy(dtype=np.float64) - gmean) <= TOL
        for layer, mask in (("a", ~is_b), ("b", is_b)):
            part = sub.loc[mask].reset_index(drop=True)
            in_p = part["gene_symbol"].isin(p_genes)
            missing = mismatch = 0
            if layer == "b":
                missing, mismatch = feature_status(dataset, part, gmean, caches)
            counts.append(
                {
                    "label": LABEL,
                    "dataset": dataset,
                    "layer": layer,
                    "n_pairs": int(len(part)),
                    "n_units": int(part["unit"].nunique()) if len(part) else 0,
                    "n_genes": int(part["gene_symbol"].nunique()) if len(part) else 0,
                    "n_genes_in_P": int(part.loc[in_p, "gene_symbol"].nunique()) if len(part) else 0,
                    "n_feature_missing": missing if layer == "b" else "",
                    "n_feature_mismatch": mismatch if layer == "b" else "",
                    "global_mean": gmean,
                }
            )
            score = {
                "f1s_lgbm": part["f1s_lgbm"].to_numpy(dtype=np.float64) if len(part) else np.array([]),
                "lpool": part["l2"].to_numpy(dtype=np.float64) if len(part) else np.array([]),
                "l1": part["l1"].to_numpy(dtype=np.float64) if len(part) else np.array([]),
            }
            for left_name, right_name in SPEC["contrasts"]:
                if len(part) == 0:
                    rows.append(
                        {
                            "label": LABEL,
                            "dataset": dataset,
                            "layer": layer,
                            "left": left_name,
                            "right": right_name,
                            "diff": np.nan,
                            "defined": 0,
                            "ci_low": np.nan,
                            "ci_high": np.nan,
                            "n_nonfinite_draws": DRAWS,
                            "n_pairs": 0,
                            "n_units": 0,
                            "n_genes": 0,
                        }
                    )
                    continue
                rows.append(contrast(part, layer, left_name, right_name, score[left_name], score[right_name]))
            print(dataset, layer, len(part), flush=True)
    write("l1_split_diff.tsv", pd.DataFrame(rows))
    write("l1_split_counts.tsv", pd.DataFrame(counts))
    g_gap = abs(float(caches["APOLLO"]["gmean_recovered"].iloc[0]) - gmean)
    print("apollo gmean abs diff", g_gap, flush=True)
    print("stage11 l1 split done", flush=True)


if __name__ == "__main__":
    main()
