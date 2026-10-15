"""Post-hoc matched-context lookup L* and the HYB rule.

The cascade and the context map are read from config/analysis_plan_B8.yaml.
That file was committed before this script is run. This script does not refit models.
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
from scipy.stats import binomtest

from scripts.b2_lib import CANCERS
from scripts.common import ROOT

_freeze = importlib.import_module("scripts.07_freeze")
extra_cptac_rho = _freeze.extra_cptac_rho
load_defined_rho = _freeze.load_defined_rho
load_lists = _freeze.load_lists

LABEL = "사후"
OUT = ROOT / "results" / "stage8"
PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B8.yaml"))
SEED = int(PLAN["seed"])
DRAWS = int(PLAN["bootstrap_draws"])
CONTRASTS = PLAN["contrasts"]


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
    if score.size == 0 or score.max() == score.min():
        return np.nan
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = average_ranks(score)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def write(name, frame):
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / name, sep="\t", index=False)


def lookup_rho(eval_genes):
    _genes, p_genes, _adc = load_lists()
    rho = load_defined_rho(p_genes, _adc)
    rho, n_added = extra_cptac_rho(sorted(eval_genes), rho)
    print(f"extra cptac rho pairs {n_added}", flush=True)
    pset = set(p_genes)
    by_gene = {g: [] for g in eval_genes}
    gvals = []
    for (cancer, gene), value in rho.items():
        if cancer not in CANCERS:
            continue
        if gene in by_gene:
            by_gene[gene].append(float(value))
        if gene in pset:
            gvals.append(float(value))
    l2 = {g: float(np.mean(vs)) for g, vs in by_gene.items() if vs}
    gmean = float(np.mean(gvals))
    return rho, l2, gmean, n_added, len(gvals)


def l1_value(rho, gene, cohorts):
    vals = [rho[(c, gene)] for c in cohorts if (c, gene) in rho]
    if not vals:
        return np.nan
    return float(np.mean(vals))


def cohorts_for(dataset, context):
    if dataset == "RPPA":
        name = PLAN["rppa_context_to_cptac"].get(context)
        return [name] if name else []
    if dataset == "ProCan":
        return list(PLAN["procan_context_to_cptac"].get(context, []))
    return list(PLAN["apollo_context_to_cptac"].get(context, []))


def load_frames():
    rppa = pd.read_csv(
        ROOT / "results" / "stage4" / "confirm" / "rppa_pairs.tsv",
        sep="\t",
        usecols=["xena_antibody", "gene_symbol", "cancer", "actual", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"],
    )
    rppa["dataset"] = "RPPA"
    rppa["unit"] = rppa["xena_antibody"].astype(str)
    rppa["context"] = rppa["cancer"].astype(str)
    pro = pd.read_csv(
        ROOT / "results" / "stage4" / "confirm" / "procan_pairs.tsv",
        sep="\t",
        usecols=["protein_id", "gene_symbol", "lineage", "actual", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"],
    )
    pro["dataset"] = "ProCan"
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
    apollo["dataset"] = "APOLLO"
    apollo["unit"] = apollo["gene_symbol"].astype(str)
    apollo["context"] = "LUAD"
    return pd.concat(
        [
            rppa[["dataset", "unit", "gene_symbol", "context", "actual", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"]],
            pro[["dataset", "unit", "gene_symbol", "context", "actual", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"]],
            apollo[["dataset", "unit", "gene_symbol", "context", "actual", "f0", "f1s_ridge", "f1s_lgbm", "b_sd"]],
        ],
        ignore_index=True,
    )


def assign_lookup(frame, rho, l2, gmean):
    pro = frame[frame["dataset"] == "ProCan"]
    l3 = {}
    for gene, sub in pro.groupby("gene_symbol", sort=False):
        vals = sub["actual"].to_numpy(dtype=np.float64)
        vals = vals[np.isfinite(vals)]
        if vals.size:
            l3[gene] = float(np.mean(vals))
    rows_l1 = []
    rows_l2 = []
    rows_l3 = []
    layers = []
    lstars = []
    hyb = []
    for rec in frame.itertuples(index=False):
        cohorts = cohorts_for(rec.dataset, rec.context)
        v1 = l1_value(rho, rec.gene_symbol, cohorts)
        v2 = l2.get(rec.gene_symbol, np.nan)
        v3 = l3.get(rec.gene_symbol, np.nan) if rec.dataset != "ProCan" else np.nan
        if np.isfinite(v1):
            layer, star = "L1", v1
        elif np.isfinite(v2):
            layer, star = "L2", v2
        elif rec.dataset != "ProCan" and np.isfinite(v3):
            layer, star = "L3", v3
        else:
            layer, star = "L4", gmean
        rows_l1.append(v1)
        rows_l2.append(v2)
        rows_l3.append(v3)
        layers.append(layer)
        lstars.append(star)
        hyb.append(v1 if np.isfinite(v1) else rec.f1s_lgbm)
    out = frame.copy()
    out["l1"] = rows_l1
    out["l2"] = rows_l2
    out["l3"] = rows_l3
    out["l4"] = gmean
    out["layer"] = layers
    out["Lstar"] = lstars
    out["HYB"] = hyb
    out.insert(0, "label", LABEL)
    return out


def groups_of(frame):
    return [np.asarray(v, dtype=np.int64) for v in frame.groupby("unit", sort=True).indices.values()]


def draw_indices(groups, rng):
    chosen = rng.integers(0, len(groups), size=len(groups))
    return np.concatenate([groups[i] for i in chosen])


def subset(frame, stratum):
    if stratum == "all":
        return frame
    return frame[frame["layer"] == stratum]


def metric_rows(pairs):
    point_rows = []
    diff_rows = []
    scores = ("f1s_lgbm", "f1s_ridge", "f0", "Lstar", "HYB", "b_sd")
    for dataset, block in pairs.groupby("dataset", sort=True):
        for stratum in CONTRASTS["strata"]:
            base = subset(block, stratum)
            finite_actual = base[np.isfinite(base["actual"])].copy()
            for score in scores:
                sub = finite_actual[np.isfinite(finite_actual[score])].reset_index(drop=True)
                rho = finite_spearman(sub[score], sub["actual"], 3) if len(sub) else np.nan
                point_rows.append(
                    {
                        "label": LABEL,
                        "dataset": dataset,
                        "stratum": stratum,
                        "score": score,
                        "spearman": rho,
                        "defined": int(np.isfinite(rho)),
                        "n_pairs": int(len(sub)),
                        "n_genes": int(sub["gene_symbol"].nunique()) if len(sub) else 0,
                        "n_units": int(sub["unit"].nunique()) if len(sub) else 0,
                    }
                )
            for left, right in CONTRASTS["spearman_differences"]:
                sub = finite_actual[np.isfinite(finite_actual[left]) & np.isfinite(finite_actual[right])].reset_index(drop=True)
                point = finite_spearman(sub[left], sub["actual"], 3) - finite_spearman(sub[right], sub["actual"], 3) if len(sub) else np.nan
                groups = groups_of(sub) if len(sub) else []
                draws = np.full(DRAWS, np.nan)
                if np.isfinite(point) and groups:
                    rng = np.random.default_rng(SEED)
                    actual = sub["actual"].to_numpy(dtype=np.float64)
                    a = sub[left].to_numpy(dtype=np.float64)
                    b = sub[right].to_numpy(dtype=np.float64)
                    for i in range(DRAWS):
                        idx = draw_indices(groups, rng)
                        draws[i] = finite_spearman(a[idx], actual[idx], 3) - finite_spearman(b[idx], actual[idx], 3)
                ok = draws[np.isfinite(draws)]
                diff_rows.append(
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
                        "n_genes": int(sub["gene_symbol"].nunique()) if len(sub) else 0,
                        "n_units": int(sub["unit"].nunique()) if len(sub) else 0,
                    }
                )
                print(dataset, stratum, left, right, point, flush=True)
    write("spearman.tsv", pd.DataFrame(point_rows))
    write("spearman_diff.tsv", pd.DataFrame(diff_rows))


def auroc_rows(pairs):
    point_rows = []
    diff_rows = []
    for dataset, block in pairs.groupby("dataset", sort=True):
        for stratum in CONTRASTS["strata"]:
            base = subset(block, stratum)
            finite_actual = base[np.isfinite(base["actual"])].copy()
            for score in CONTRASTS["auroc_scores"]:
                sub = finite_actual[np.isfinite(finite_actual[score])].reset_index(drop=True)
                auc = auroc_high(sub[score], sub["actual"] >= 0.5) if len(sub) else np.nan
                point_rows.append(
                    {
                        "label": LABEL,
                        "dataset": dataset,
                        "stratum": stratum,
                        "score": score,
                        "auroc": auc,
                        "defined": int(np.isfinite(auc)),
                        "n_pairs": int(len(sub)),
                        "n_positive": int((sub["actual"] >= 0.5).sum()) if len(sub) else 0,
                        "n_genes": int(sub["gene_symbol"].nunique()) if len(sub) else 0,
                        "n_units": int(sub["unit"].nunique()) if len(sub) else 0,
                    }
                )
            for left, right in CONTRASTS["auroc_differences"]:
                sub = finite_actual[np.isfinite(finite_actual[left]) & np.isfinite(finite_actual[right])].reset_index(drop=True)
                if len(sub) == 0:
                    point = np.nan
                else:
                    point = auroc_high(sub[left], sub["actual"] >= 0.5) - auroc_high(sub[right], sub["actual"] >= 0.5)
                groups = groups_of(sub) if len(sub) else []
                draws = np.full(DRAWS, np.nan)
                if np.isfinite(point) and groups:
                    rng = np.random.default_rng(SEED)
                    actual = sub["actual"].to_numpy(dtype=np.float64)
                    a = sub[left].to_numpy(dtype=np.float64)
                    b = sub[right].to_numpy(dtype=np.float64)
                    for i in range(DRAWS):
                        idx = draw_indices(groups, rng)
                        draws[i] = auroc_high(a[idx], actual[idx] >= 0.5) - auroc_high(b[idx], actual[idx] >= 0.5)
                ok = draws[np.isfinite(draws)]
                diff_rows.append(
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
                    }
                )
                print("auroc", dataset, stratum, left, right, point, flush=True)
    write("auroc.tsv", pd.DataFrame(point_rows))
    write("auroc_diff.tsv", pd.DataFrame(diff_rows))


def within_context(pairs):
    rows = []
    for dataset, block in pairs.groupby("dataset", sort=True):
        finite_actual = block[np.isfinite(block["actual"])]
        for context, sub in finite_actual.groupby("context", sort=True):
            rec = {"label": LABEL, "dataset": dataset, "context": context, "n_pairs": int(len(sub))}
            for score in ("f1s_lgbm", "f1s_ridge", "f0", "Lstar", "b_sd"):
                rec[score] = finite_spearman(sub[score], sub["actual"], 3)
            for left, right in CONTRASTS["within_context_differences"]:
                rec[f"{left}_minus_{right}"] = rec[left] - rec[right] if np.isfinite(rec[left]) and np.isfinite(rec[right]) else np.nan
            rows.append(rec)
    detail = pd.DataFrame(rows)
    write("within_context.tsv", detail)
    summary = []
    rppa = detail[detail["dataset"] == "RPPA"]
    for left, right in CONTRASTS["within_context_differences"]:
        col = f"{left}_minus_{right}"
        vals = rppa[col].to_numpy(dtype=np.float64)
        vals = vals[np.isfinite(vals)]
        nonzero = vals[vals != 0]
        n_pos = int((nonzero > 0).sum())
        n_tie = int(vals.size - nonzero.size)
        if nonzero.size:
            p = float(binomtest(n_pos, int(nonzero.size), 0.5, alternative="two-sided").pvalue)
        else:
            p = np.nan
        rng = np.random.default_rng(SEED)
        draws = np.empty(DRAWS)
        if vals.size:
            for i in range(DRAWS):
                draws[i] = float(np.median(vals[rng.integers(0, vals.size, size=vals.size)]))
            ok = draws[np.isfinite(draws)]
            lo = float(np.quantile(ok, 0.025))
            hi = float(np.quantile(ok, 0.975))
        else:
            lo = hi = np.nan
        summary.append(
            {
                "label": LABEL,
                "dataset": "RPPA",
                "contrast": col,
                "n_contexts_defined": int(vals.size),
                "n_positive": n_pos,
                "n_negative": int((nonzero < 0).sum()),
                "n_zero": n_tie,
                "sign_test_p": p,
                "median_diff": float(np.median(vals)) if vals.size else np.nan,
                "ci_low": lo,
                "ci_high": hi,
                "bootstrap_draws": DRAWS,
                "bootstrap_seed": SEED,
            }
        )
    write("within_context_sign.tsv", pd.DataFrame(summary))


def both_f0_counts():
    fill = pd.read_csv(ROOT / "results" / "stage7" / "f0_fill.tsv", sep="\t")
    pairs = pd.read_csv(OUT / "lookup_assignment.tsv", sep="\t", usecols=["dataset", "gene_symbol", "actual"])
    finite_genes = pairs[np.isfinite(pairs["actual"])].groupby("dataset")["gene_symbol"].apply(set).to_dict()
    rows = []
    for (dataset, kind), sub in fill.groupby(["dataset", "fill"], sort=True):
        genes = set(sub["gene_symbol"])
        rows.append(
            {
                "label": LABEL,
                "dataset": dataset,
                "fill": kind,
                "n_genes_including_undefined_actual": int(len(genes)),
                "n_genes_with_finite_actual": int(len(genes & finite_genes.get(dataset, set()))),
            }
        )
    write("f0_fill_both_counts.tsv", pd.DataFrame(rows))


def main():
    frame = load_frames()
    rho, l2, gmean, n_added, n_g = lookup_rho(set(frame["gene_symbol"]))
    pairs = assign_lookup(frame, rho, l2, gmean)
    write("lookup_assignment.tsv", pairs)
    write(
        "lookup_build.tsv",
        pd.DataFrame(
            [
                {
                    "label": LABEL,
                    "global_mean": gmean,
                    "n_p_gene_rho_values": n_g,
                    "n_extra_cptac_rho_pairs": n_added,
                    "plan": "config/analysis_plan_B8.yaml",
                }
            ]
        ),
    )
    use = pairs[np.isfinite(pairs["actual"])].copy()
    metric_rows(use)
    auroc_rows(use)
    within_context(use)
    both_f0_counts()
    print("stage8 lookup done", flush=True)


if __name__ == "__main__":
    main()
