"""Independent check of a B6 confirmation directory.

Uses pandas, numpy, and scipy only. Does not import scripts/17_confirm_B6.py.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


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


def finite_spearman(a, b, min_n):
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


def wilcoxon_greater(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    result = wilcoxon(values, alternative="greater", zero_method="wilcox")
    return int(values.size), int(np.sum(values == 0)), float(result.statistic), float(result.pvalue)


def mean_median_ci(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(20261001)
    drawn = rng.choice(values, size=(2000, values.size), replace=True)
    means = drawn.mean(axis=1)
    medians = np.median(drawn, axis=1)
    return {
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "mean_ci_low": float(np.quantile(means, 0.025)),
        "mean_ci_high": float(np.quantile(means, 0.975)),
        "median_ci_low": float(np.quantile(medians, 0.025)),
        "median_ci_high": float(np.quantile(medians, 0.975)),
    }


def paired_spearman(actual, hat, base):
    s_f = finite_spearman(hat, actual, 3)
    s_b = finite_spearman(base, actual, 3)
    rng = np.random.default_rng(20261001)
    n = actual.size
    draws_f = np.empty(2000)
    draws_b = np.empty(2000)
    for i in range(2000):
        chosen = rng.integers(0, n, size=n)
        draws_f[i] = finite_spearman(hat[chosen], actual[chosen], 3)
        draws_b[i] = finite_spearman(base[chosen], actual[chosen], 3)
    diff = draws_f - draws_b
    ok = np.isfinite(diff)
    return {
        "S_F": s_f,
        "S_B": s_b,
        "diff": s_f - s_b,
        "pvalue": (1.0 + int(np.sum(diff[ok] <= 0))) / 2001.0,
        "diff_ci_low": float(np.quantile(diff[ok], 0.025)),
        "diff_ci_high": float(np.quantile(diff[ok], 0.975)),
        "n": int(n),
    }


def holm(rows, order):
    ranked = sorted(rows, key=lambda row: ((row["pvalue"] if np.isfinite(row["pvalue"]) else np.inf), order.index(row["hypothesis"])))
    stopped = False
    calls = {}
    n = len(order)
    for k, row in enumerate(ranked, start=1):
        threshold = 0.05 / (n - k + 1)
        fails = (not np.isfinite(row["pvalue"])) or row["pvalue"] > threshold
        if stopped or fails:
            stopped = True
            calls[row["hypothesis"]] = "not_rejected"
        else:
            calls[row["hypothesis"]] = "rejected"
    return calls


def close(a, b):
    if not np.isfinite(a) and not np.isfinite(b):
        return True
    return abs(float(a) - float(b)) <= 1e-9


def main():
    dest = Path(sys.argv[1])
    with np.load(dest / "protein.npz", allow_pickle=False) as protein_raw:
        protein_genes = np.array(protein_raw["genes"])
        protein_values = np.array(protein_raw["values"])
    with np.load(dest / "used_predictions.npz", allow_pickle=False) as used_raw:
        used = {key: np.array(used_raw[key]) for key in used_raw.files}
    level = pd.read_csv(dest / "gene_level.tsv", sep="\t").set_index("gene_symbol")
    saved = pd.read_csv(dest / "genes.tsv", sep="\t")
    holm_saved = pd.read_csv(dest / "holm.tsv", sep="\t")
    if not np.array_equal(protein_genes, used["genes"]):
        raise RuntimeError("gene order differs")
    rows = []
    for i, symbol in enumerate(used["genes"]):
        yy = protein_values[i]
        actual = finite_spearman(used["log2"][i], yy, 30)
        rho_s0 = finite_spearman(used["rna_z"][i], yy, 30)
        rho_ar = finite_spearman(used["ar"][i], yy, 30)
        rho_r1 = finite_spearman(used["r1"][i], yy, 30)
        rho_r2 = finite_spearman(used["r2"][i], yy, 30)
        rows.append(
            {
                "gene_symbol": str(symbol),
                "actual": actual,
                "d_p1": rho_ar - rho_s0 if np.isfinite(rho_ar) and np.isfinite(rho_s0) else np.nan,
                "d_q1": rho_ar - rho_r1 if np.isfinite(rho_ar) and np.isfinite(rho_r1) else np.nan,
                "d_q2": rho_ar - rho_r2 if np.isfinite(rho_ar) and np.isfinite(rho_r2) else np.nan,
            }
        )
    fresh = pd.DataFrame(rows)
    merged = saved.merge(fresh, on="gene_symbol", suffixes=("_saved", "_fresh"))
    worst = 0.0
    for col in ("actual", "d_p1", "d_q1", "d_q2"):
        a = merged[f"{col}_saved"].to_numpy(dtype=np.float64)
        b = merged[f"{col}_fresh"].to_numpy(dtype=np.float64)
        gap = np.nanmax(np.abs(a - b))
        worst = max(worst, float(gap))
    if worst > 1e-9:
        raise RuntimeError(f"rho gap {worst}")
    p1 = fresh.loc[saved["in_U_P1"].to_numpy() == 1, "d_p1"].to_numpy()
    q1 = fresh.loc[saved["in_U_Q"].to_numpy() == 1, "d_q1"].to_numpy()
    q2 = fresh.loc[saved["in_U_Q"].to_numpy() == 1, "d_q2"].to_numpy()
    mask = (saved["in_U_P2"].to_numpy() == 1) & np.isfinite(fresh["actual"]) & fresh["gene_symbol"].isin(level.index)
    sub = fresh.loc[mask].merge(level, left_on="gene_symbol", right_index=True)
    sub = sub[np.isfinite(sub["f1s_lgbm"]) & np.isfinite(sub["b_sd"])]
    p2 = paired_spearman(sub["actual"].to_numpy(), sub["f1s_lgbm"].to_numpy(), sub["b_sd"].to_numpy())
    checks = []
    for name, values in (("P1", p1), ("Q1", q1), ("Q2", q2)):
        n, _z, stat, pvalue = wilcoxon_greater(values)
        effect = mean_median_ci(values)
        saved_row = pd.read_csv(dest / "hypotheses.tsv", sep="\t")
        saved_row = saved_row[saved_row["hypothesis"] == name].iloc[0]
        for key, value in (("n", n), ("statistic", stat), ("pvalue", pvalue), ("mean", effect["mean"]), ("median", effect["median"])):
            if not close(saved_row[key], value):
                raise RuntimeError(f"{name} {key} {saved_row[key]} != {value}")
        checks.append({"hypothesis": name, "pvalue": pvalue})
    saved_p2 = pd.read_csv(dest / "hypotheses.tsv", sep="\t")
    saved_p2 = saved_p2[saved_p2["hypothesis"] == "P2"].iloc[0]
    for key in ("n", "S_F", "S_B", "diff", "pvalue"):
        if not close(saved_p2[key], p2[key]):
            raise RuntimeError(f"P2 {key} {saved_p2[key]} != {p2[key]}")
    calls = holm([{"hypothesis": "P1", "pvalue": checks[0]["pvalue"]}, {"hypothesis": "P2", "pvalue": p2["pvalue"]}], ["P1", "P2"])
    calls.update(holm([{"hypothesis": "Q1", "pvalue": checks[1]["pvalue"]}, {"hypothesis": "Q2", "pvalue": checks[2]["pvalue"]}], ["Q1", "Q2"]))
    for rec in holm_saved.itertuples(index=False):
        if calls[rec.hypothesis] != rec.call:
            raise RuntimeError(f"Holm {rec.hypothesis} {rec.call} != {calls[rec.hypothesis]}")
    print("check ok", dest, "rho_gap", worst, flush=True)


if __name__ == "__main__":
    main()
