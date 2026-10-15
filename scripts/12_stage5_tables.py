"""Stage-5 post-hoc tables that use only the stage-4 confirmation pair files.

Does not read locked matrices, paths.yaml, or n_threads.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.b2_lib import RES, ROOT
from scripts.common import sha256_file

LABEL = "사후"
SEED = 20261014
BOOT_SEED = 20261001
BOOT_N = 2000
WITHIN_MIN = 10
KS = (10, 20, 50, 100)
N_REPEATS = 20
N_FOLDS = 5
OUT = RES / "stage5" / "posthoc"
PREDICTORS = ("f1s_lgbm", "f1s_ridge")


def finite_spearman(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < 3:
        return np.nan
    x = a[m]
    y = b[m]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if den == 0:
        return np.nan
    return float((rx * ry).sum() / den)


def write(path, frame):
    frame = frame.copy()
    frame.insert(0, "label", LABEL)
    frame.to_csv(path, sep="\t", index=False)
    print("wrote", path, len(frame), flush=True)


def load_pairs():
    rppa = pd.read_csv(RES / "stage4" / "confirm" / "rppa_pairs.tsv", sep="\t")
    pro = pd.read_csv(RES / "stage4" / "confirm" / "procan_pairs.tsv", sep="\t")
    return {"rppa": ("xena_antibody", "cancer", rppa), "procan": ("protein_id", "lineage", pro)}


def bias_table(bundles):
    rows = []
    for platform, (_unit, context, frame) in bundles.items():
        for pred in PREDICTORS:
            m = np.isfinite(frame["actual"]) & np.isfinite(frame[pred])
            err = frame.loc[m, "actual"].to_numpy() - frame.loc[m, pred].to_numpy()
            ctx = []
            for _name, sub in frame.loc[m].groupby(context, sort=True):
                d = sub["actual"].to_numpy() - sub[pred].to_numpy()
                d = d[np.isfinite(d)]
                if d.size:
                    ctx.append(float(d.mean()))
            ctx = np.asarray(ctx, dtype=np.float64)
            rows.append(
                {
                    "platform": platform,
                    "predictor": pred,
                    "n_pairs": int(err.size),
                    "mean_actual_minus_pred": float(err.mean()) if err.size else np.nan,
                    "n_contexts": int(ctx.size),
                    "context_mean_sd_ddof0": float(ctx.std(ddof=0)) if ctx.size else np.nan,
                }
            )
    write(OUT / "posthoc_24_bias.tsv", pd.DataFrame(rows))


def calibration(bundles):
    rng = np.random.default_rng(SEED)
    draws = []
    skips = []
    for platform, (unit, _context, frame) in bundles.items():
        for pred in PREDICTORS:
            half80 = f"{pred}_half_80"
            half90 = f"{pred}_half_90"
            ok = frame[np.isfinite(frame["actual"]) & np.isfinite(frame[pred])].copy()
            units = np.array(sorted(ok[unit].unique()))
            rng.shuffle(units)
            folds = np.array_split(units, N_FOLDS)
            for fold_i, eval_units in enumerate(folds):
                eval_set = set(eval_units.tolist())
                train_units = np.array([u for u in units if u not in eval_set])
                eval_rows = ok[ok[unit].isin(eval_set)]
                for k in KS:
                    if train_units.size < k:
                        skips.append(
                            {
                                "platform": platform,
                                "predictor": pred,
                                "fold": fold_i,
                                "K": k,
                                "n_train_units": int(train_units.size),
                                "n_skipped_reps": N_REPEATS,
                            }
                        )
                        continue
                    for rep in range(N_REPEATS):
                        chosen = rng.choice(train_units, size=k, replace=False)
                        cal = ok[ok[unit].isin(set(chosen.tolist()))]
                        actual = cal["actual"].to_numpy(dtype=np.float64)
                        hat = cal[pred].to_numpy(dtype=np.float64)
                        delta = float((actual - hat).mean())
                        shifted = hat + delta
                        abs_err = np.abs(actual - shifted)
                        q80 = float(np.quantile(abs_err, 0.80))
                        q90 = float(np.quantile(abs_err, 0.90))
                        ev_a = eval_rows["actual"].to_numpy(dtype=np.float64)
                        ev_p = eval_rows[pred].to_numpy(dtype=np.float64)
                        h80 = eval_rows[half80].to_numpy(dtype=np.float64)
                        h90 = eval_rows[half90].to_numpy(dtype=np.float64)
                        use80 = np.isfinite(ev_a) & np.isfinite(ev_p) & np.isfinite(h80)
                        use90 = np.isfinite(ev_a) & np.isfinite(ev_p) & np.isfinite(h90)
                        sp = finite_spearman(ev_p, ev_a)
                        ss = finite_spearman(ev_p + delta, ev_a)
                        draws.append(
                            {
                                "platform": platform,
                                "predictor": pred,
                                "fold": fold_i,
                                "K": k,
                                "rep": rep,
                                "n_train_units": int(train_units.size),
                                "n_cal_pairs": int(len(cal)),
                                "n_eval_pairs": int(len(eval_rows)),
                                "delta": delta,
                                "q80": q80,
                                "q90": q90,
                                "cov80_none": float(np.mean(np.abs(ev_a[use80] - ev_p[use80]) <= h80[use80])),
                                "cov90_none": float(np.mean(np.abs(ev_a[use90] - ev_p[use90]) <= h90[use90])),
                                "cov80_bias": float(np.mean(np.abs(ev_a[use80] - ev_p[use80] - delta) <= h80[use80])),
                                "cov90_bias": float(np.mean(np.abs(ev_a[use90] - ev_p[use90] - delta) <= h90[use90])),
                                "cov80_both": float(np.mean(np.abs(ev_a[use80] - ev_p[use80] - delta) <= q80)),
                                "cov90_both": float(np.mean(np.abs(ev_a[use90] - ev_p[use90] - delta) <= q90)),
                                "width80_none": float(np.mean(2.0 * h80[use80])),
                                "width90_none": float(np.mean(2.0 * h90[use90])),
                                "width80_bias": float(np.mean(2.0 * h80[use80])),
                                "width90_bias": float(np.mean(2.0 * h90[use90])),
                                "width80_both": float(2.0 * q80),
                                "width90_both": float(2.0 * q90),
                                "spearman_pred": sp,
                                "spearman_shifted": ss,
                                "spearman_abs_diff": abs(sp - ss) if np.isfinite(sp) and np.isfinite(ss) else np.nan,
                            }
                        )
    draw_df = pd.DataFrame(draws)
    write(OUT / "posthoc_24_draws.tsv", draw_df)
    if skips:
        write(OUT / "posthoc_24_skips.tsv", pd.DataFrame(skips))
    value_cols = [
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
        "delta",
        "spearman_abs_diff",
    ]
    summary = draw_df.groupby(["platform", "predictor", "K"], sort=True)[value_cols].mean().reset_index()
    counts = draw_df.groupby(["platform", "predictor", "K"], sort=True).size().rename("n_draws").reset_index()
    summary = summary.merge(counts, on=["platform", "predictor", "K"])
    write(OUT / "posthoc_24_summary.tsv", summary)
    print("max spearman shift", float(draw_df["spearman_abs_diff"].max()), flush=True)


def within_antibody(rppa):
    rows = []
    for pred in list(PREDICTORS) + ["b_sd"]:
        for antibody, sub in rppa.groupby("xena_antibody", sort=True):
            m = np.isfinite(sub[pred]) & np.isfinite(sub["actual"])
            n = int(m.sum())
            if n < WITHIN_MIN:
                continue
            rows.append(
                {
                    "xena_antibody": antibody,
                    "gene_symbol": sub["gene_symbol"].iloc[0],
                    "predictor": pred,
                    "n_cancers": n,
                    "spearman": finite_spearman(sub.loc[m, pred], sub.loc[m, "actual"]),
                }
            )
    detail = pd.DataFrame(rows)
    write(OUT / "posthoc_25_within_antibody.tsv", detail)
    wide = detail.pivot(index="xena_antibody", columns="predictor", values="spearman")
    summaries = []
    for pred in PREDICTORS:
        if pred not in wide.columns or "b_sd" not in wide.columns:
            continue
        both = wide[[pred, "b_sd"]].dropna()
        diff = (both[pred] - both["b_sd"]).to_numpy(dtype=np.float64)
        rng = np.random.default_rng(BOOT_SEED)
        drawn = rng.choice(diff, size=(BOOT_N, diff.size), replace=True)
        medians = np.median(drawn, axis=1)
        summaries.append(
            {
                "predictor": pred,
                "n_antibodies": int(diff.size),
                "median_diff": float(np.median(diff)),
                "ci_low": float(np.quantile(medians, 0.025)),
                "ci_high": float(np.quantile(medians, 0.975)),
                "bootstrap_draws": BOOT_N,
                "bootstrap_seed": BOOT_SEED,
            }
        )
    write(OUT / "posthoc_25_diff.tsv", pd.DataFrame(summaries))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    caps = {k: os.environ.get(k, "") for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_MAX_THREADS")}
    print("thread env", caps, "B_CPU_CAP", os.environ.get("B_CPU_CAP"), flush=True)
    bundles = load_pairs()
    bias_table(bundles)
    within_antibody(bundles["rppa"][2])
    calibration(bundles)
    print("sha rppa", sha256_file(RES / "stage4" / "confirm" / "rppa_pairs.tsv"), flush=True)
    print("sha procan", sha256_file(RES / "stage4" / "confirm" / "procan_pairs.tsv"), flush=True)


if __name__ == "__main__":
    main()
