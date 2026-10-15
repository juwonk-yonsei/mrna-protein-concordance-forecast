"""Development-only power and CPTAC LUAD reference values for the B6 preregistration.

Does not read APOLLO protein values and does not build APOLLO predictions.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

from scripts.common import ROOT

RES = ROOT / "results"
PREP = RES / "stage6" / "prep"
N_REP = 10000
ALPHA = 0.025
SEED = 20261001


def batch_pvalues(source, idx):
    out = np.empty(idx.shape[0], dtype=np.float64)
    with threadpool_limits(limits=1):
        for i in range(idx.shape[0]):
            sample = source[idx[i]]
            try:
                out[i] = float(wilcoxon(sample, alternative="greater", zero_method="wilcox").pvalue)
            except ValueError:
                out[i] = np.nan
    return out


def rejection_rate(source, n, rng):
    source = np.asarray(source, dtype=np.float64)
    idx = rng.integers(0, source.size, size=(N_REP, n))
    n_jobs = os.cpu_count() or 1
    chunks = np.array_split(np.arange(N_REP), n_jobs)
    parts = Parallel(n_jobs=n_jobs, backend="loky")(delayed(batch_pvalues)(source, idx[chunk]) for chunk in chunks)
    pvalues = np.concatenate(parts)
    errors = int(np.sum(~np.isfinite(pvalues)))
    rejects = int(np.sum(pvalues < ALPHA))
    return rejects / N_REP, errors


def scenario_rows(name, source, n, rng):
    source = np.asarray(source, dtype=np.float64)
    source = source[np.isfinite(source)]
    half = float(source.mean()) / 2.0
    shifted = source - half
    rows = []
    for label, vector in (("as_is", source), ("half_effect", shifted)):
        rate, errors = rejection_rate(vector, n, rng)
        rows.append(
            {
                "hypothesis": name,
                "scenario": label,
                "n_source": int(source.size),
                "n_draw": n,
                "source_mean": float(vector.mean()),
                "source_median": float(np.median(vector)),
                "rejection_rate": rate,
                "n_wilcoxon_errors": errors,
                "threshold": ALPHA,
            }
        )
    noise = rng.normal(0.0, 0.05, size=shifted.size)
    noisy = shifted + noise
    rate, errors = rejection_rate(noisy, n, rng)
    rows.append(
        {
            "hypothesis": name,
            "scenario": "half_effect_plus_noise",
            "n_source": int(source.size),
            "n_draw": n,
            "source_mean": float(noisy.mean()),
            "source_median": float(np.median(noisy)),
            "rejection_rate": rate,
            "n_wilcoxon_errors": errors,
            "threshold": ALPHA,
        }
    )
    return rows


def p2_reference():
    pred = pd.read_csv(RES / "f1s_predictions.tsv", sep="\t")
    sub = pred[pred["heldout_cancer"] == "LUAD"].copy()
    use = sub[np.isfinite(sub["actual"]) & np.isfinite(sub["f1s_lgbm"]) & np.isfinite(sub["holdout_rna_sd_percentile"])]
    actual = use["actual"].to_numpy(dtype=np.float64)
    hat = use["f1s_lgbm"].to_numpy(dtype=np.float64)
    base = use["holdout_rna_sd_percentile"].to_numpy(dtype=np.float64)

    def spearman(a, b):
        if np.unique(a).size < 2 or np.unique(b).size < 2:
            return np.nan
        ra = pd.Series(a).rank(method="average").to_numpy()
        rb = pd.Series(b).rank(method="average").to_numpy()
        ra = ra - ra.mean()
        rb = rb - rb.mean()
        den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
        if den == 0:
            return np.nan
        return float((ra * rb).sum() / den)

    s_f = spearman(hat, actual)
    s_b = spearman(base, actual)
    rng = np.random.default_rng(SEED)
    draws_f = np.empty(2000)
    draws_b = np.empty(2000)
    n = actual.size
    for i in range(2000):
        chosen = rng.integers(0, n, size=n)
        draws_f[i] = spearman(hat[chosen], actual[chosen])
        draws_b[i] = spearman(base[chosen], actual[chosen])
    diff = draws_f - draws_b

    def ci(v):
        v = v[np.isfinite(v)]
        return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))

    return {
        "n_proteins": int(n),
        "S_F": s_f,
        "S_F_ci_low": ci(draws_f)[0],
        "S_F_ci_high": ci(draws_f)[1],
        "S_B": s_b,
        "S_B_ci_low": ci(draws_b)[0],
        "S_B_ci_high": ci(draws_b)[1],
        "diff": s_f - s_b,
        "diff_ci_low": ci(diff)[0],
        "diff_ci_high": ci(diff)[1],
        "n_nonfinite_draws": int(np.sum(~np.isfinite(diff))),
    }


def main():
    print("thread env", {k: os.environ.get(k, "") for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, flush=True)
    ar = pd.read_csv(RES / "ar_pairs.tsv", sep="\t")
    loco = pd.read_csv(RES / "loco_pairs.tsv", sep="\t")
    luad = ar[ar["heldout_cancer"] == "LUAD"].copy()
    loco_luad = loco.loc[loco["heldout_cancer"] == "LUAD", ["gene_symbol", "r1", "r2"]]
    joined = luad.merge(loco_luad, on="gene_symbol", how="left", suffixes=("", "_loco"))
    both = joined[np.isfinite(joined["r1"]) & np.isfinite(joined["r1_loco"])]
    if not np.allclose(both["r1"], both["r1_loco"]):
        raise RuntimeError("ar_pairs r1 does not match loco_pairs r1")
    p1 = luad.loc[np.isfinite(luad["delta_ar"]), "delta_ar"].to_numpy()
    q1_rows = joined[np.isfinite(joined["ar"]) & np.isfinite(joined["r1_loco"])]
    q1 = (q1_rows["ar"] - q1_rows["r1_loco"]).to_numpy()
    q2_rows = joined[np.isfinite(joined["ar"]) & np.isfinite(joined["r2"])]
    q2 = (q2_rows["ar"] - q2_rows["r2"]).to_numpy()
    n_p1 = sum(1 for line in (PREP / "U_P1_genes.tsv").read_text().splitlines() if line and line != "gene_symbol")
    n_q = sum(1 for line in (PREP / "U_Q_genes.tsv").read_text().splitlines() if line and line != "gene_symbol")
    rng = np.random.default_rng(SEED)
    rows = []
    rows.extend(scenario_rows("P1", p1, n_p1, rng))
    rows.extend(scenario_rows("Q1", q1, n_q, rng))
    rows.extend(scenario_rows("Q2", q2, n_q, rng))
    power = pd.DataFrame(rows)
    power.to_csv(PREP / "power.tsv", sep="\t", index=False)
    ref = pd.DataFrame([p2_reference()])
    ref.to_csv(PREP / "p2_cptac_luad_reference.tsv", sep="\t", index=False)
    print(power.to_string(index=False), flush=True)
    print(ref.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
