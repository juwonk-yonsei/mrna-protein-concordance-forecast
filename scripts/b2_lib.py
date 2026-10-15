"""Shared stage-2 helpers. Importing this module removes the stage-1 CPU cap."""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.common import ROOT

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
RES = ROOT / "results"
PROC = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
CANCERS = PLAN["evaluation"]["cancers"]
MIN_N = int(PLAN["evaluation"]["min_n"])
INNER_N = int(PLAN["evaluation"]["inner_min_n"])
DDOF = int(PLAN["evaluation"]["std_ddof"])
K_EXTRA = int(PLAN["ar"]["extra_genes"])
ALPHAS = [float(a) for a in PLAN["ar"]["alphas"]]
LAMBDAS = [float(v) for v in PLAN["ar"]["lambdas"]]
FEATURE_NAMES = list(PLAN["reliability"]["features"])


def spearman(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < INNER_N:
        return np.nan
    x = a[m]
    y = b[m]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    rx = rankdata(x)
    ry = rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if den == 0:
        return np.nan
    return float((rx * ry).sum() / den)


def zscore_cols(arr):
    arr = np.asarray(arr, dtype=np.float64)
    with np.errstate(all="ignore"):
        mu = np.nanmean(arr, axis=0)
        sd = np.nanstd(arr, axis=0, ddof=DDOF)
    z = arr - mu
    bad = (~np.isfinite(sd)) | (sd == 0)
    sd = sd.copy()
    sd[bad] = 1.0
    z = z / sd
    z[:, bad] = 0.0
    z[~np.isfinite(arr)] = np.nan
    return z.astype(np.float32)


def fit_slope(x_self, y):
    """No-intercept least squares. x_self, y are (n, p)."""
    m = np.isfinite(y) & np.isfinite(x_self)
    prod = np.where(m, x_self * y, 0.0)
    sq = np.where(m, x_self * x_self, 0.0)
    num = prod.sum(axis=0)
    den = sq.sum(axis=0)
    b = np.divide(num, den, out=np.zeros(y.shape[1], dtype=np.float64), where=den > 0)
    b[~np.isfinite(b)] = 0.0
    b[b < 0] = 0.0
    return b


def residuals(x_self, y, b):
    e = y - b.reshape(1, -1) * x_self
    e[~np.isfinite(y)] = np.nan
    return e


def select_residual_features(x, e, self_index):
    """50 RNA columns with largest |Pearson| to the residual. Self excluded."""
    y0 = np.nan_to_num(e.astype(np.float64), nan=0.0)
    xc = x.astype(np.float64)
    xc = xc - xc.mean(axis=0, keepdims=True)
    yc = y0 - y0.mean(axis=0, keepdims=True)
    sx = np.sqrt((xc * xc).sum(axis=0))
    sy = np.sqrt((yc * yc).sum(axis=0))
    denom = sx[:, None] * sy[None, :]
    corr = np.zeros((x.shape[1], e.shape[1]), dtype=np.float64)
    np.divide(xc.T @ yc, denom, out=corr, where=denom > 0)
    score = np.abs(corr)
    score = score + (1e-12 * (1.0 - np.arange(score.shape[0]) / max(score.shape[0], 1)))[:, None]
    score[~np.isfinite(score)] = -np.inf
    for j, gi in enumerate(self_index):
        score[gi, j] = -np.inf
    k = min(K_EXTRA, score.shape[0] - 1)
    picked = np.argpartition(-score, kth=k - 1, axis=0)[:k, :]
    feats = np.empty((k, e.shape[1]), dtype=np.int32)
    for j in range(e.shape[1]):
        cols = picked[:, j]
        order = np.argsort(-score[cols, j])
        feats[:, j] = cols[order]
    return feats


def score_chunk(payload):
    """Inner-test Spearman for every alpha and lambda. BLAS stays at 1 thread."""
    xself_tr, xf_tr, e_tr, xself_te, xf_te, y_te, b = payload
    n_prot = b.shape[0]
    acc = np.zeros((len(ALPHAS), len(LAMBDAS), n_prot), dtype=np.float64)
    cnt = np.zeros((len(ALPHAS), len(LAMBDAS), n_prot), dtype=np.int32)
    with threadpool_limits(limits=1):
        for j in range(n_prot):
            row_ok = np.isfinite(xf_tr[j]).all(axis=1)
            mtr = np.isfinite(e_tr[j]) & np.isfinite(xself_tr[j]) & row_ok
            if int(mtr.sum()) < 1:
                continue
            mte = np.isfinite(y_te[j]) & np.isfinite(xself_te[j]) & np.isfinite(xf_te[j]).all(axis=1)
            if int(mte.sum()) < INNER_N:
                continue
            xtr = np.asarray(xf_tr[j, mtr], dtype=np.float64)
            etr = np.asarray(e_tr[j, mtr], dtype=np.float64)
            xte = np.asarray(xf_te[j, mte], dtype=np.float64)
            yj = np.asarray(y_te[j, mte], dtype=np.float64)
            base = float(b[j]) * np.asarray(xself_te[j, mte], dtype=np.float64)
            k = xtr.shape[1]
            gram0 = xtr.T @ xtr
            xty = xtr.T @ etr
            for a_i, alpha in enumerate(ALPHAS):
                gram = gram0.copy()
                gram.flat[:: k + 1] += alpha
                try:
                    coef = np.linalg.solve(gram, xty)
                except np.linalg.LinAlgError:
                    continue
                ghat = xte @ coef
                for l_i, lam in enumerate(LAMBDAS):
                    rho = spearman(base + lam * ghat, yj)
                    if np.isfinite(rho):
                        acc[a_i, l_i, j] += rho
                        cnt[a_i, l_i, j] += 1
    return acc, cnt


def choose_alpha_lambda(acc, cnt):
    mean = np.full_like(acc, np.nan)
    np.divide(acc, cnt, out=mean, where=cnt > 0)
    n_prot = acc.shape[2]
    alpha = np.empty(n_prot, dtype=np.float64)
    lam = np.empty(n_prot, dtype=np.float64)
    fallback = 0
    fb_l = float(PLAN["ar"]["fallback_if_no_valid_inner_fold"]["lambda"])
    fb_a = float(PLAN["ar"]["fallback_if_no_valid_inner_fold"]["alpha"])
    for j in range(n_prot):
        block = mean[:, :, j]
        if not np.isfinite(block).any():
            alpha[j] = fb_a
            lam[j] = fb_l
            fallback += 1
            continue
        best = np.nanmax(block)
        found = False
        for l_i, lam_v in enumerate(LAMBDAS):
            for a_i in range(len(ALPHAS) - 1, -1, -1):
                val = block[a_i, l_i]
                if np.isfinite(val) and val == best:
                    alpha[j] = ALPHAS[a_i]
                    lam[j] = lam_v
                    found = True
                    break
            if found:
                break
        if not found:
            alpha[j] = fb_a
            lam[j] = fb_l
            fallback += 1
    return alpha, lam, fallback


def predict_ar(x_tr, y_tr, x_te, self_index, feats, alphas, lams):
    xself_tr = x_tr[:, self_index]
    b = fit_slope(xself_tr, y_tr)
    e = residuals(xself_tr, y_tr, b)
    n_te = x_te.shape[0]
    pred = np.full((n_te, y_tr.shape[1]), np.nan, dtype=np.float64)
    with threadpool_limits(limits=1):
        for j in range(y_tr.shape[1]):
            cols = feats[:, j]
            mtr = np.isfinite(e[:, j]) & np.isfinite(x_tr[:, cols]).all(axis=1)
            if int(mtr.sum()) < 1:
                continue
            xtr = np.asarray(x_tr[mtr][:, cols], dtype=np.float64)
            etr = np.asarray(e[mtr, j], dtype=np.float64)
            k = xtr.shape[1]
            gram = xtr.T @ xtr
            gram.flat[:: k + 1] += alphas[j]
            try:
                coef = np.linalg.solve(gram, xtr.T @ etr)
            except np.linalg.LinAlgError:
                continue
            ghat = np.asarray(x_te[:, cols], dtype=np.float64) @ coef
            pred[:, j] = b[j] * np.asarray(x_te[:, self_index[j]], dtype=np.float64) + lams[j] * ghat
    return pred, b


def conformal_half_width(abs_resid, level):
    r = np.sort(np.asarray(abs_resid, dtype=np.float64))
    r = r[np.isfinite(r)]
    n = r.size
    if n == 0:
        return np.inf
    alpha = float(PLAN["reliability"]["intervals"][int(level)])
    k = int(np.ceil((n + 1) * (1.0 - alpha)))
    if k > n:
        return np.inf
    return float(r[k - 1])


def self_check():
    rng = np.random.default_rng(20261001)
    x = rng.normal(size=(20, 1))
    y = 0.5 * x
    y[0, 0] = np.nan
    b = fit_slope(x, y)
    if abs(b[0] - 0.5) > 1e-8:
        raise RuntimeError("slope self-check failed")
    b_neg = fit_slope(x, -x)
    if b_neg[0] != 0:
        raise RuntimeError("negative slope was not clipped")
    # Tie: equal scores, smaller lambda, then larger alpha.
    acc = np.zeros((len(ALPHAS), len(LAMBDAS), 1))
    cnt = np.ones_like(acc)
    acc[:] = 0.2
    alpha, lam, fb = choose_alpha_lambda(acc, cnt)
    if lam[0] != LAMBDAS[0] or alpha[0] != ALPHAS[-1] or fb != 0:
        raise RuntimeError("tie break self-check failed")
    width = conformal_half_width(np.arange(10, dtype=float), 80)
    # alpha 0.2, k = ceil(11*0.8) = 9, 0-based index 8, value 8
    if width != 8:
        raise RuntimeError(f"conformal self-check failed: {width}")
    print("b2 self-check ok", flush=True)
