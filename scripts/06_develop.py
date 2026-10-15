"""Section 3. AR and RNA-proxy reliability. CPTAC only."""
import os
import sys
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import time

import numpy as np
import pandas as pd
from lightgbm import LGBMRegressor
from scipy.stats import rankdata
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score

from scripts.b2_lib import (
    ALPHAS,
    CANCERS,
    FEATURE_NAMES,
    INNER_N,
    K_EXTRA,
    LAMBDAS,
    MIN_N,
    PLAN,
    PROC,
    RES,
    choose_alpha_lambda,
    conformal_half_width,
    fit_slope,
    predict_ar,
    residuals,
    score_chunk,
    select_residual_features,
    self_check,
    spearman,
    zscore_cols,
)


def load_targets():
    genes = [g for g in (RES / "G_RNA_genes.txt").read_text().split() if g]
    p_genes = [g for g in (RES / "P_genes.txt").read_text().split() if g]
    adc = pd.read_csv(RES.parent / "config" / "adc_targets.tsv", sep="\t", dtype=str)["gene_symbol"].tolist()
    gene_index = {g: i for i, g in enumerate(genes)}
    extra = [g for g in adc if g not in set(p_genes) and g in gene_index]
    missing = [g for g in adc if g not in gene_index]
    targets = p_genes + extra
    return genes, p_genes, adc, missing, targets, gene_index


def load_tumor(genes, targets):
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "tumor"]
    data = {}
    rna_stats = []
    all_ids = []
    for cancer in CANCERS:
        rna = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__rna.parquet")
        prot = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__protein.parquet")
        sids = paired.loc[paired["cancer"] == cancer, "sample_id"].tolist()
        sids = [s for s in sids if s in rna.columns and s in prot.columns]
        raw = rna.reindex(index=genes, columns=sids).to_numpy(dtype=np.float64)
        with np.errstate(all="ignore"):
            mu = np.nanmean(raw, axis=1)
            sd = np.nanstd(raw, axis=1, ddof=0)
        for i, gene in enumerate(genes):
            rna_stats.append({"cancer": cancer, "gene_symbol": gene, "rna_mean": mu[i], "rna_sd": sd[i]})
        x = zscore_cols(raw.T)
        y_raw = prot.reindex(index=targets, columns=sids).to_numpy(dtype=np.float64).T
        data[cancer] = {"sample_id": sids, "x": x, "y": zscore_cols(y_raw)}
        all_ids.extend(sids)
        print(f"loaded {cancer} n={len(sids)}", flush=True)
    if len(all_ids) != len(set(all_ids)):
        raise RuntimeError("tumor sample ids overlap across cancers")
    stats = pd.DataFrame(rna_stats)
    stats.to_csv(RES / "cptac_rna_gene_stats.tsv", sep="\t", index=False)
    return data, stats


def stack(data, cancers):
    x = np.vstack([data[c]["x"] for c in cancers])
    y = np.vstack([data[c]["y"] for c in cancers])
    labels = np.concatenate(
        [np.full(data[c]["x"].shape[0], i, dtype=np.int32) for i, c in enumerate(cancers)]
    )
    return x, y, labels


def pack_split(x, y, e, b, feats, self_index, rows, j0, j1):
    idx = np.arange(j0, j1)
    xs = x[rows]
    n = xs.shape[0]
    k = feats.shape[0]
    xf = np.empty((len(idx), n, k), dtype=np.float32)
    xself = np.empty((len(idx), n), dtype=np.float32)
    ee = np.empty((len(idx), n), dtype=np.float32)
    yy = np.empty((len(idx), n), dtype=np.float32)
    for t, j in enumerate(idx):
        xf[t] = xs[:, feats[:, j]]
        xself[t] = xs[:, self_index[j]]
        ee[t] = e[rows, j]
        yy[t] = y[rows, j]
    return xself, xf, ee, yy, np.asarray(b[idx], dtype=np.float64)


def inner_scores(pool, x, y, labels, self_index, n_workers):
    n_prot = y.shape[1]
    acc = np.zeros((len(ALPHAS), len(LAMBDAS), n_prot), dtype=np.float64)
    cnt = np.zeros((len(ALPHAS), len(LAMBDAS), n_prot), dtype=np.int32)
    n_groups = int(labels.max()) + 1
    chunk = max(16, n_prot // (n_workers * 4))
    for g in range(n_groups):
        te = np.flatnonzero(labels == g)
        tr = np.flatnonzero(labels != g)
        if tr.size < 1 or te.size < INNER_N:
            continue
        b = fit_slope(x[tr][:, self_index], y[tr])
        e = residuals(x[tr][:, self_index], y[tr], b)
        feats = select_residual_features(x[tr], e, self_index)
        # Residuals were computed on the training-row subset. Map them back
        # to the full outer-sample axis so pack_split can use row indices.
        e_full = np.full_like(y, np.nan, dtype=np.float32)
        e_full[tr] = e
        payloads = []
        spans = []
        for j0 in range(0, n_prot, chunk):
            j1 = min(n_prot, j0 + chunk)
            tr_self, tr_xf, tr_e, _, b_sub = pack_split(x, y, e_full, b, feats, self_index, tr, j0, j1)
            te_self, te_xf, _, te_y, _ = pack_split(x, y, e_full, b, feats, self_index, te, j0, j1)
            payloads.append((tr_self, tr_xf, tr_e, te_self, te_xf, te_y, b_sub))
            spans.append((j0, j1))
        parts = pool.map(score_chunk, payloads, chunksize=1)
        for (j0, j1), (part_acc, part_cnt) in zip(spans, parts):
            acc[:, :, j0:j1] += part_acc
            cnt[:, :, j0:j1] += part_cnt
    return choose_alpha_lambda(acc, cnt)


def pair_scores(pred, y):
    out = np.full(pred.shape[1], np.nan)
    n_obs = np.zeros(pred.shape[1], dtype=np.int32)
    for j in range(pred.shape[1]):
        m = np.isfinite(pred[:, j]) & np.isfinite(y[:, j])
        n_obs[j] = int(m.sum())
        if n_obs[j] >= MIN_N:
            out[j] = spearman(pred[m, j], y[m, j])
    return out, n_obs


def run_ar(data, targets, p_set, self_index, pool, n_workers):
    stage1 = pd.read_csv(RES / "loco_pairs.tsv", sep="\t")
    adc_stage1 = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    r1_lookup = {}
    for rec in stage1.itertuples(index=False):
        r1_lookup[(rec.heldout_cancer, rec.gene_symbol)] = rec.r1
    for rec in adc_stage1.itertuples(index=False):
        r1_lookup[(rec.heldout_cancer, rec.gene_symbol)] = rec.r1
    rows = []
    for hold in CANCERS:
        t0 = time.perf_counter()
        train = [c for c in CANCERS if c != hold]
        train_ids = {s for c in train for s in data[c]["sample_id"]}
        assert hold not in train and train_ids.isdisjoint(data[hold]["sample_id"])
        x_tr, y_tr, labels = stack(data, train)
        alpha, lam, n_fb = inner_scores(pool, x_tr, y_tr, labels, self_index, n_workers)
        b_full = fit_slope(x_tr[:, self_index], y_tr)
        e_full = residuals(x_tr[:, self_index], y_tr, b_full)
        feats = select_residual_features(x_tr, e_full, self_index)
        pred, _ = predict_ar(x_tr, y_tr, data[hold]["x"], self_index, feats, alpha, lam)
        ar, n_obs = pair_scores(pred, data[hold]["y"])
        s0, _ = pair_scores(data[hold]["x"][:, self_index], data[hold]["y"])
        for j, gene in enumerate(targets):
            r1 = r1_lookup.get((hold, gene), np.nan)
            rows.append(
                {
                    "heldout_cancer": hold,
                    "gene_symbol": gene,
                    "in_P": gene in p_set,
                    "n": int(n_obs[j]),
                    "s0": s0[j],
                    "r1": r1,
                    "ar": ar[j],
                    "b": b_full[j],
                    "alpha": alpha[j],
                    "lambda": lam[j],
                    "fallback": n_fb,
                    "delta_ar": ar[j] - s0[j] if np.isfinite(ar[j]) and np.isfinite(s0[j]) else np.nan,
                    "delta_ar_r1": ar[j] - r1 if np.isfinite(ar[j]) and np.isfinite(r1) else np.nan,
                }
            )
        print(f"AR {hold} {time.perf_counter() - t0:.1f}s fallback_proteins={n_fb}", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "ar_pairs.tsv", sep="\t", index=False)
    return out


def quant_row(name, series):
    s = series.dropna()
    return {
        "metric": name,
        "n": int(s.size),
        "median": float(s.median()) if len(s) else np.nan,
        "q1": float(s.quantile(0.25)) if len(s) else np.nan,
        "q3": float(s.quantile(0.75)) if len(s) else np.nan,
        "frac_gt_0.05": float((s > 0.05).mean()) if len(s) else np.nan,
        "frac_lt_m0.05": float((s < -0.05).mean()) if len(s) else np.nan,
    }


def summarize_ar(pairs, adc):
    ev = pairs[(pairs["in_P"]) & (pairs["n"] >= MIN_N)]
    rows = [quant_row("delta_ar", ev["delta_ar"]), quant_row("delta_ar_r1", ev["delta_ar_r1"])]
    # Stage-1 delta_r1 on the same P pairs, for the summary line.
    stage1 = pd.read_csv(RES / "loco_delta_summary.tsv", sep="\t")
    pd.DataFrame(rows).to_csv(RES / "ar_delta_summary.tsv", sep="\t", index=False)
    by = []
    for hold, sub in ev.groupby("heldout_cancer"):
        by.append(
            {
                "heldout_cancer": hold,
                "n": int(sub.shape[0]),
                "delta_ar_median": float(sub["delta_ar"].median()),
                "delta_ar_r1_median": float(sub["delta_ar_r1"].median()),
                "s0_median": float(sub["s0"].median()),
                "ar_median": float(sub["ar"].median()),
                "r1_median": float(sub["r1"].median()),
            }
        )
    pd.DataFrame(by).to_csv(RES / "ar_by_heldout.tsv", sep="\t", index=False)
    lam_counts = ev["lambda"].value_counts().rename_axis("lambda").reset_index(name="n_pairs")
    lam_counts["parameter"] = "lambda"
    alpha_counts = ev["alpha"].value_counts().rename_axis("alpha").reset_index(name="n_pairs")
    alpha_counts["parameter"] = "alpha"
    lam_counts.to_csv(RES / "ar_lambda_counts.tsv", sep="\t", index=False)
    alpha_counts.to_csv(RES / "ar_alpha_counts.tsv", sep="\t", index=False)
    adc_rows = []
    sub_all = pairs[pairs["gene_symbol"].isin(adc)]
    for gene in adc:
        sub = sub_all[sub_all["gene_symbol"] == gene]
        use = sub[sub["n"] >= MIN_N]
        lams = use["lambda"].dropna()
        if len(lams):
            counts = lams.value_counts()
            mode = float(sorted(counts[counts == counts.max()].index)[0])
        else:
            mode = np.nan
        adc_rows.append(
            {
                "gene_symbol": gene,
                "in_P": bool(sub["in_P"].iloc[0]) if len(sub) else False,
                "n_cancers_n_ge_30": int((sub["n"] >= MIN_N).sum()),
                "s0_median": use["s0"].median(),
                "r1_median": use["r1"].median(),
                "ar_median": use["ar"].median(),
                "lambda_mode": mode,
            }
        )
    pd.DataFrame(adc_rows).to_csv(RES / "ar_adc.tsv", sep="\t", index=False)
    stage1.to_csv(RES / "ar_stage1_delta_reference.tsv", sep="\t", index=False)


def defined_rho():
    s0 = pd.read_csv(RES / "s0_pairs.tsv", sep="\t")
    use = s0[(s0["n"] >= MIN_N) & s0["s0"].notna()]
    rho = {(r.cancer, r.gene_symbol): float(r.s0) for r in use.itertuples(index=False)}
    adc = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    adc_rho = {}
    for r in adc.itertuples(index=False):
        if r.n >= MIN_N and np.isfinite(r.s0):
            adc_rho[(r.heldout_cancer, r.gene_symbol)] = float(r.s0)
            if not r.in_P:
                rho[(r.heldout_cancer, r.gene_symbol)] = float(r.s0)
    return rho, adc_rho


def load_rates():
    rates = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    return {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}


def sd_percentile(stats):
    out = {}
    for cancer, sub in stats.groupby("cancer"):
        sd = sub["rna_sd"].to_numpy(dtype=np.float64)
        genes = sub["gene_symbol"].tolist()
        finite = np.isfinite(sd)
        pct = np.full(len(sd), np.nan)
        if finite.any():
            pct[finite] = 100.0 * rankdata(sd[finite], method="average") / float(finite.sum())
        for gene, val, s in zip(genes, pct, sd):
            out[(cancer, gene)] = (float(s) if np.isfinite(s) else np.nan, float(val) if np.isfinite(val) else np.nan)
    return out


def feature_row(gene, source, query, rho, rna_mean, pct, rate, global_mean):
    vals = [rho[(c, gene)] for c in source if (c, gene) in rho]
    if vals:
        rho_mean = float(np.mean(vals))
        rho_sd = float(np.std(vals, ddof=0)) if len(vals) >= 2 else 0.0
        rho_min = float(np.min(vals))
        rho_max = float(np.max(vals))
        n_def = len(vals)
    else:
        rho_mean = global_mean
        rho_sd = 0.0
        rho_min = global_mean
        rho_max = global_mean
        n_def = 0
    means = [rna_mean.get((c, gene), np.nan) for c in source]
    means = [v for v in means if np.isfinite(v)]
    rates = [rate.get((c, gene), 0.0) for c in source]
    hold_mean = rna_mean.get((query, gene), np.nan)
    hold_sd, hold_pct = pct.get((query, gene), (np.nan, np.nan))
    return {
        "rho_mean": rho_mean,
        "rho_sd": rho_sd,
        "rho_min": rho_min,
        "rho_max": rho_max,
        "n_defined_cancers": n_def,
        "holdout_rna_mean": hold_mean,
        "holdout_rna_sd": hold_sd,
        "holdout_rna_sd_percentile": hold_pct,
        "train_rna_mean": float(np.mean(means)) if means else np.nan,
        "train_protein_rate_mean": float(np.mean(rates)) if rates else np.nan,
    }


def global_mean_rho(source, rho, p_genes):
    vals = [rho[(c, g)] for c in source for g in p_genes if (c, g) in rho]
    return float(np.mean(vals)) if vals else np.nan


def design(genes, source, query, rho, rna_mean, pct, rate, p_genes):
    gmean = global_mean_rho(source, rho, p_genes)
    rows = []
    for gene in genes:
        rec = feature_row(gene, source, query, rho, rna_mean, pct, rate, gmean)
        rec["gene_symbol"] = gene
        rec["query_cancer"] = query
        actual = rho.get((query, gene), np.nan)
        rec["actual"] = actual
        rows.append(rec)
    return pd.DataFrame(rows)


def fill_median(train_x, apply_x):
    med = np.nanmedian(train_x, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)

    def _fill(x):
        out = np.array(x, dtype=np.float64, copy=True)
        inds = np.where(~np.isfinite(out))
        out[inds] = np.take(med, inds[1])
        return out

    return _fill(train_x), _fill(apply_x), med


def standardize(train_x, apply_x):
    mu = train_x.mean(axis=0)
    sd = train_x.std(axis=0, ddof=0)
    sd = np.where((sd == 0) | ~np.isfinite(sd), 1.0, sd)
    bad = (train_x.std(axis=0, ddof=0) == 0) | ~np.isfinite(train_x.std(axis=0, ddof=0))

    def _z(x):
        z = (x - mu) / sd
        z[:, bad] = 0.0
        return z

    return _z(train_x), _z(apply_x)


def ridge_alpha_and_oof(frame):
    alphas = [float(a) for a in PLAN["reliability"]["ridge_alphas"]]
    cancers = list(dict.fromkeys(frame["query_cancer"]))
    x_all = frame[FEATURE_NAMES].to_numpy(dtype=np.float64)
    y_all = frame["actual"].to_numpy(dtype=np.float64)
    groups = [frame["query_cancer"].to_numpy() == c for c in cancers]
    maes = []
    for alpha in alphas:
        per = []
        for i, mask in enumerate(groups):
            tr = ~mask & np.isfinite(y_all)
            te = mask & np.isfinite(y_all)
            if tr.sum() < 1 or te.sum() < 1:
                continue
            xtr, xte, _ = fill_median(x_all[tr], x_all[te])
            ztr, zte = standardize(xtr, xte)
            model = Ridge(alpha=alpha, fit_intercept=True)
            model.fit(ztr, y_all[tr])
            pred = model.predict(zte)
            per.append(float(np.mean(np.abs(pred - y_all[te]))))
        maes.append(float(np.mean(per)) if per else np.nan)
    maes = np.asarray(maes, dtype=np.float64)
    if not np.isfinite(maes).any():
        chosen = float(PLAN["reliability"]["ridge_fallback_alpha"])
    else:
        best = np.nanmin(maes)
        chosen = alphas[0]
        for alpha, val in zip(alphas, maes):
            if np.isfinite(val) and val == best:
                chosen = alpha
                break
    oof = np.full(len(frame), np.nan)
    for mask in groups:
        tr = ~mask & np.isfinite(y_all)
        te = mask & np.isfinite(y_all)
        if tr.sum() < 1 or te.sum() < 1:
            continue
        xtr, xte, _ = fill_median(x_all[tr], x_all[te])
        ztr, zte = standardize(xtr, xte)
        model = Ridge(alpha=chosen, fit_intercept=True)
        model.fit(ztr, y_all[tr])
        oof[np.flatnonzero(te)] = model.predict(zte)
    return chosen, oof, maes


def lgbm_params():
    cfg = dict(PLAN["reliability"]["lightgbm"])
    cfg["verbosity"] = -1
    return cfg


def lgbm_oof(frame):
    cancers = list(dict.fromkeys(frame["query_cancer"]))
    x_all = frame[FEATURE_NAMES].to_numpy(dtype=np.float64)
    y_all = frame["actual"].to_numpy(dtype=np.float64)
    oof = np.full(len(frame), np.nan)
    for cancer in cancers:
        mask = frame["query_cancer"].to_numpy() == cancer
        tr = ~mask & np.isfinite(y_all)
        te = mask & np.isfinite(y_all)
        if tr.sum() < 1 or te.sum() < 1:
            continue
        xtr, xte, _ = fill_median(x_all[tr], x_all[te])
        model = LGBMRegressor(**lgbm_params())
        model.fit(xtr, y_all[tr])
        oof[np.flatnonzero(te)] = model.predict(xte)
    return oof


def fit_predict_apply(train, apply):
    train = train[np.isfinite(train["actual"])].copy()
    alpha, oof_ridge, _ = ridge_alpha_and_oof(train)
    oof_lgb = lgbm_oof(train)
    # Attach OOF predictions onto the training frame for conformal scores.
    train = train.copy()
    train["oof_ridge"] = oof_ridge
    train["oof_lgbm"] = oof_lgb
    train["oof_f0"] = train["rho_mean"]
    x_tr = train[FEATURE_NAMES].to_numpy(dtype=np.float64)
    y_tr = train["actual"].to_numpy(dtype=np.float64)
    x_ap = apply[FEATURE_NAMES].to_numpy(dtype=np.float64)
    x_tr_f, x_ap_f, _ = fill_median(x_tr, x_ap)
    z_tr, z_ap = standardize(x_tr_f, x_ap_f)
    ridge = Ridge(alpha=alpha, fit_intercept=True)
    ridge.fit(z_tr, y_tr)
    lgbm = LGBMRegressor(**lgbm_params())
    lgbm.fit(x_tr_f, y_tr)
    pred = apply.copy()
    pred["f0"] = pred["rho_mean"]
    pred["f1_ridge"] = ridge.predict(z_ap)
    pred["f1_lgbm"] = lgbm.predict(x_ap_f)
    pred["ridge_alpha"] = alpha
    widths = {}
    for name, col in (("f0", "oof_f0"), ("ridge", "oof_ridge"), ("lgbm", "oof_lgbm")):
        resid = np.abs(train[col] - train["actual"])
        widths[(name, 80)] = conformal_half_width(resid, 80)
        widths[(name, 90)] = conformal_half_width(resid, 90)
    return pred, train, widths


def metric_block(actual, pred, q80, q90):
    m = np.isfinite(actual) & np.isfinite(pred)
    a = actual[m]
    p = pred[m]
    rec = {
        "n": int(m.sum()),
        "spearman": spearman(p, a) if m.sum() else np.nan,
        "pearson": pearson(p, a) if m.sum() else np.nan,
        "mae": float(np.mean(np.abs(p - a))) if m.sum() else np.nan,
        "coverage_80": float(np.mean(np.abs(p - a) <= q80)) if m.sum() and np.isfinite(q80) else (1.0 if m.sum() and not np.isfinite(q80) else np.nan),
        "width_80": float(2 * q80) if np.isfinite(q80) else np.inf,
        "coverage_90": float(np.mean(np.abs(p - a) <= q90)) if m.sum() and np.isfinite(q90) else (1.0 if m.sum() and not np.isfinite(q90) else np.nan),
        "width_90": float(2 * q90) if np.isfinite(q90) else np.inf,
        "auroc": auroc(a, p) if m.sum() else np.nan,
    }
    return rec


def pearson(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if int(m.sum()) < INNER_N:
        return np.nan
    x = np.asarray(a, dtype=np.float64)[m]
    y = np.asarray(b, dtype=np.float64)[m]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def auroc(actual, pred):
    label = np.asarray(actual) >= 0.5
    if label.size < 1 or np.unique(label).size < 2:
        return np.nan
    return float(roc_auc_score(label, pred))


def run_reliability(p_genes, adc, stats):
    rho, _ = defined_rho()
    # Training universe stays P. ADC-only genes keep their rho for evaluation.
    p_rho = {k: v for k, v in rho.items() if k[1] in set(p_genes)}
    rate = load_rates()
    rna_mean = {(r.cancer, r.gene_symbol): float(r.rna_mean) for r in stats.itertuples(index=False)}
    pct = sd_percentile(stats)
    feature_rows = []
    pred_rows = []
    for hold in CANCERS:
        t0 = time.perf_counter()
        source = [c for c in CANCERS if c != hold]
        pieces = []
        for inner in source:
            base = [c for c in source if c != inner]
            part = design(p_genes, base, inner, p_rho, rna_mean, pct, rate, p_genes)
            part = part[np.isfinite(part["actual"])]
            pieces.append(part)
        train = pd.concat(pieces, ignore_index=True)
        apply_genes = list(dict.fromkeys(p_genes + [g for g in adc if g not in set(p_genes)]))
        apply = design(apply_genes, source, hold, rho, rna_mean, pct, rate, p_genes)
        pred, train_oof, widths = fit_predict_apply(train, apply)
        pred["heldout_cancer"] = hold
        pred["in_P"] = pred["gene_symbol"].isin(set(p_genes))
        for name, key in (("f0", "f0"), ("ridge", "ridge"), ("lgbm", "lgbm")):
            pred[f"q80_{name}"] = widths[(key, 80)]
            pred[f"q90_{name}"] = widths[(key, 90)]
        pred_rows.append(pred)
        train_oof["heldout_for_model"] = hold
        feature_rows.append(train_oof)
        print(f"REL {hold} {time.perf_counter() - t0:.1f}s rows={len(train)}", flush=True)
    preds = pd.concat(pred_rows, ignore_index=True)
    feats = pd.concat(feature_rows, ignore_index=True)
    preds.to_csv(RES / "reliability_predictions.tsv", sep="\t", index=False)
    feats.to_csv(RES / "reliability_training_rows.tsv", sep="\t", index=False)
    summarize_reliability(preds, adc)
    return preds


def summarize_reliability(preds, adc):
    models = (
        ("f0", "f0", "q80_f0", "q90_f0"),
        ("f1_ridge", "f1_ridge", "q80_ridge", "q90_ridge"),
        ("f1_lgbm", "f1_lgbm", "q80_lgbm", "q90_lgbm"),
    )
    by = []
    overall = []
    adc_rows = []
    base = preds[preds["in_P"] & np.isfinite(preds["actual"])]
    for hold, sub in base.groupby("heldout_cancer"):
        for name, col, q80, q90 in models:
            rec = metric_block(sub["actual"].to_numpy(), sub[col].to_numpy(), float(sub[q80].iloc[0]), float(sub[q90].iloc[0]))
            rec.update({"heldout_cancer": hold, "model": name})
            by.append(rec)
    for name, col, q80, q90 in models:
        # Width can differ by holdout. Coverage uses each row's own half-width.
        a = base["actual"].to_numpy()
        p = base[col].to_numpy()
        w80 = base[q80].to_numpy()
        w90 = base[q90].to_numpy()
        m = np.isfinite(a) & np.isfinite(p)
        overall.append(
            {
                "model": name,
                "n": int(m.sum()),
                "spearman": spearman(p[m], a[m]),
                "pearson": pearson(p[m], a[m]),
                "mae": float(np.mean(np.abs(p[m] - a[m]))),
                "coverage_80": float(np.mean(np.abs(p[m] - a[m]) <= w80[m])),
                "width_80": float(np.mean(2 * w80[m])),
                "coverage_90": float(np.mean(np.abs(p[m] - a[m]) <= w90[m])),
                "width_90": float(np.mean(2 * w90[m])),
                "auroc": auroc(a[m], p[m]),
            }
        )
    adc_set = set(adc)
    adc_base = preds[preds["gene_symbol"].isin(adc_set) & np.isfinite(preds["actual"])]
    for name, col, q80, q90 in models:
        a = adc_base["actual"].to_numpy()
        p = adc_base[col].to_numpy()
        w80 = adc_base[q80].to_numpy()
        w90 = adc_base[q90].to_numpy()
        m = np.isfinite(a) & np.isfinite(p)
        rec = {
            "model": name,
            "n": int(m.sum()),
            "spearman": spearman(p[m], a[m]) if m.sum() else np.nan,
            "pearson": pearson(p[m], a[m]) if m.sum() else np.nan,
            "mae": float(np.mean(np.abs(p[m] - a[m]))) if m.sum() else np.nan,
            "coverage_80": float(np.mean(np.abs(p[m] - a[m]) <= w80[m])) if m.sum() else np.nan,
            "width_80": float(np.mean(2 * w80[m])) if m.sum() else np.nan,
            "coverage_90": float(np.mean(np.abs(p[m] - a[m]) <= w90[m])) if m.sum() else np.nan,
            "width_90": float(np.mean(2 * w90[m])) if m.sum() else np.nan,
            "auroc": auroc(a[m], p[m]) if m.sum() else np.nan,
        }
        adc_rows.append(rec)
        for hold, sub in adc_base.groupby("heldout_cancer"):
            block = metric_block(
                sub["actual"].to_numpy(),
                sub[col].to_numpy(),
                float(sub[q80].iloc[0]) if len(sub) else np.nan,
                float(sub[q90].iloc[0]) if len(sub) else np.nan,
            )
            block.update({"heldout_cancer": hold, "model": name})
            adc_rows.append(block) if False else None
    # Per-cancer ADC metrics in their own rows.
    adc_by = []
    for hold, sub in adc_base.groupby("heldout_cancer"):
        for name, col, q80, q90 in models:
            block = metric_block(sub["actual"].to_numpy(), sub[col].to_numpy(), float(sub[q80].iloc[0]), float(sub[q90].iloc[0]))
            block.update({"heldout_cancer": hold, "model": name})
            adc_by.append(block)
    pd.DataFrame(by).to_csv(RES / "reliability_by_heldout.tsv", sep="\t", index=False)
    pd.DataFrame(overall).to_csv(RES / "reliability_overall.tsv", sep="\t", index=False)
    pd.DataFrame(adc_rows).to_csv(RES / "reliability_adc_overall.tsv", sep="\t", index=False)
    pd.DataFrame(adc_by).to_csv(RES / "reliability_adc_by_heldout.tsv", sep="\t", index=False)


def main():
    import multiprocessing as mp

    self_check()
    genes, p_genes, adc, missing, targets, gene_index = load_targets()
    (RES / "adc_missing_grna.txt").write_text("\n".join(missing) + ("\n" if missing else ""))
    data, stats = load_tumor(genes, targets)
    self_index = np.array([gene_index[g] for g in targets], dtype=np.int32)
    n_workers = os.cpu_count() or 1
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers) as pool:
        pairs = run_ar(data, targets, set(p_genes), self_index, pool, n_workers)
    summarize_ar(pairs, adc)
    run_reliability(p_genes, adc, stats)
    print("develop done", flush=True)


if __name__ == "__main__":
    main()
