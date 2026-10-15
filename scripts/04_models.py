import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import warnings

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="Degrees of freedom")

from scripts.common import load_paths, load_plan

PATHS = load_paths()
PLAN = load_plan()
RES = PATHS["results_dir"]
PROC = PATHS["processed_dir"]
SEED = int(PLAN["seed"])
RNG = np.random.default_rng(SEED)
MIN_N = int(PLAN["evaluation"]["min_measured_tumors"])
INNER_N = int(PLAN["r1"]["inner_min_n"])
DDOF = int(PLAN["standardization"]["std_ddof"])
ALPHAS = [float(a) for a in PLAN["r1"]["alphas"]]
K_EXTRA = int(PLAN["r1"]["extra_genes"])


def spearman(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < INNER_N:
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


def ridge_fit_predict(x_tr, y_tr, x_te, alpha):
    x_tr = np.asarray(x_tr, dtype=np.float64)
    y_tr = np.asarray(y_tr, dtype=np.float64)
    k = x_tr.shape[1]
    gram = x_tr.T @ x_tr
    gram.flat[:: k + 1] += alpha
    coef = np.linalg.solve(gram, x_tr.T @ y_tr)
    return np.asarray(x_te, dtype=np.float64) @ coef


def dual_predict(x_tr, y_tr, x_te, alpha):
    x_tr = np.asarray(x_tr, dtype=np.float64)
    y_tr = np.asarray(y_tr, dtype=np.float64)
    n = x_tr.shape[0]
    gram = x_tr @ x_tr.T
    gram.flat[:: n + 1] += alpha
    coef = np.linalg.solve(gram, y_tr)
    return (np.asarray(x_te, dtype=np.float64) @ x_tr.T) @ coef


def self_check():
    x = RNG.normal(size=(30, 8))
    y = RNG.normal(size=(30, 4))
    alpha = 2.5
    primal = np.linalg.solve(x.T @ x + alpha * np.eye(8), x.T @ y)
    dual = x.T @ np.linalg.solve(x @ x.T + alpha * np.eye(30), y)
    if not np.allclose(primal, dual, atol=1e-6):
        raise RuntimeError("dual ridge does not match primal")
    base = np.array([1.0, 3.0, 2.0, 5.0, 4.0])
    if abs(spearman(base, base * -4.0 + 9.0) + 1) > 1e-8:
        raise RuntimeError("spearman self-check failed")
    print("self-check ok", flush=True)


def load_gene_list(name):
    text = (RES / name).read_text().strip().splitlines()
    return [g for g in text if g]


def load_tumor():
    genes = load_gene_list("G_RNA_genes.txt")
    p_genes = load_gene_list("P_genes.txt")
    gene_index = {g: i for i, g in enumerate(genes)}
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "tumor"]
    data = {}
    all_ids = []
    for cancer in PLAN["cancers"]:
        rna = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__rna.parquet")
        prot = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__protein.parquet")
        sids = paired.loc[paired["cancer"] == cancer, "sample_id"].tolist()
        sids = [s for s in sids if s in rna.columns and s in prot.columns]
        x_raw = rna.reindex(index=genes, columns=sids).to_numpy(dtype=np.float64).T
        y_raw = prot.reindex(index=p_genes, columns=sids).to_numpy(dtype=np.float64).T
        data[cancer] = {
            "sample_id": sids,
            "x": zscore_cols(x_raw),
            "y": zscore_cols(y_raw),
        }
        all_ids.extend(sids)
        print(f"loaded {cancer} n={len(sids)}", flush=True)
    if len(all_ids) != len(set(all_ids)):
        raise RuntimeError("tumor sample ids overlap across cancers")
    return genes, p_genes, gene_index, data


def stack(data, cancers, p_index):
    xs = [data[c]["x"] for c in cancers]
    ys = [data[c]["y"][:, p_index] if p_index is not None else data[c]["y"] for c in cancers]
    labels = np.concatenate([np.full(data[c]["x"].shape[0], i, dtype=np.int32) for i, c in enumerate(cancers)])
    return np.vstack(xs), np.vstack(ys), labels


def select_features(x, y, self_index):
    # x n x G, y n x P with nan. Missing y filled with 0 for Pearson only.
    y0 = np.nan_to_num(y.astype(np.float64), nan=0.0)
    xc = x.astype(np.float64)
    xc = xc - xc.mean(axis=0, keepdims=True)
    yc = y0 - y0.mean(axis=0, keepdims=True)
    sx = np.sqrt((xc * xc).sum(axis=0))
    sy = np.sqrt((yc * yc).sum(axis=0))
    denom = sx[:, None] * sy[None, :]
    corr = np.zeros((x.shape[1], y.shape[1]), dtype=np.float64)
    np.divide(xc.T @ yc, denom, out=corr, where=denom > 0)
    score = np.abs(corr)
    score = score + (1e-12 * (1.0 - np.arange(score.shape[0]) / score.shape[0]))[:, None]
    score[~np.isfinite(score)] = -np.inf
    for j, gi in enumerate(self_index):
        score[gi, j] = -np.inf
    k = min(K_EXTRA, score.shape[0] - 1)
    picked = np.argpartition(-score, kth=k - 1, axis=0)[:k, :]
    feats = np.empty((k + 1, y.shape[1]), dtype=np.int32)
    for j in range(y.shape[1]):
        cols = picked[:, j]
        cols = cols[np.argsort(-score[cols, j])]
        feats[0, j] = self_index[j]
        feats[1:, j] = cols
    return feats


def r1_alphas(x, y, labels, feats):
    """Return chosen alpha per protein and the per-alpha mean spearman.

    Feature columns are gathered once. RNA matrices have no missing values
    in this dataset, so the fit mask is the protein-measurement mask.
    """
    n_alpha = len(ALPHAS)
    p = y.shape[1]
    n = x.shape[0]
    k = feats.shape[0]
    n_groups = int(labels.max()) + 1 if len(labels) else 0
    acc = np.zeros((n_alpha, p), dtype=np.float64)
    cnt = np.zeros((n_alpha, p), dtype=np.int32)
    x = np.asarray(x)
    y = np.asarray(y)
    xf = np.empty((p, n, k), dtype=np.float32)
    for j in range(p):
        xf[j] = x[:, feats[:, j]]
    y_ok = np.isfinite(y)
    # 51-by-51 solves do not benefit from 16 BLAS threads.
    with threadpool_limits(limits=1):
        for g in range(n_groups):
            te = labels == g
            tr = ~te
            if int(tr.sum()) < 1 or int(te.sum()) < INNER_N:
                continue
            for j in range(p):
                mtr = tr & y_ok[:, j]
                mte = te & y_ok[:, j]
                if int(mtr.sum()) < 1 or int(mte.sum()) < INNER_N:
                    continue
                xtr = xf[j, mtr]
                ytr = y[mtr, j]
                xte = xf[j, mte]
                yte = y[mte, j]
                for a_i, alpha in enumerate(ALPHAS):
                    pred = ridge_fit_predict(xtr, ytr, xte, alpha)
                    rho = spearman(pred, yte)
                    if np.isfinite(rho):
                        acc[a_i, j] += rho
                        cnt[a_i, j] += 1
    mean = np.full((n_alpha, p), np.nan)
    np.divide(acc, cnt, out=mean, where=cnt > 0)
    chosen = np.empty(p, dtype=np.float64)
    fallback = 0
    for j in range(p):
        col = mean[:, j]
        if not np.isfinite(col).any():
            chosen[j] = float(PLAN["r1"]["fallback_alpha_if_no_valid_inner_fold"])
            fallback += 1
            continue
        best = np.nanmax(col)
        # smallest alpha among ties
        for a_i, alpha in enumerate(ALPHAS):
            if np.isfinite(col[a_i]) and col[a_i] == best:
                chosen[j] = alpha
                break
    return chosen, fallback


def r1_predict(x_tr, y_tr, x_te, feats, alphas):
    p = y_tr.shape[1]
    pred = np.full((x_te.shape[0], p), np.nan, dtype=np.float64)
    with threadpool_limits(limits=1):
        for j in range(p):
            cols = feats[:, j]
            mtr = np.isfinite(y_tr[:, j]) & np.isfinite(x_tr[:, cols]).all(axis=1)
            if mtr.sum() < 1:
                continue
            pred[:, j] = ridge_fit_predict(
                x_tr[mtr][:, cols], y_tr[mtr, j], x_te[:, cols], alphas[j]
            )
    return pred


def r2_choose_alpha(x, y, labels):
    n_groups = int(labels.max()) + 1
    means = []
    for alpha in ALPHAS:
        per_group = []
        for g in range(n_groups):
            te = labels == g
            tr = ~te
            if tr.sum() < 1 or te.sum() < INNER_N:
                continue
            y_fill = np.nan_to_num(y[tr].astype(np.float64), nan=0.0)
            pred = dual_predict(x[tr], y_fill, x[te], alpha)
            rhos = []
            for j in range(y.shape[1]):
                rho = spearman(pred[:, j], y[te, j])
                if np.isfinite(rho):
                    rhos.append(rho)
            if rhos:
                per_group.append(float(np.mean(rhos)))
        means.append(float(np.mean(per_group)) if per_group else np.nan)
    means = np.array(means, dtype=np.float64)
    if not np.isfinite(means).any():
        return float(PLAN["r2"]["fallback_alpha_if_no_valid_inner_fold"]), means, True
    best = np.nanmax(means)
    for alpha, val in zip(ALPHAS, means):
        if np.isfinite(val) and val == best:
            return float(alpha), means, False
    return float(ALPHAS[0]), means, True


def pair_spearman(pred, y):
    out = np.full(pred.shape[1], np.nan)
    n_obs = np.zeros(pred.shape[1], dtype=np.int32)
    for j in range(pred.shape[1]):
        m = np.isfinite(pred[:, j]) & np.isfinite(y[:, j])
        n_obs[j] = int(m.sum())
        if n_obs[j] >= MIN_N:
            out[j] = spearman(pred[m, j], y[m, j])
    return out, n_obs


def quant_frame(df, value_col, keys):
    rows = []
    if keys:
        groups = df.groupby(keys)
    else:
        groups = [((), df)]
    for key, sub in groups:
        s = sub[value_col].dropna()
        if not isinstance(key, tuple):
            key = (key,)
        rec = {k: v for k, v in zip(keys, key)}
        rec.update(
            {
                "n": int(s.size),
                "median": float(s.median()) if len(s) else np.nan,
                "q1": float(s.quantile(0.25)) if len(s) else np.nan,
                "q3": float(s.quantile(0.75)) if len(s) else np.nan,
            }
        )
        rows.append(rec)
    return pd.DataFrame(rows)


def run_cptac(genes, p_genes, gene_index, data):
    self_index = np.array([gene_index[g] for g in p_genes], dtype=np.int32)
    s0_rows = []
    loco_rows = []
    alpha_rows = []
    for hold in PLAN["cancers"]:
        t0 = time.perf_counter()
        train = [c for c in PLAN["cancers"] if c != hold]
        assert hold not in train
        train_ids = {s for c in train for s in data[c]["sample_id"]}
        test_ids = set(data[hold]["sample_id"])
        assert train_ids.isdisjoint(test_ids)
        x_tr, y_tr, labels = stack(data, train, None)
        x_te = data[hold]["x"]
        y_te = data[hold]["y"]
        # S0 on the holdout, from rna z of the same gene.
        s0, n_obs = pair_spearman(x_te[:, self_index], y_te)
        for j, gene in enumerate(p_genes):
            s0_rows.append(
                {"cancer": hold, "gene_symbol": gene, "n": int(n_obs[j]), "s0": s0[j], "in_P": True}
            )
        feats = select_features(x_tr, y_tr, self_index)
        alphas, n_fb = r1_alphas(x_tr, y_tr, labels, feats)
        pred_r1 = r1_predict(x_tr, y_tr, x_te, feats, alphas)
        r1, n_r1 = pair_spearman(pred_r1, y_te)
        a_r2, r2_means, r2_fb = r2_choose_alpha(x_tr, y_tr, labels)
        y_fill = np.nan_to_num(y_tr.astype(np.float64), nan=0.0)
        pred_r2 = dual_predict(x_tr, y_fill, x_te, a_r2)
        r2, n_r2 = pair_spearman(pred_r2, y_te)
        for j, gene in enumerate(p_genes):
            loco_rows.append(
                {
                    "heldout_cancer": hold,
                    "gene_symbol": gene,
                    "n": int(n_obs[j]),
                    "s0": s0[j],
                    "r1": r1[j],
                    "r2": r2[j],
                    "alpha_r1": float(alphas[j]),
                    "alpha_r2": a_r2,
                    "delta_r1": r1[j] - s0[j] if np.isfinite(r1[j]) and np.isfinite(s0[j]) else np.nan,
                    "delta_r2": r2[j] - s0[j] if np.isfinite(r2[j]) and np.isfinite(s0[j]) else np.nan,
                }
            )
        alpha_rows.append(
            {
                "heldout_cancer": hold,
                "r1_n_fallback": n_fb,
                "r2_alpha": a_r2,
                "r2_fallback": int(r2_fb),
                "r2_inner_means": ",".join(str(x) for x in r2_means),
            }
        )
        print(
            f"LOCO {hold} {time.perf_counter() - t0:.1f}s r2_alpha={a_r2} r1_fallback={n_fb}",
            flush=True,
        )
    s0_df = pd.DataFrame(s0_rows)
    loco = pd.DataFrame(loco_rows)
    # Evaluated pairs are those with n >= MIN_N. Spearman may still be NA if constant.
    s0_df.to_csv(RES / "s0_pairs.tsv", sep="\t", index=False)
    loco.to_csv(RES / "loco_pairs.tsv", sep="\t", index=False)
    pd.DataFrame(alpha_rows).to_csv(RES / "loco_alpha_by_heldout.tsv", sep="\t", index=False)
    evaluated = loco[loco["n"] >= MIN_N].copy()
    s0_eval = s0_df[s0_df["n"] >= MIN_N].copy()
    rows = []
    for cancer, sub in s0_eval.groupby("cancer"):
        s = sub["s0"].dropna()
        rows.append(
            {
                "cancer": cancer,
                "n": int(s.size),
                "median": float(s.median()) if len(s) else np.nan,
                "q1": float(s.quantile(0.25)) if len(s) else np.nan,
                "q3": float(s.quantile(0.75)) if len(s) else np.nan,
                "frac_gt_0.5": float((s > 0.5).mean()) if len(s) else np.nan,
            }
        )
    pd.DataFrame(rows).to_csv(RES / "s0_by_cancer.tsv", sep="\t", index=False)
    delta_rows = []
    for col in ("delta_r1", "delta_r2"):
        s = evaluated[col].dropna()
        delta_rows.append(
            {
                "metric": col,
                "n": int(s.size),
                "median": float(s.median()) if len(s) else np.nan,
                "q1": float(s.quantile(0.25)) if len(s) else np.nan,
                "q3": float(s.quantile(0.75)) if len(s) else np.nan,
                "frac_gt_0.05": float((s > 0.05).mean()) if len(s) else np.nan,
                "frac_lt_m0.05": float((s < -0.05).mean()) if len(s) else np.nan,
            }
        )
    pd.DataFrame(delta_rows).to_csv(RES / "loco_delta_summary.tsv", sep="\t", index=False)
    by_h = []
    for hold, sub in evaluated.groupby("heldout_cancer"):
        rec = {"heldout_cancer": hold, "n": int(sub.shape[0])}
        for col in ("delta_r1", "delta_r2", "s0", "r1", "r2"):
            s = sub[col].dropna()
            rec[f"{col}_median"] = float(s.median()) if len(s) else np.nan
        by_h.append(rec)
    pd.DataFrame(by_h).to_csv(RES / "loco_by_heldout.tsv", sep="\t", index=False)
    gene_mean = (
        evaluated.groupby("gene_symbol")[["delta_r1", "delta_r2", "s0"]]
        .mean(numeric_only=True)
        .reset_index()
    )
    topbot = []
    for col in ("delta_r1", "delta_r2"):
        ordered = gene_mean.dropna(subset=[col]).sort_values([col, "gene_symbol"], ascending=[False, True])
        for rank, rec in enumerate(ordered.head(20).to_dict("records"), start=1):
            topbot.append({"which": f"top_{col}", "rank": rank, **rec})
        ordered_b = gene_mean.dropna(subset=[col]).sort_values([col, "gene_symbol"], ascending=[True, True])
        for rank, rec in enumerate(ordered_b.head(20).to_dict("records"), start=1):
            topbot.append({"which": f"bottom_{col}", "rank": rank, **rec})
    pd.DataFrame(topbot).to_csv(RES / "loco_top_bottom.tsv", sep="\t", index=False)
    alpha_counts = evaluated["alpha_r1"].value_counts().rename_axis("alpha").reset_index(name="n_pairs")
    alpha_counts["model"] = "R1"
    r2_counts = evaluated.drop_duplicates("heldout_cancer")["alpha_r2"].value_counts().rename_axis("alpha").reset_index(name="n_heldout")
    r2_counts["model"] = "R2"
    alpha_counts.to_csv(RES / "alpha_r1_counts.tsv", sep="\t", index=False)
    r2_counts.to_csv(RES / "alpha_r2_counts.tsv", sep="\t", index=False)
    return evaluated, self_index


def adc_tables(genes, p_genes, gene_index, data):
    targets = pd.read_csv(PATHS["root"] / "config" / "adc_targets.tsv", sep="\t")
    pset = set(p_genes)
    rows = []
    # S0 and, where the gene is in P, copy LOCO metrics. Genes not in P are fit
    # only for S0 here; R1 for non-P ADC targets is fit below per holdout.
    loco = pd.read_csv(RES / "loco_pairs.tsv", sep="\t")
    for gene in targets["gene_symbol"]:
        if gene in pset:
            sub = loco[loco["gene_symbol"] == gene]
            for rec in sub.to_dict("records"):
                rows.append(
                    {
                        "gene_symbol": gene,
                        "heldout_cancer": rec["heldout_cancer"],
                        "in_P": True,
                        "n": rec["n"],
                        "s0": rec["s0"],
                        "r1": rec["r1"],
                        "r2": rec["r2"],
                    }
                )
            continue
        if gene not in gene_index:
            for cancer in PLAN["cancers"]:
                rows.append(
                    {
                        "gene_symbol": gene,
                        "heldout_cancer": cancer,
                        "in_P": False,
                        "n": 0,
                        "s0": np.nan,
                        "r1": np.nan,
                        "r2": np.nan,
                    }
                )
            continue
        # Not in P: S0 and R1. R2 targets are only P.
        gi = gene_index[gene]
        y_by = {}
        for cancer in PLAN["cancers"]:
            raw = load_one_protein(cancer, gene, data[cancer]["sample_id"])
            if raw is None:
                y_by[cancer] = np.full(len(data[cancer]["sample_id"]), np.nan, dtype=np.float32)
            else:
                y_by[cancer] = zscore_cols(raw.reshape(-1, 1)).ravel()
        for hold in PLAN["cancers"]:
            train = [c for c in PLAN["cancers"] if c != hold]
            x_tr = np.vstack([data[c]["x"] for c in train])
            y_tr = np.concatenate([y_by[c] for c in train]).reshape(-1, 1)
            labels = np.concatenate(
                [np.full(data[c]["x"].shape[0], i, dtype=np.int32) for i, c in enumerate(train)]
            )
            feats = select_features(x_tr, y_tr, np.array([gi], dtype=np.int32))
            alphas, _ = r1_alphas(x_tr, y_tr, labels, feats)
            pred = r1_predict(x_tr, y_tr, data[hold]["x"], feats, alphas)
            y_te = y_by[hold]
            n = int(np.isfinite(y_te).sum())
            s0 = spearman(data[hold]["x"][:, gi], y_te) if n >= MIN_N else np.nan
            r1 = spearman(pred[:, 0], y_te) if n >= MIN_N else np.nan
            rows.append(
                {
                    "gene_symbol": gene,
                    "heldout_cancer": hold,
                    "in_P": False,
                    "n": n,
                    "s0": s0,
                    "r1": r1,
                    "r2": np.nan,
                }
            )
    pd.DataFrame(rows).to_csv(RES / "adc_loco.tsv", sep="\t", index=False)


def load_one_protein(cancer, gene, sample_ids):
    path = PROC / "cptac" / f"{cancer}__tumor__protein.parquet"
    prot = pd.read_parquet(path)
    if gene not in prot.index:
        return None
    return prot.reindex(columns=sample_ids).loc[gene].to_numpy(dtype=np.float64)


def run_normal(genes, p_genes, gene_index):
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "normal"]
    rows = []
    adc_genes = pd.read_csv(PATHS["root"] / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    adc_rows = []
    self_index = [gene_index[g] for g in p_genes if g in gene_index]
    p_use = [g for g in p_genes if g in gene_index]
    for cancer in PLAN["cancers"]:
        rna_path = PROC / "cptac" / f"{cancer}__normal__rna.parquet"
        prot_path = PROC / "cptac" / f"{cancer}__normal__protein.parquet"
        if not rna_path.exists() or not prot_path.exists():
            rows.append({"cancer": cancer, "n_samples": 0, "n": 0, "s0_median": np.nan, "status": "file_absent"})
            continue
        rna = pd.read_parquet(rna_path)
        prot = pd.read_parquet(prot_path)
        sids = paired.loc[paired["cancer"] == cancer, "sample_id"].tolist()
        sids = [s for s in sids if s in rna.columns and s in prot.columns]
        if len(sids) == 0:
            rows.append({"cancer": cancer, "n_samples": 0, "n": 0, "s0_median": np.nan, "status": "no_paired"})
            continue
        x = zscore_cols(rna.reindex(index=genes, columns=sids).to_numpy(dtype=np.float64).T)
        y = zscore_cols(prot.reindex(index=p_use, columns=sids).to_numpy(dtype=np.float64).T)
        rhos, ns = pair_spearman(x[:, self_index], y)
        s = pd.Series(rhos).dropna()
        rows.append(
            {
                "cancer": cancer,
                "n_samples": len(sids),
                "n": int(s.size),
                "s0_median": float(s.median()) if len(s) else np.nan,
                "status": "ok",
            }
        )
        for gene in adc_genes:
            if gene not in prot.index or gene not in gene_index:
                adc_rows.append({"cancer": cancer, "gene_symbol": gene, "n": 0, "s0": np.nan})
                continue
            yy = zscore_cols(prot.reindex(columns=sids).loc[gene].to_numpy(dtype=np.float64).reshape(-1, 1)).ravel()
            xx = x[:, gene_index[gene]]
            n = int(np.isfinite(yy).sum())
            adc_rows.append(
                {
                    "cancer": cancer,
                    "gene_symbol": gene,
                    "n": n,
                    "s0": spearman(xx, yy) if n >= MIN_N else np.nan,
                }
            )
    pd.DataFrame(rows).to_csv(RES / "normal_s0_by_cancer.tsv", sep="\t", index=False)
    pd.DataFrame(adc_rows).to_csv(RES / "normal_adc.tsv", sep="\t", index=False)


def run_ccle(genes, p_genes, gene_index, data):
    rna = pd.read_parquet(PROC / "ccle" / "rna.parquet")
    prot = pd.read_parquet(PROC / "ccle" / "protein.parquet")
    meta = pd.read_csv(PROC / "ccle" / "models.tsv", sep="\t", dtype=str)
    model_ids = [m for m in meta["ModelID"] if m in rna.columns and m in prot.columns]
    cptac_ids = {s for c in data.values() for s in c["sample_id"]}
    assert set(model_ids).isdisjoint(cptac_ids)
    present = set(rna.index)
    n_missing = sum(1 for g in genes if g not in present)
    (RES / "ccle_missing_grna_n.txt").write_text(f"{n_missing}\n")
    meta = meta[meta["ModelID"].isin(model_ids)].copy()
    meta["OncotreeLineage"] = meta["OncotreeLineage"].fillna("NA")
    counts = meta["OncotreeLineage"].value_counts()
    keep_lin = [g for g, n in counts.items() if n >= int(PLAN["ccle"]["min_cell_lines"]) and g not in {"", "NA", "nan"}]
    # Train on all CPTAC cancers.
    self_index = np.array([gene_index[g] for g in p_genes], dtype=np.int32)
    x_tr, y_tr, labels = stack(data, PLAN["cancers"], None)
    feats = select_features(x_tr, y_tr, self_index)
    alphas, n_fb = r1_alphas(x_tr, y_tr, labels, feats)
    a_r2, _, r2_fb = r2_choose_alpha(x_tr, y_tr, labels)
    # Prefit R1 weights.
    weights = []
    for j in range(len(p_genes)):
        cols = feats[:, j]
        mtr = np.isfinite(y_tr[:, j]) & np.isfinite(x_tr[:, cols]).all(axis=1)
        xj = np.asarray(x_tr[mtr][:, cols], dtype=np.float64)
        yj = np.asarray(y_tr[mtr, j], dtype=np.float64)
        k = xj.shape[1]
        gram = xj.T @ xj
        gram.flat[:: k + 1] += alphas[j]
        weights.append(np.linalg.solve(gram, xj.T @ yj))
    y_fill = np.nan_to_num(y_tr.astype(np.float64), nan=0.0)
    lineage_rows = []
    adc_rows = []
    targets = pd.read_csv(PATHS["root"] / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    gene_pos = {g: i for i, g in enumerate(genes)}
    for lin in keep_lin:
        ids = meta.loc[meta["OncotreeLineage"] == lin, "ModelID"].tolist()
        x_raw = np.zeros((len(ids), len(genes)), dtype=np.float64)
        for j, g in enumerate(genes):
            if g in rna.index:
                x_raw[:, j] = pd.to_numeric(rna.loc[g, ids], errors="coerce").to_numpy(dtype=np.float64)
        y_raw = prot.reindex(index=p_genes, columns=ids).apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float64).T
        x = zscore_cols(x_raw)
        y = zscore_cols(y_raw)
        # Genes absent from CCLE were left as 0 in x_raw, so their z-score is 0.
        s0, n_obs = pair_spearman(x[:, self_index], y)
        pred_r1 = np.full_like(y, np.nan, dtype=np.float64)
        for j in range(len(p_genes)):
            pred_r1[:, j] = np.asarray(x[:, feats[:, j]], dtype=np.float64) @ weights[j]
        r1, _ = pair_spearman(pred_r1, y)
        pred_r2 = dual_predict(x_tr, y_fill, x, a_r2)
        r2, _ = pair_spearman(pred_r2, y)
        def med(a):
            s = pd.Series(a).dropna()
            return float(s.median()) if len(s) else np.nan
        d1 = np.asarray(r1, dtype=np.float64) - np.asarray(s0, dtype=np.float64)
        d2 = np.asarray(r2, dtype=np.float64) - np.asarray(s0, dtype=np.float64)

        def frac_gt(d):
            m = np.isfinite(d)
            if not m.any():
                return np.nan
            return float(np.mean(d[m] > 0.05))

        def qv(d, q):
            s = pd.Series(d).dropna()
            if s.empty:
                return np.nan
            return float(s.quantile(q))

        lineage_rows.append(
            {
                "lineage": lin,
                "n_cell_lines": len(ids),
                "n_s0": int(np.isfinite(s0).sum()),
                "s0_median": med(s0),
                "r1_median": med(r1),
                "r2_median": med(r2),
                "delta_r1_median": med(d1),
                "delta_r2_median": med(d2),
                "delta_r1_q1": qv(d1, 0.25),
                "delta_r1_q3": qv(d1, 0.75),
                "delta_r2_q1": qv(d2, 0.25),
                "delta_r2_q3": qv(d2, 0.75),
                "frac_delta_r1_gt_0.05": frac_gt(d1),
                "frac_delta_r2_gt_0.05": frac_gt(d2),
            }
        )
        for gene in targets:
            if gene not in prot.index:
                adc_rows.append({"lineage": lin, "gene_symbol": gene, "n": 0, "s0": np.nan, "r1": np.nan, "r2": np.nan})
                continue
            yy = zscore_cols(pd.to_numeric(prot.loc[gene, ids], errors="coerce").to_numpy(dtype=np.float64).reshape(-1, 1)).ravel()
            n = int(np.isfinite(yy).sum())
            if gene in gene_pos:
                xx = x[:, gene_pos[gene]]
                s0g = spearman(xx, yy) if n >= MIN_N else np.nan
            else:
                s0g = np.nan
            if gene in p_genes:
                j = p_genes.index(gene)
                r1g = spearman(pred_r1[:, j], yy) if n >= MIN_N else np.nan
                r2g = spearman(pred_r2[:, j], yy) if n >= MIN_N else np.nan
            else:
                r1g = np.nan
                r2g = np.nan
            adc_rows.append({"lineage": lin, "gene_symbol": gene, "n": n, "s0": s0g, "r1": r1g, "r2": r2g})
        print(f"CCLE {lin} n={len(ids)}", flush=True)
    pd.DataFrame(lineage_rows).to_csv(RES / "ccle_lineage_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(adc_rows).to_csv(RES / "ccle_adc_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"r1_n_fallback": n_fb, "r2_alpha": a_r2, "r2_fallback": int(r2_fb), "n_g_rna_missing_in_ccle": n_missing}]
    ).to_csv(RES / "ccle_model_alpha.tsv", sep="\t", index=False)


def main():
    self_check()
    genes, p_genes, gene_index, data = load_tumor()
    run_cptac(genes, p_genes, gene_index, data)
    adc_tables(genes, p_genes, gene_index, data)
    run_normal(genes, p_genes, gene_index)
    run_ccle(genes, p_genes, gene_index, data)
    print("models done", flush=True)


if __name__ == "__main__":
    main()
