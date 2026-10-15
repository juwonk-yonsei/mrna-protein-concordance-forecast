"""Refit stage-1 R1 and R2 on all ten CPTAC cancers and score RPPA and ProCan.

The stage-1 script sets a 16-thread cap at import, so this file does not import it.
B1's n_threads field is not applied. Gene-level ridge fits run under an outer
joblib pool with n_jobs=-1; each fit uses threadpool_limits(limits=1).
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml
from joblib import Parallel, delayed
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

from scripts.b2_lib import PROC, RAW, RES, ROOT, zscore_cols
from scripts.loader import load_locked_matrix

LABEL = "사후"
MIN_N = 30
INNER_N = 3
K_EXTRA = 50
OUT = RES / "stage5" / "posthoc"


def finite_spearman(a, b, min_n):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < min_n:
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


def load_plan():
    plan = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B1.yaml"))
    if "n_threads" in plan:
        print("B1 n_threads field present and not applied", plan["n_threads"], flush=True)
    return plan


def gene_list(name):
    return [g for g in (RES / name).read_text().split() if g]


def load_tumor(genes, p_genes, cancers):
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "tumor"]
    data = {}
    for cancer in cancers:
        rna = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__rna.parquet")
        prot = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__protein.parquet")
        sids = [s for s in paired.loc[paired["cancer"] == cancer, "sample_id"] if s in rna.columns and s in prot.columns]
        x_raw = rna.reindex(index=genes, columns=sids).to_numpy(dtype=np.float64).T
        y_raw = prot.reindex(index=p_genes, columns=sids).to_numpy(dtype=np.float64).T
        data[cancer] = {"x": zscore_cols(x_raw), "y": zscore_cols(y_raw)}
        print(f"loaded {cancer} n={len(sids)}", flush=True)
        del rna, prot
    return data


def stack(data, cancers):
    x = np.vstack([data[c]["x"] for c in cancers]).astype(np.float64, copy=False)
    y = np.vstack([data[c]["y"] for c in cancers]).astype(np.float64, copy=False)
    labels = np.concatenate([np.full(data[c]["x"].shape[0], i, dtype=np.int32) for i, c in enumerate(cancers)])
    return x, y, labels


def select_features(x, y, self_index):
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


def fit_one(xf, ycol, labels, alphas):
    y_ok = np.isfinite(ycol)
    acc = np.zeros(len(alphas))
    cnt = np.zeros(len(alphas))
    chosen = alphas[0]
    with threadpool_limits(limits=1):
        for g in range(int(labels.max()) + 1):
            te = labels == g
            tr = ~te
            mtr = tr & y_ok
            mte = te & y_ok
            if int(mtr.sum()) < 1 or int(mte.sum()) < INNER_N:
                continue
            xtr = xf[mtr]
            ytr = ycol[mtr]
            xte = xf[mte]
            yte = ycol[mte]
            for a_i, alpha in enumerate(alphas):
                k = xtr.shape[1]
                gram = xtr.T @ xtr
                gram.flat[:: k + 1] += alpha
                coef = np.linalg.solve(gram, xtr.T @ ytr)
                rho = finite_spearman(xte @ coef, yte, INNER_N)
                if np.isfinite(rho):
                    acc[a_i] += rho
                    cnt[a_i] += 1
    with np.errstate(invalid="ignore"):
        mean = np.divide(acc, cnt, out=np.full(len(alphas), np.nan), where=cnt > 0)
    if not np.isfinite(mean).any():
        return 1.0, True
    best = np.nanmax(mean)
    for alpha, val in zip(alphas, mean):
        if np.isfinite(val) and val == best:
            chosen = float(alpha)
            break
    return chosen, False


def final_weights(xf, ycol, alpha):
    m = np.isfinite(ycol)
    xtr = xf[m]
    ytr = ycol[m]
    k = xtr.shape[1]
    gram = xtr.T @ xtr
    gram.flat[:: k + 1] += alpha
    with threadpool_limits(limits=1):
        return np.linalg.solve(gram, xtr.T @ ytr)


def r2_alpha(x, y, labels, alphas):
    means = []
    n_groups = int(labels.max()) + 1
    for alpha in alphas:
        per_group = []
        for g in range(n_groups):
            te = labels == g
            tr = ~te
            if int(tr.sum()) < 1 or int(te.sum()) < INNER_N:
                continue
            y_fill = np.nan_to_num(y[tr], nan=0.0)
            pred = dual_predict(x[tr], y_fill, x[te], alpha)
            rhos = [finite_spearman(pred[:, j], y[te, j], INNER_N) for j in range(y.shape[1])]
            rhos = [v for v in rhos if np.isfinite(v)]
            if rhos:
                per_group.append(float(np.mean(rhos)))
        means.append(float(np.mean(per_group)) if per_group else np.nan)
    means = np.asarray(means, dtype=np.float64)
    if not np.isfinite(means).any():
        return 1.0, means, True
    best = np.nanmax(means)
    for alpha, val in zip(alphas, means):
        if np.isfinite(val) and val == best:
            return float(alpha), means, False
    return float(alphas[0]), means, True


def dual_predict(x_tr, y_tr, x_te, alpha):
    x_tr = np.asarray(x_tr, dtype=np.float64)
    y_tr = np.asarray(y_tr, dtype=np.float64)
    n = x_tr.shape[0]
    gram = x_tr @ x_tr.T
    gram.flat[:: n + 1] += alpha
    coef = np.linalg.solve(gram, y_tr)
    return (np.asarray(x_te, dtype=np.float64) @ x_tr.T) @ coef


def self_check():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(30, 8))
    y = rng.normal(size=(30, 4))
    alpha = 2.5
    primal = np.linalg.solve(x.T @ x + alpha * np.eye(8), x.T @ y)
    dual = x.T @ np.linalg.solve(x @ x.T + alpha * np.eye(30), y)
    if not np.allclose(primal, dual, atol=1e-6):
        raise RuntimeError("dual ridge does not match primal")
    print("self-check ok", flush=True)


def rppa_header():
    import gzip

    with gzip.open(RAW / "tcga" / "TCGA-RPPA-pancan-clean.xena.gz", "rt", encoding="utf-8", errors="replace") as handle:
        return handle.readline().rstrip("\n").split("\t")[1:]


def load_rppa():
    frame = load_locked_matrix("tcga_rppa", unlock=True)
    header = rppa_header()
    if len(frame.columns) != len(header):
        raise RuntimeError("RPPA width does not match the header")
    frame.columns = header
    return frame


def z_matrix(raw):
    raw = np.array(raw, dtype=np.float64, copy=True)
    absent = ~np.isfinite(raw).any(axis=0)
    raw[:, absent] = 0.0
    z = np.asarray(zscore_cols(raw), dtype=np.float64)
    return np.nan_to_num(z, nan=0.0)


def summarize(rows):
    frame = pd.DataFrame(rows)
    unit_rows = []
    context_rows = []
    summary = []
    for platform, sub in frame.groupby("platform", sort=True):
        for model, col in (("R1", "delta_r1"), ("R2", "delta_r2"), ("AR", "delta_ar")):
            for context, ctx in sub.groupby("context", sort=True):
                d = ctx[col].to_numpy(dtype=np.float64)
                d = d[np.isfinite(d)]
                context_rows.append(
                    {
                        "platform": platform,
                        "context": context,
                        "model": model,
                        "n": int(d.size),
                        "median_delta": float(np.median(d)) if d.size else np.nan,
                        "frac_delta_lt_neg_0.05": float(np.mean(d < -0.05)) if d.size else np.nan,
                    }
                )
            means = []
            for _unit, u in sub.groupby("unit", sort=True):
                d = u[col].to_numpy(dtype=np.float64)
                d = d[np.isfinite(d)]
                if d.size:
                    means.append(float(d.mean()))
                    unit_rows.append({"platform": platform, "unit": _unit, "model": model, "n_contexts": int(d.size), "mean_delta": float(d.mean())})
            means = np.asarray(means, dtype=np.float64)
            pair = sub[col].to_numpy(dtype=np.float64)
            pair = pair[np.isfinite(pair)]
            try:
                tested = wilcoxon(means, alternative="greater", zero_method="wilcox")
                stat, pvalue, error = float(tested.statistic), float(tested.pvalue), ""
            except ValueError as exc:
                stat, pvalue, error = np.nan, np.nan, str(exc)
            summary.append(
                {
                    "platform": platform,
                    "model": model,
                    "n_units": int(means.size),
                    "mean_of_unit_means": float(means.mean()) if means.size else np.nan,
                    "median_of_unit_means": float(np.median(means)) if means.size else np.nan,
                    "wilcoxon_statistic": stat,
                    "wilcoxon_pvalue": pvalue,
                    "wilcoxon_error": error,
                    "n_pairs": int(pair.size),
                    "frac_delta_lt_neg_0.05": float(np.mean(pair < -0.05)) if pair.size else np.nan,
                }
            )
    write(OUT / "posthoc_23_summary.tsv", pd.DataFrame(summary))
    write(OUT / "posthoc_23_by_context.tsv", pd.DataFrame(context_rows))
    write(OUT / "posthoc_23_units.tsv", pd.DataFrame(unit_rows))


def eval_rppa(genes, p_index, feats, weights, x_tr, y_fill, alpha_r2, stage3):
    samples, _ = stage3.disease_samples()
    mapped = pd.read_csv(RES / "tcga_rna_id_map_b3.tsv", sep="\t", dtype=str)
    mat = stage3.load_tcga_matrix(mapped, samples)
    rppa = load_rppa()
    header_index = {}
    for j, sample in enumerate(rppa.columns):
        header_index.setdefault(sample, j)
    ab_index = {a: i for i, a in enumerate(rppa.index)}
    values = rppa.to_numpy(dtype=np.float64)
    del rppa
    pairs = pd.read_csv(RES / "stage4" / "confirm" / "rppa_pairs.tsv", sep="\t")
    h1 = set(pd.read_csv(RES / "stage4" / "confirm" / "h1_units.tsv", sep="\t")["xena_antibody"])
    pairs = pairs[pairs["xena_antibody"].isin(h1)].copy()
    rows = []
    for cancer, ids in samples.items():
        used = [s for s in ids if s in header_index]
        if not used:
            continue
        raw = mat.reindex(index=genes, columns=used).to_numpy(dtype=np.float64).T
        z = z_matrix(raw)
        pred_r2 = dual_predict(x_tr, y_fill, z, alpha_r2)
        sub = pairs[pairs["cancer"] == cancer]
        for rec in sub.itertuples(index=False):
            if rec.xena_antibody not in ab_index:
                continue
            prot = values[ab_index[rec.xena_antibody], [header_index[s] for s in used]]
            j = p_index.get(rec.gene_symbol)
            if j is None:
                r1 = r2 = np.nan
            else:
                r1 = finite_spearman(z[:, feats[:, j]] @ weights[:, j], prot, MIN_N)
                r2 = finite_spearman(pred_r2[:, j], prot, MIN_N)
            s0 = float(rec.rho_s0) if np.isfinite(rec.rho_s0) else np.nan
            rows.append(
                {
                    "platform": "rppa",
                    "unit": rec.xena_antibody,
                    "context": cancer,
                    "gene_symbol": rec.gene_symbol,
                    "delta_r1": r1 - s0 if np.isfinite(r1) and np.isfinite(s0) else np.nan,
                    "delta_r2": r2 - s0 if np.isfinite(r2) and np.isfinite(s0) else np.nan,
                    "delta_ar": float(rec.delta) if np.isfinite(rec.delta) else np.nan,
                }
            )
        print("rppa", cancer, len(used), flush=True)
    del mat, values
    return rows


def scan_procan_rna(genes, model_ids):
    import csv

    mindex = {m: i for i, m in enumerate(model_ids)}
    gindex = {g: i for i, g in enumerate(genes)}
    arr = np.full((len(model_ids), len(genes)), np.nan, dtype=np.float32)
    filled = np.zeros((len(model_ids), len(genes)), dtype=bool)
    path = RAW / "cmp" / "rnaseq_merged" / "rnaseq_merged_20260323.csv"
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        i_model = header.index("model_id")
        i_gene = header.index("gene_symbol")
        i_tpm = header.index("rsem_tpm")
        i_source = header.index("data_source")
        for parts in reader:
            if len(parts) <= i_source or parts[i_source] != "Sanger":
                continue
            model_id = parts[i_model]
            if model_id not in mindex or parts[i_gene] not in gindex:
                continue
            r = mindex[model_id]
            c = gindex[parts[i_gene]]
            if filled[r, c]:
                continue
            try:
                tpm = float(parts[i_tpm])
            except ValueError:
                continue
            if not np.isfinite(tpm):
                continue
            arr[r, c] = np.float32(np.log2(tpm + 1.0))
            filled[r, c] = True
    return arr


def eval_procan(genes, p_index, feats, weights, x_tr, y_fill, alpha_r2):
    ar = pd.read_parquet(RES / "frozen_ar_procan.parquet", columns=["lineage", "model_id"])
    by_lin = {lin: sub["model_id"].drop_duplicates().tolist() for lin, sub in ar.groupby("lineage", sort=True)}
    model_ids = list(dict.fromkeys(s for ids in by_lin.values() for s in ids))
    print("scanning procan rna", len(model_ids), flush=True)
    arr = scan_procan_rna(genes, model_ids)
    mapping = pd.read_csv(RAW / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t", dtype=str)
    project_of = dict(zip(mapping["SIDM"], mapping["Project_Identifier"]))
    protein = load_locked_matrix("procan_protein", unlock=True)
    sample_index = {s: i for i, s in enumerate(protein.index.astype(str))}
    protein_ids = pd.read_csv(RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    col_of = {p: i for i, p in enumerate(protein.columns.astype(str))}
    if list(protein.columns.astype(str)) != protein_ids["protein_id"].tolist():
        raise RuntimeError("ProCan protein ids do not match the saved list")
    values = protein.to_numpy(dtype=np.float64)
    del protein
    pairs = pd.read_csv(RES / "stage4" / "confirm" / "procan_pairs.tsv", sep="\t")
    h2 = set(pd.read_csv(RES / "stage4" / "confirm" / "h2_units.tsv", sep="\t")["protein_id"])
    pairs = pairs[pairs["protein_id"].isin(h2)].copy()
    mindex = {m: i for i, m in enumerate(model_ids)}
    rows = []
    for lin, ids in by_lin.items():
        used = [m for m in ids if m in mindex and project_of.get(m) in sample_index]
        raw = np.asarray(arr[[mindex[m] for m in used]], dtype=np.float64)
        z = z_matrix(raw)
        pred_r2 = dual_predict(x_tr, y_fill, z, alpha_r2)
        row_index = [sample_index[project_of[m]] for m in used]
        sub = pairs[pairs["lineage"] == lin]
        for rec in sub.itertuples(index=False):
            if rec.protein_id not in col_of:
                continue
            prot = values[row_index, col_of[rec.protein_id]]
            j = p_index.get(rec.gene_symbol)
            if j is None:
                r1 = r2 = np.nan
            else:
                r1 = finite_spearman(z[:, feats[:, j]] @ weights[:, j], prot, MIN_N)
                r2 = finite_spearman(pred_r2[:, j], prot, MIN_N)
            s0 = float(rec.rho_s0) if np.isfinite(rec.rho_s0) else np.nan
            rows.append(
                {
                    "platform": "procan",
                    "unit": rec.protein_id,
                    "context": lin,
                    "gene_symbol": rec.gene_symbol,
                    "delta_r1": r1 - s0 if np.isfinite(r1) and np.isfinite(s0) else np.nan,
                    "delta_r2": r2 - s0 if np.isfinite(r2) and np.isfinite(s0) else np.nan,
                    "delta_ar": float(rec.delta) if np.isfinite(rec.delta) else np.nan,
                }
            )
        print("procan", lin, len(used), flush=True)
    del values, arr
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("thread env", {k: os.environ.get(k, "") for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, flush=True)
    self_check()
    plan = load_plan()
    cancers = list(plan["cancers"])
    alphas = [float(a) for a in plan["r1"]["alphas"]]
    genes = gene_list("G_RNA_genes.txt")
    p_genes = gene_list("P_genes.txt")
    gindex = {g: i for i, g in enumerate(genes)}
    missing = [g for g in p_genes if g not in gindex]
    if missing:
        raise RuntimeError(f"P genes missing from G_RNA: {len(missing)}")
    self_index = np.array([gindex[g] for g in p_genes], dtype=np.int32)
    data = load_tumor(genes, p_genes, cancers)
    x, y, labels = stack(data, cancers)
    del data
    print("stacked", x.shape, y.shape, flush=True)
    feats = select_features(x, y, self_index)
    print("features selected", feats.shape, flush=True)
    p = y.shape[1]
    xf = np.empty((p, x.shape[0], feats.shape[0]), dtype=np.float64)
    for j in range(p):
        xf[j] = x[:, feats[:, j]]
    fitted = Parallel(n_jobs=-1, backend="loky")(delayed(fit_one)(xf[j], y[:, j], labels, alphas) for j in range(p))
    chosen = np.array([rec[0] for rec in fitted], dtype=np.float64)
    n_fallback = int(sum(rec[1] for rec in fitted))
    print("r1 alphas", n_fallback, flush=True)
    weights = np.vstack(Parallel(n_jobs=-1, backend="loky")(delayed(final_weights)(xf[j], y[:, j], chosen[j]) for j in range(p))).T
    del xf
    print("weights", weights.shape, flush=True)
    alpha_r2, r2_means, r2_fallback = r2_alpha(x, y, labels, alphas)
    print("r2", alpha_r2, r2_fallback, r2_means, flush=True)
    y_fill = np.nan_to_num(y, nan=0.0)
    pd.DataFrame(
        {
            "label": LABEL,
            "gene_symbol": p_genes,
            "alpha": chosen,
            "fallback": [rec[1] for rec in fitted],
        }
    ).to_csv(OUT / "posthoc_23_r1_alpha.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"label": LABEL, "item": "r2_alpha", "value": alpha_r2}, {"label": LABEL, "item": "r2_fallback", "value": int(r2_fallback)}, {"label": LABEL, "item": "r1_n_fallback", "value": n_fallback}]
    ).to_csv(OUT / "posthoc_23_r2_alpha.tsv", sep="\t", index=False)
    import importlib.util

    spec = importlib.util.spec_from_file_location("b3stage", ROOT / "scripts" / "08_stage3.py")
    stage3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage3)
    p_index = {g: i for i, g in enumerate(p_genes)}
    rows = eval_rppa(genes, p_index, feats, weights, x, y_fill, alpha_r2, stage3)
    rows.extend(eval_procan(genes, p_index, feats, weights, x, y_fill, alpha_r2))
    summarize(rows)
    print("r1r2 done", flush=True)


if __name__ == "__main__":
    main()
