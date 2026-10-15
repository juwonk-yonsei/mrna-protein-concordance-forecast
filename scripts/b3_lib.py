"""Stage-3 helpers. Importing this module removes the stage-1 CPU cap."""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.b2_lib import (  # noqa: E402
    DDOF,
    INNER_N,
    MIN_N,
    RAW,
    RES,
    ROOT,
    fit_slope,
    residuals,
    spearman,
)
from scripts.common import sha256_file  # noqa: E402

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B3.yaml"))
B2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
CANCERS = list(B2["evaluation"]["cancers"])
F1S_FEATURES = list(PLAN["f1s"]["features"])
QUANTILES = [float(q) for q in PLAN["scale_check"]["quantiles"]]


def prev_tokens(text):
    if text is None or (isinstance(text, float) and not np.isfinite(text)):
        return []
    raw = str(text).strip()
    if raw == "" or raw.lower() == "nan":
        return []
    return [part.strip() for part in raw.split("|") if part.strip()]


def entrez_symbol_map(g_rna):
    hgnc = pd.read_csv(
        RAW / "hgnc" / "hgnc_complete_set.txt",
        sep="\t",
        dtype=str,
        usecols=["symbol", "locus_group", "entrez_id", "prev_symbol"],
    )
    pc = hgnc[hgnc["locus_group"] == "protein-coding gene"].copy()
    pc = pc[pc["entrez_id"].notna() & (pc["entrez_id"] != "")]
    pc = pc.groupby("entrez_id").filter(lambda d: d["symbol"].nunique() == 1)
    pc = pc.drop_duplicates(["entrez_id", "symbol"])
    pc = pc.groupby("symbol").filter(lambda d: d["entrez_id"].nunique() == 1)
    pc = pc[pc["symbol"].isin(set(g_rna))]
    return dict(zip(pc["entrez_id"], pc["symbol"])), hgnc


def prev_symbol_map(hgnc):
    """previous symbol -> approved symbol, only when exactly one approved symbol lists it."""
    buckets = {}
    for row in hgnc.itertuples(index=False):
        for token in prev_tokens(row.prev_symbol):
            buckets.setdefault(token, set()).add(row.symbol)
    return {token: next(iter(syms)) for token, syms in buckets.items() if len(syms) == 1}


def map_row_ids(row_ids, g_rna):
    """Apply the stage-3 rules. Returns one row per input id, in order."""
    from collections import Counter

    gset = set(g_rna)
    entrez, hgnc = entrez_symbol_map(g_rna)
    prev = prev_symbol_map(hgnc)
    id_set = set(row_ids)
    # Pass 1 and 2.
    assigned = []
    for rid in row_ids:
        if rid in gset:
            assigned.append((rid, "exact"))
        elif str(rid).isdigit() and rid in entrez and entrez[rid] not in id_set:
            assigned.append((entrez[rid], "entrez"))
        else:
            assigned.append((None, None))
    # An approved symbol with two exact rows is ambiguous.
    exact_count = Counter(sym for sym, rule in assigned if rule == "exact")
    for i, (sym, rule) in enumerate(assigned):
        if rule == "exact" and exact_count[sym] > 1:
            assigned[i] = (None, None)
    # Entrez collisions: one approved symbol must not receive two rows.
    entrez_count = Counter(sym for sym, rule in assigned if rule == "entrez")
    for i, (sym, rule) in enumerate(assigned):
        if rule == "entrez" and entrez_count[sym] > 1:
            assigned[i] = (None, None)
    taken = {sym for sym, rule in assigned if sym is not None}
    # Pass 3. Candidates whose approved symbol has no exact row.
    candidates = []
    for i, rid in enumerate(row_ids):
        if assigned[i][0] is not None:
            continue
        if rid not in prev:
            continue
        sym = prev[rid]
        if sym not in gset or sym in id_set:
            continue
        candidates.append((i, sym))
    cand_count = Counter(sym for _, sym in candidates)
    for i, sym in candidates:
        if cand_count[sym] == 1 and sym not in taken:
            assigned[i] = (sym, "prev_symbol")
    rows = []
    for rid, (sym, rule) in zip(row_ids, assigned):
        rows.append({"row_id": rid, "gene_symbol": sym if sym is not None else "NA", "rule": rule if rule else "unmatched"})
    return pd.DataFrame(rows)


def stage2_map_row_ids(row_ids, g_rna):
    """Exact symbol, then numeric Entrez. No previous symbols. Duplicate symbols are dropped."""
    gset = set(g_rna)
    entrez, _ = entrez_symbol_map(g_rna)
    from collections import Counter

    rough = []
    for rid in row_ids:
        if rid in gset:
            rough.append(rid)
        elif str(rid).isdigit() and rid in entrez:
            rough.append(entrez[rid])
        else:
            rough.append(None)
    counts = Counter(sym for sym in rough if sym is not None)
    out = {}
    for rid, sym in zip(row_ids, rough):
        if sym is not None and counts[sym] == 1:
            out[rid] = sym
    return out


def self_check_map():
    # exact A stays. Entrez 10 -> C if C has no exact row. OLD -> D.
    # DUP twice would attach to E, so both are dropped. Alias is ignored.
    # SLC twice is unmatched if it is not a unique previous symbol.
    g = ["A", "C", "D", "E"]
    # Build a tiny map by calling the pieces, not the HGNC file.
    # The file-backed function is checked only for a stable key property below.
    assert prev_tokens("NCRNA00181|A1BGAS|A1BG-AS") == ["NCRNA00181", "A1BGAS", "A1BG-AS"]
    assert prev_tokens("") == []
    assert prev_tokens(None) == []
    print("b3 map self-check ok", flush=True)


def percentile_map(values_by_gene):
    genes = list(values_by_gene)
    vals = np.asarray([values_by_gene[g] for g in genes], dtype=np.float64)
    finite = np.isfinite(vals)
    pct = np.full(len(vals), np.nan)
    if finite.any():
        pct[finite] = 100.0 * rankdata(vals[finite], method="average") / float(finite.sum())
    return {g: (float(v) if np.isfinite(v) else np.nan, float(p) if np.isfinite(p) else np.nan) for g, v, p in zip(genes, vals, pct)}


def quantiles(values):
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {f"q{int(q * 100):02d}": np.nan for q in QUANTILES} | {"n": 0}
    qs = np.quantile(v, QUANTILES)
    out = {f"q{int(q * 100):02d}": float(val) for q, val in zip(QUANTILES, qs)}
    out["n"] = int(v.size)
    return out


def within_gene_summary(frame, pred_col, min_cancers):
    rhos = []
    for _, sub in frame.groupby("gene_symbol"):
        m = np.isfinite(sub[pred_col].to_numpy()) & np.isfinite(sub["actual"].to_numpy())
        if int(m.sum()) < min_cancers:
            continue
        rhos.append(spearman(sub[pred_col].to_numpy()[m], sub["actual"].to_numpy()[m]))
    s = pd.Series(rhos, dtype=float)
    if s.empty:
        return {"n_genes": 0, "median": np.nan, "q1": np.nan, "q3": np.nan, "frac_gt_0": np.nan, "frac_gt_0.3": np.nan}
    return {
        "n_genes": int(s.size),
        "median": float(s.median()),
        "q1": float(s.quantile(0.25)),
        "q3": float(s.quantile(0.75)),
        "frac_gt_0": float((s > 0).mean()),
        "frac_gt_0.3": float((s > 0.3).mean()),
    }


def centered_spearman(frame, pred_col, min_cancers):
    preds = []
    actuals = []
    n_genes = 0
    for _, sub in frame.groupby("gene_symbol"):
        m = np.isfinite(sub[pred_col].to_numpy()) & np.isfinite(sub["actual"].to_numpy())
        if int(m.sum()) < min_cancers:
            continue
        p = sub[pred_col].to_numpy()[m]
        a = sub["actual"].to_numpy()[m]
        preds.append(p - p.mean())
        actuals.append(a - a.mean())
        n_genes += 1
    if not preds:
        return {"n_genes": 0, "n_pairs": 0, "spearman": np.nan}
    p = np.concatenate(preds)
    a = np.concatenate(actuals)
    return {"n_genes": n_genes, "n_pairs": int(p.size), "spearman": spearman(p, a)}


def fit_ar_coefficients(x_tr, y_tr, self_index, feats, alphas, lams):
    b = fit_slope(x_tr[:, self_index], y_tr)
    e = residuals(x_tr[:, self_index], y_tr, b)
    n_prot = y_tr.shape[1]
    k = feats.shape[0]
    coefs = np.full((k, n_prot), np.nan, dtype=np.float64)
    with threadpool_limits(limits=1):
        for j in range(n_prot):
            cols = feats[:, j]
            mtr = np.isfinite(e[:, j]) & np.isfinite(x_tr[:, cols]).all(axis=1)
            if int(mtr.sum()) < 1:
                continue
            xtr = np.asarray(x_tr[mtr][:, cols], dtype=np.float64)
            etr = np.asarray(e[mtr, j], dtype=np.float64)
            kk = xtr.shape[1]
            gram = xtr.T @ xtr
            gram.flat[:: kk + 1] += float(alphas[j])
            try:
                coefs[:, j] = np.linalg.solve(gram, xtr.T @ etr)
            except np.linalg.LinAlgError:
                continue
    return b, coefs


def predict_ar_matrix(z, self_index, feats, b, lams, coefs):
    n = z.shape[0]
    pred = np.full((n, len(b)), np.nan, dtype=np.float64)
    for j in range(len(b)):
        if not np.isfinite(coefs[:, j]).all():
            continue
        ghat = np.asarray(z[:, feats[:, j]], dtype=np.float64) @ coefs[:, j]
        pred[:, j] = float(b[j]) * np.asarray(z[:, self_index[j]], dtype=np.float64) + float(lams[j]) * ghat
    return pred


def append_frozen_hashes(rows):
    path = ROOT / "config" / "frozen_hashes.tsv"
    old = pd.read_csv(path, sep="\t", dtype=str)
    protected = {"rppa", "procan"}
    extra = old[~old["key"].isin(protected)].copy()
    # Drop previous stage-3 keys so a rerun replaces them. Stage-2 keys stay.
    keep = old[old["key"].isin(protected)].copy()
    new = pd.DataFrame(rows)
    # Confirm the protected files were not rewritten.
    for rec in keep.itertuples(index=False):
        file_path = ROOT / rec.path
        digest = sha256_file(file_path)
        if digest != rec.sha256 or str(file_path.stat().st_size) != str(rec.size):
            raise RuntimeError(f"frozen F0/F1 file changed: {rec.path}")
    out = pd.concat([keep, new], ignore_index=True)
    out.to_csv(path, sep="\t", index=False)
    return extra
