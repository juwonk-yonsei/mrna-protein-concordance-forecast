"""APOLLO LUAD confirmation.

No arguments: exit 2 and open nothing.
--fake planted or --fake null: synthetic protein, no locked values.
--real: exit 2 unless config/unlock_B6.md exists. That file is created only
at the preregistered unlock step.
"""
import argparse
import os
import sys
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml
from scipy.stats import wilcoxon

from scripts.common import ROOT, sha256_file

PLAN_PATH = ROOT / "config" / "analysis_plan_B6.yaml"
UNLOCK = ROOT / "config" / "unlock_B6.md"
FROZEN = ROOT / "results" / "stage6" / "frozen"
PREP = ROOT / "results" / "stage6" / "prep"
OPENED = []


def log_open(path):
    OPENED.append(str(path))


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


def finite_pearson(a, b, min_n):
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
    x = x - x.mean()
    y = y - y.mean()
    den = np.sqrt((x * x).sum() * (y * y).sum())
    if den == 0:
        return np.nan
    return float((x * y).sum() / den)


def auroc_high(score, positive):
    score = np.asarray(score, dtype=np.float64)
    positive = np.asarray(positive, dtype=bool)
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = average_ranks(score)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def column_z(values):
    values = np.asarray(values, dtype=np.float64)
    out = np.zeros(values.shape, dtype=np.float64)
    m = np.isfinite(values)
    if int(m.sum()) == 0:
        return out
    mu = float(values[m].mean())
    sd = float(values[m].std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        return out
    out[m] = (values[m] - mu) / sd
    return out


def wilcoxon_greater(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    n_zero = int(np.sum(values == 0))
    rec = {"n": int(values.size), "n_zero": n_zero, "statistic": np.nan, "pvalue": np.nan, "error": ""}
    if values.size == 0:
        rec["error"] = "n=0"
        return rec
    try:
        result = wilcoxon(values, alternative="greater", zero_method="wilcox")
        rec["statistic"] = float(result.statistic)
        rec["pvalue"] = float(result.pvalue)
    except ValueError as exc:
        rec["error"] = f"ValueError: {exc}"
    return rec


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
        "frac_gt_0": float(np.mean(values > 0)),
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

    def ci(v):
        v = v[np.isfinite(v)]
        return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))

    return {
        "n": int(n),
        "S_F": s_f,
        "S_B": s_b,
        "diff": s_f - s_b if np.isfinite(s_f) and np.isfinite(s_b) else np.nan,
        "S_F_ci_low": ci(draws_f)[0],
        "S_F_ci_high": ci(draws_f)[1],
        "S_B_ci_low": ci(draws_b)[0],
        "S_B_ci_high": ci(draws_b)[1],
        "diff_ci_low": ci(diff)[0],
        "diff_ci_high": ci(diff)[1],
        "n_diff_le_0": int(np.sum(diff[ok] <= 0)),
        "n_nonfinite_draws": int((~ok).sum()),
        "pvalue": (1.0 + int(np.sum(diff[ok] <= 0))) / 2001.0,
    }


def holm(rows, order):
    ranked = sorted(rows, key=lambda row: ((row["pvalue"] if np.isfinite(row["pvalue"]) else np.inf), order.index(row["hypothesis"])))
    stopped = False
    out = []
    n = len(order)
    for k, row in enumerate(ranked, start=1):
        threshold = 0.05 / (n - k + 1)
        fails = (not np.isfinite(row["pvalue"])) or row["pvalue"] > threshold
        if stopped or fails:
            stopped = True
            call = "not_rejected"
        else:
            call = "rejected"
        out.append({"hypothesis": row["hypothesis"], "pvalue": row["pvalue"], "rank": k, "threshold": threshold, "call": call})
    return out


def load_predictions():
    with np.load(FROZEN / "apollo_sample_pred.npz", allow_pickle=False) as raw:
        log_open(FROZEN / "apollo_sample_pred.npz")
        sample = {key: np.array(raw[key]) for key in raw.files}
    gene = pd.read_csv(FROZEN / "apollo_gene_pred.tsv", sep="\t")
    log_open(FROZEN / "apollo_gene_pred.tsv")
    s_tab = pd.read_csv(FROZEN / "apollo_s.tsv", sep="\t")
    log_open(FROZEN / "apollo_s.tsv")
    rna_path = ROOT / "data" / "processed" / "apollo_luad" / "rna_log2tpm.parquet"
    log_open(rna_path)
    rna = pd.read_parquet(rna_path)
    return sample, gene, float(s_tab["value"].iloc[0]), rna


def header_symbols():
    from scripts.loader import structure_locked

    path = ROOT / "data" / "raw" / "apollo_luad" / "APOLLO1_Level3_Gprot_v061418.n87.csv"
    log_open(path)
    info = structure_locked("apollo_luad_protein")
    kept = []
    seen = {}
    for name in info["row_names"]:
        text = "" if name is None else str(name)
        if any(sep in text for sep in (";", "|", ",")):
            continue
        seen[text] = seen.get(text, 0) + 1
    dropped = {symbol for symbol, n in seen.items() if n > 1}
    present = set()
    for name in info["row_names"]:
        text = "" if name is None else str(name)
        if any(sep in text for sep in (";", "|", ",")) or text in dropped or text in present:
            continue
        present.add(text)
        kept.append(text)
    return kept, list(info["col_names"])


def align_vectors(sample, gene, rna, symbols, protein_columns):
    sample_map = pd.read_csv(ROOT / "config" / "apollo_sample_map.tsv", sep="\t")
    log_open(ROOT / "config" / "apollo_sample_map.tsv")
    mid_of = dict(zip(sample_map["protein_column"], sample_map["middle_id"]))
    mids = [mid_of[c] for c in protein_columns]
    saved_samples = [str(s) for s in sample["samples"]]
    if mids != saved_samples:
        raise RuntimeError("sample map order does not match frozen predictions")
    pred = gene.set_index("gene_symbol")
    log2 = rna.reindex(index=symbols, columns=mids).to_numpy(dtype=np.float64)
    rna_z = np.vstack([column_z(log2[i]) for i in range(len(symbols))])
    ar = np.full((len(symbols), len(mids)), np.nan)
    r1 = np.full_like(ar, np.nan)
    r2 = np.full_like(ar, np.nan)
    p1 = {str(g): i for i, g in enumerate(sample["p1_genes"])}
    q = {str(g): i for i, g in enumerate(sample["q_genes"])}
    for i, symbol in enumerate(symbols):
        if symbol in p1:
            j = p1[symbol]
            rna_z[i] = np.asarray(sample["rna_z"][:, j], dtype=np.float64)
            ar[i] = np.asarray(sample["ar_pred"][:, j], dtype=np.float64)
        if symbol in q:
            j = q[symbol]
            r1[i] = np.asarray(sample["r1_pred"][:, j], dtype=np.float64)
            r2[i] = np.asarray(sample["r2_pred"][:, j], dtype=np.float64)
    return {
        "symbols": np.array(symbols),
        "samples": np.array(mids),
        "log2": log2,
        "rna_z": rna_z,
        "ar": ar,
        "r1": r1,
        "r2": r2,
        "gene": pred,
    }


def score_genes(bundle, protein):
    u_p1 = set(pd.read_csv(PREP / "U_P1_genes.tsv", sep="\t")["gene_symbol"])
    u_p2 = set(pd.read_csv(PREP / "U_P2_genes.tsv", sep="\t")["gene_symbol"])
    u_q = set(pd.read_csv(PREP / "U_Q_genes.tsv", sep="\t")["gene_symbol"])
    log_open(PREP / "U_P1_genes.tsv")
    log_open(PREP / "U_P2_genes.tsv")
    log_open(PREP / "U_Q_genes.tsv")
    rows = []
    for i, symbol in enumerate(bundle["symbols"]):
        yy = protein[i]
        actual = finite_spearman(bundle["log2"][i], yy, 30)
        rho_s0 = finite_spearman(bundle["rna_z"][i], yy, 30) if symbol in u_p1 else np.nan
        rho_ar = finite_spearman(bundle["ar"][i], yy, 30) if symbol in u_p1 else np.nan
        rho_r1 = finite_spearman(bundle["r1"][i], yy, 30) if symbol in u_q else np.nan
        rho_r2 = finite_spearman(bundle["r2"][i], yy, 30) if symbol in u_q else np.nan
        rec = bundle["gene"].loc[symbol] if symbol in bundle["gene"].index else None
        rows.append(
            {
                "gene_symbol": symbol,
                "in_U_P1": int(symbol in u_p1),
                "in_U_P2": int(symbol in u_p2),
                "in_U_Q": int(symbol in u_q),
                "actual": actual,
                "rho_s0": rho_s0,
                "rho_ar": rho_ar,
                "rho_r1": rho_r1,
                "rho_r2": rho_r2,
                "d_p1": rho_ar - rho_s0 if np.isfinite(rho_ar) and np.isfinite(rho_s0) else np.nan,
                "d_q1": rho_ar - rho_r1 if np.isfinite(rho_ar) and np.isfinite(rho_r1) else np.nan,
                "d_q2": rho_ar - rho_r2 if np.isfinite(rho_ar) and np.isfinite(rho_r2) else np.nan,
                "f1s_lgbm": float(rec["f1s_lgbm"]) if rec is not None else np.nan,
                "b_sd": float(rec["b_sd"]) if rec is not None else np.nan,
                "f1s_lgbm_half_80": float(rec["f1s_lgbm_half_80"]) if rec is not None else np.nan,
                "f1s_lgbm_half_90": float(rec["f1s_lgbm_half_90"]) if rec is not None else np.nan,
            }
        )
    return pd.DataFrame(rows)


def location_row(name, values):
    tested = wilcoxon_greater(values)
    effect = mean_median_ci(values) if tested["n"] else {}
    row = {"hypothesis": name, "family": "primary" if name == "P1" else "secondary"}
    row.update(tested)
    row.update(effect)
    return row


def hypothesis_tables(scored):
    p1 = scored.loc[scored["in_U_P1"] == 1, "d_p1"].to_numpy(dtype=np.float64)
    q1 = scored.loc[scored["in_U_Q"] == 1, "d_q1"].to_numpy(dtype=np.float64)
    q2 = scored.loc[scored["in_U_Q"] == 1, "d_q2"].to_numpy(dtype=np.float64)
    p2_rows = scored[(scored["in_U_P2"] == 1) & np.isfinite(scored["actual"]) & np.isfinite(scored["f1s_lgbm"]) & np.isfinite(scored["b_sd"])]
    p2 = paired_spearman(p2_rows["actual"].to_numpy(), p2_rows["f1s_lgbm"].to_numpy(), p2_rows["b_sd"].to_numpy())
    rows = [location_row("P1", p1), location_row("Q1", q1), location_row("Q2", q2)]
    p2_row = {"hypothesis": "P2", "family": "primary", "n": p2["n"], "pvalue": p2["pvalue"]}
    p2_row.update(p2)
    primary = holm([rows[0], p2_row], ["P1", "P2"])
    secondary = holm(rows[1:], ["Q1", "Q2"])
    return pd.DataFrame(rows + [p2_row]), pd.DataFrame(primary + secondary)


def make_fake(kind, bundle):
    rng = np.random.default_rng(20261022)
    symbols = list(bundle["symbols"])
    n = len(symbols)
    p = bundle["log2"].shape[1]
    values = np.empty((n, p), dtype=np.float64)
    gene = bundle["gene"]
    p2 = [g for g in symbols if g in gene.index]
    rank = gene.loc[p2, "f1s_lgbm"].rank(method="average")
    a_of = {g: 0.2 + 0.7 * (float(rank.loc[g]) / len(p2)) for g in p2}
    p1 = {g for g, row in zip(symbols, bundle["ar"]) if np.isfinite(row).any()}
    for i, symbol in enumerate(symbols):
        if kind == "null" or symbol not in a_of:
            values[i] = rng.normal(size=p)
            continue
        a = a_of[symbol]
        scale = np.sqrt(max(0.0, 1.0 - a * a))
        noise = rng.normal(size=p)
        if symbol in p1:
            z = column_z(bundle["ar"][i])
        else:
            z = np.nan_to_num(bundle["rna_z"][i], nan=0.0)
        values[i] = a * z + scale * noise
    values[rng.random(values.shape) < 0.05] = np.nan
    return values


def check_hashes():
    table = pd.read_csv(ROOT / "config" / "frozen_hashes_B6.tsv", sep="\t")
    log_open(ROOT / "config" / "frozen_hashes_B6.tsv")
    bad = []
    for rec in table.itertuples(index=False):
        path = ROOT / rec.path
        log_open(path)
        digest = sha256_file(path)
        if digest != rec.sha256 or path.stat().st_size != int(rec.size_bytes):
            bad.append(rec.path)
    if bad:
        raise RuntimeError(f"frozen hash mismatch: {bad}")


def write_outputs(dest, scored, hypotheses, holm_tab, bundle, protein, s_value):
    dest.mkdir(parents=True, exist_ok=True)
    scored.to_csv(dest / "genes.tsv", sep="\t", index=False)
    hypotheses.to_csv(dest / "hypotheses.tsv", sep="\t", index=False)
    holm_tab.to_csv(dest / "holm.tsv", sep="\t", index=False)
    np.savez_compressed(dest / "protein.npz", genes=bundle["symbols"], samples=bundle["samples"], values=protein)
    np.savez_compressed(
        dest / "used_predictions.npz",
        genes=bundle["symbols"],
        samples=bundle["samples"],
        log2=bundle["log2"],
        rna_z=bundle["rna_z"],
        ar=bundle["ar"],
        r1=bundle["r1"],
        r2=bundle["r2"],
    )
    pd.DataFrame([{"item": "s", "value": s_value}]).to_csv(dest / "s.tsv", sep="\t", index=False)
    bundle["gene"].reset_index()[["gene_symbol", "f1s_lgbm", "b_sd"]].to_csv(dest / "gene_level.tsv", sep="\t", index=False)


def refresh_hashes(dest):
    rows = []
    for path in sorted(dest.iterdir()):
        if path.name in {"output_hashes.tsv", "opened_paths.tsv"}:
            continue
        rows.append({"file": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    pd.DataFrame(rows).to_csv(dest / "output_hashes.tsv", sep="\t", index=False)
    pd.DataFrame({"path": OPENED}).to_csv(dest / "opened_paths.tsv", sep="\t", index=False)


def descriptive(dest, scored, s_value, gene):
    p_path = ROOT / "results" / "P_genes.txt"
    model_path = ROOT / "results" / "ar_final_models.parquet"
    log_open(p_path)
    log_open(model_path)
    p_genes = set(pd.read_csv(p_path, sep="\t", header=None)[0].astype(str))
    model = pd.read_parquet(model_path, columns=["gene_symbol", "lambda"])
    lam = dict(zip(model["gene_symbol"], model["lambda"]))
    d = scored.loc[np.isfinite(scored["d_p1"]), "d_p1"].to_numpy(dtype=np.float64)
    q = np.quantile(d, [0.25, 0.5, 0.75]) if d.size else [np.nan, np.nan, np.nan]
    rows = [
        {"item": "n", "value": int(d.size)},
        {"item": "q1", "value": float(q[0])},
        {"item": "median", "value": float(q[1])},
        {"item": "q3", "value": float(q[2])},
        {"item": "frac_gt_0.05", "value": float(np.mean(d > 0.05)) if d.size else np.nan},
        {"item": "frac_lt_m0.05", "value": float(np.mean(d < -0.05)) if d.size else np.nan},
    ]
    use = scored[np.isfinite(scored["d_p1"])].copy()
    use["lam"] = use["gene_symbol"].map(lam)
    use["in_P"] = use["gene_symbol"].isin(p_genes)
    for value, sub in use.groupby("lam"):
        rows.append({"item": f"median_lambda_{value}", "value": float(sub["d_p1"].median())})
    for flag, sub in use.groupby("in_P"):
        rows.append({"item": f"median_in_P_{int(flag)}", "value": float(sub["d_p1"].median())})
    n_lambda_na = int(use["lam"].isna().sum())
    rows.append({"item": "n_lambda_na", "value": n_lambda_na})
    pd.DataFrame(rows).to_csv(dest / "descriptive_ar.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "item": "s",
                "value": s_value,
                "range_low": 0.8830509478307224,
                "range_high": 0.9788744329350166,
                "inside": int(0.8830509478307224 <= s_value <= 0.9788744329350166),
            }
        ]
    ).to_csv(dest / "descriptive_s.tsv", sep="\t", index=False)
    write_predictor_metrics(dest, scored, gene)
    write_calibration(dest, scored, gene)
    write_consistency(dest, scored)
    write_adc(dest, scored, gene)


def write_predictor_metrics(dest, scored, gene):
    predictors = ["f0", "f1_ridge", "f1_lgbm", "f1s_ridge", "f1s_lgbm", "b_sd"]
    extra = [f"{name}_half_{level}" for name in predictors if name != "b_sd" for level in (80, 90)]
    pred = gene.reset_index()[["gene_symbol", *predictors, *extra]]
    base = scored.loc[scored["in_U_P2"] == 1, ["gene_symbol", "actual"]].merge(pred, on="gene_symbol", how="left")
    rows = []
    actual = base["actual"].to_numpy(dtype=np.float64)
    n_actual = int(np.isfinite(actual).sum())
    for name in predictors:
        hat = base[name].to_numpy(dtype=np.float64)
        m = np.isfinite(actual) & np.isfinite(hat)
        n_used = int(m.sum())
        dropped = int(base.shape[0] - n_used)
        row = {
            "predictor": name,
            "n_U_P2": int(base.shape[0]),
            "n_actual_defined": n_actual,
            "n_used": n_used,
            "n_excluded": dropped,
            "spearman": finite_spearman(hat, actual, 3),
            "auroc_actual_ge_0.5": auroc_high(hat[m], actual[m] >= 0.5) if n_used else np.nan,
        }
        if name != "b_sd":
            row["pearson"] = finite_pearson(hat, actual, 3)
            err = actual[m] - hat[m]
            row["mae"] = float(np.mean(np.abs(err))) if n_used else np.nan
            row["mean_bias"] = float(np.mean(err)) if n_used else np.nan
            for level in (80, 90):
                half = base[f"{name}_half_{level}"].to_numpy(dtype=np.float64)
                use = m & np.isfinite(half)
                row[f"n_interval_{level}"] = int(use.sum())
                row[f"coverage_{level}"] = float(np.mean(np.abs(actual[use] - hat[use]) <= half[use])) if use.any() else np.nan
                row[f"mean_width_{level}"] = float(np.mean(2.0 * half[use])) if use.any() else np.nan
        rows.append(row)
    pd.DataFrame(rows).to_csv(dest / "descriptive_predictors.tsv", sep="\t", index=False)


def write_calibration(dest, scored, gene):
    keep = ["gene_symbol", "f1s_lgbm", "f1s_ridge", "f1s_lgbm_half_80", "f1s_lgbm_half_90", "f1s_ridge_half_80", "f1s_ridge_half_90"]
    frame = scored[["gene_symbol", "actual"]].merge(gene.reset_index()[keep], on="gene_symbol", how="inner")
    rng = np.random.default_rng(20261014)
    draws = []
    for pred in ("f1s_lgbm", "f1s_ridge"):
        half80 = f"{pred}_half_80"
        half90 = f"{pred}_half_90"
        ok = frame[np.isfinite(frame["actual"]) & np.isfinite(frame[pred])].copy()
        units = np.array(sorted(ok["gene_symbol"].astype(str).unique()))
        rng.shuffle(units)
        folds = np.array_split(units, 5)
        for fold_i, eval_units in enumerate(folds):
            eval_set = set(eval_units.tolist())
            train_units = np.array([u for u in units if u not in eval_set])
            eval_rows = ok[ok["gene_symbol"].astype(str).isin(eval_set)]
            for k in (10, 20, 50, 100):
                if train_units.size < k:
                    continue
                for rep in range(20):
                    chosen = rng.choice(train_units, size=k, replace=False)
                    cal = ok[ok["gene_symbol"].astype(str).isin(set(chosen.tolist()))]
                    actual = cal["actual"].to_numpy(dtype=np.float64)
                    hat = cal[pred].to_numpy(dtype=np.float64)
                    delta = float((actual - hat).mean())
                    abs_err = np.abs(actual - (hat + delta))
                    q80 = float(np.quantile(abs_err, 0.80))
                    q90 = float(np.quantile(abs_err, 0.90))
                    ev_a = eval_rows["actual"].to_numpy(dtype=np.float64)
                    ev_p = eval_rows[pred].to_numpy(dtype=np.float64)
                    h80 = eval_rows[half80].to_numpy(dtype=np.float64)
                    h90 = eval_rows[half90].to_numpy(dtype=np.float64)
                    use80 = np.isfinite(ev_a) & np.isfinite(ev_p) & np.isfinite(h80)
                    use90 = np.isfinite(ev_a) & np.isfinite(ev_p) & np.isfinite(h90)
                    sp = finite_spearman(ev_p, ev_a, 3)
                    ss = finite_spearman(ev_p + delta, ev_a, 3)
                    draws.append(
                        {
                            "predictor": pred,
                            "fold": fold_i,
                            "K": k,
                            "rep": rep,
                            "n_train_units": int(train_units.size),
                            "n_cal": int(len(cal)),
                            "n_eval": int(len(eval_rows)),
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
                            "spearman_abs_diff": abs(sp - ss) if np.isfinite(sp) and np.isfinite(ss) else np.nan,
                        }
                    )
    draw_df = pd.DataFrame(draws)
    draw_df.to_csv(dest / "descriptive_calibration_draws.tsv", sep="\t", index=False)
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
    summary = draw_df.groupby(["predictor", "K"], sort=True)[value_cols].mean().reset_index()
    counts = draw_df.groupby(["predictor", "K"], sort=True).size().rename("n_draws").reset_index()
    summary = summary.merge(counts, on=["predictor", "K"])
    summary.to_csv(dest / "descriptive_calibration.tsv", sep="\t", index=False)


def write_consistency(dest, scored):
    path = ROOT / "results" / "ar_pairs.tsv"
    log_open(path)
    pairs = pd.read_csv(path, sep="\t", usecols=["heldout_cancer", "gene_symbol", "s0", "delta_ar"])
    loco = pairs[pairs["heldout_cancer"] == "LUAD"][["gene_symbol", "s0", "delta_ar"]]
    merged = scored.merge(loco, on="gene_symbol", how="inner")
    d_m = np.isfinite(merged["d_p1"]) & np.isfinite(merged["delta_ar"])
    a_m = np.isfinite(merged["actual"]) & np.isfinite(merged["s0"])
    pd.DataFrame(
        [
            {
                "item": "d",
                "n": int(d_m.sum()),
                "spearman": finite_spearman(merged.loc[d_m, "d_p1"], merged.loc[d_m, "delta_ar"], 3),
            },
            {
                "item": "actual_rho",
                "n": int(a_m.sum()),
                "spearman": finite_spearman(merged.loc[a_m, "actual"], merged.loc[a_m, "s0"], 3),
            },
        ]
    ).to_csv(dest / "descriptive_consistency.tsv", sep="\t", index=False)


def write_adc(dest, scored, gene):
    path = ROOT / "config" / "adc_targets.tsv"
    log_open(path)
    targets = pd.read_csv(path, sep="\t")
    pred = gene.reset_index()[["gene_symbol", "f1s_lgbm", "f1s_lgbm_half_80", "f1s_lgbm_low_80", "f1s_lgbm_high_80"]]
    table = targets.merge(scored[["gene_symbol", "actual", "d_p1"]], on="gene_symbol", how="left")
    table = table.merge(pred, on="gene_symbol", how="left")
    rows = []
    for rec in table.itertuples(index=False):
        actual = float(rec.actual) if np.isfinite(rec.actual) else np.nan
        hat = float(rec.f1s_lgbm) if np.isfinite(rec.f1s_lgbm) else np.nan
        half = float(rec.f1s_lgbm_half_80) if np.isfinite(rec.f1s_lgbm_half_80) else np.nan
        low = float(rec.f1s_lgbm_low_80) if np.isfinite(rec.f1s_lgbm_low_80) else np.nan
        high = float(rec.f1s_lgbm_high_80) if np.isfinite(rec.f1s_lgbm_high_80) else np.nan
        covered = int(abs(actual - hat) <= half) if np.isfinite(actual) and np.isfinite(hat) and np.isfinite(half) else np.nan
        d_value = float(rec.d_p1) if np.isfinite(rec.d_p1) else np.nan
        rows.append(
            {
                "gene_symbol": rec.gene_symbol,
                "actual": actual,
                "rho_hat": hat,
                "low_80": low,
                "high_80": high,
                "covered_80": covered,
                "d": d_value,
            }
        )
    pd.DataFrame(rows).to_csv(dest / "descriptive_adc.tsv", sep="\t", index=False)


def run(kind, dest):
    log_open(PLAN_PATH)
    yaml.safe_load(PLAN_PATH.read_text())
    check_hashes()
    sample, gene, s_value, rna = load_predictions()
    if kind == "real":
        from scripts.loader import load_locked_matrix

        path = ROOT / "data" / "raw" / "apollo_luad" / "APOLLO1_Level3_Gprot_v061418.n87.csv"
        log_open(path)
        frame = load_locked_matrix("apollo_luad_protein", unlock=True)
        symbols, columns = header_symbols()
        counts = frame.index.astype(str).value_counts()
        keep = ~frame.index.astype(str).isin(set(counts[counts > 1].index))
        keep &= ~frame.index.astype(str).str.contains(r"[;|,]", regex=True)
        frame = frame.loc[keep]
        frame = frame[~frame.index.duplicated(keep="first")]
        protein = frame.reindex(index=symbols, columns=columns).to_numpy(dtype=np.float64)
    else:
        symbols, columns = header_symbols()
        bundle = align_vectors(sample, gene, rna, symbols, columns)
        protein = make_fake(kind, bundle)
    bundle = align_vectors(sample, gene, rna, symbols, columns)
    print(f"scoring {kind} genes={len(symbols)}", flush=True)
    scored = score_genes(bundle, protein)
    print(f"hypotheses {kind}", flush=True)
    hypotheses, holm_tab = hypothesis_tables(scored)
    write_outputs(dest, scored, hypotheses, holm_tab, bundle, protein, s_value)
    descriptive(dest, scored, s_value, bundle["gene"])
    refresh_hashes(dest)
    print(holm_tab.to_string(index=False), flush=True)
    print(hypotheses[["hypothesis", "n", "pvalue"]].to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fake", choices=["planted", "null"])
    parser.add_argument("--real", action="store_true")
    args = parser.parse_args()
    if args.fake is None and not args.real:
        print("refusing to run without --fake or --real", flush=True)
        refusal = ROOT / "results" / "stage6" / "locktest_noargs_opened_paths.tsv"
        refusal.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"path": OPENED}).to_csv(refusal, sep="\t", index=False)
        sys.exit(2)
    if args.real and not UNLOCK.exists():
        print("refusing --real without config/unlock_B6.md", flush=True)
        refusal = ROOT / "results" / "stage6" / "locktest_real_opened_paths.tsv"
        refusal.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"path": OPENED}).to_csv(refusal, sep="\t", index=False)
        sys.exit(2)
    if args.real:
        run("real", ROOT / "results" / "stage6" / "confirm")
    else:
        run(args.fake, ROOT / "results" / "stage6" / f"fake_{args.fake}")


if __name__ == "__main__":
    main()
