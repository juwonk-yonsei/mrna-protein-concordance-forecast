"""Stage-4 confirmation.

Sets B_CPU_CAP=none before any project import so loader.py and common.py do
not export a thread cap. Does not read config/paths.yaml. Does not set
OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, NUMEXPR_*, NUMBA_NUM_THREADS,
n_jobs, n_threads, or taskset.

Refuses to run unless --freeze-bsd, --fake, or --unlock is passed.
--freeze-bsd and --fake do not read locked values or locked missingness.
--unlock is the single preregistered run.
"""
import argparse
import importlib.util
import os
import sys
import time
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml
from scipy.stats import wilcoxon

from scripts.b3_lib import percentile_map
from scripts.common import ROOT, sha256_file


def load_confirm():
    spec = importlib.util.spec_from_file_location("confirm_b3", ROOT / "scripts" / "09_confirm.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


C09 = load_confirm()
PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B4.yaml"))
RES = ROOT / "results"
BOOT_N = int(PLAN["hypotheses"]["bootstrap_draws"])
BOOT_SEED = int(PLAN["hypotheses"]["bootstrap_seed"])
ALPHA = float(PLAN["hypotheses"]["family_alpha_one_sided"])
MIN_N = int(PLAN["definitions"]["min_n"])
WITHIN_MIN = int(PLAN["descriptive"]["within_antibody"]["min_defined_contexts"])
CONTEXT_MIN = int(PLAN["descriptive"]["within_context_spearman_min_pairs"])
ADC_GENES = list(PLAN["descriptive"]["adc_rppa_genes"])
HYP_ORDER = list(PLAN["hypotheses"]["family"])


def out_dir(mode):
    key = "fake_dir" if mode == "fake" else "confirm_dir"
    path = ROOT / PLAN["outputs"][key]
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_table(path, frame):
    frame.to_csv(path, sep="\t", index=False)


def verify_hash_table(path):
    table = pd.read_csv(path, sep="\t", dtype=str)
    for rec in table.itertuples(index=False):
        file_path = ROOT / rec.path
        digest = sha256_file(file_path)
        size = str(file_path.stat().st_size)
        if digest != rec.sha256 or size != str(rec.size):
            raise RuntimeError(f"hash mismatch {rec.key} {rec.path}")
    return table


def append_bsd_hash(bsd_path):
    hash_path = ROOT / "config" / "frozen_hashes.tsv"
    table = verify_hash_table(hash_path)
    row = {
        "key": "bsd",
        "path": str(bsd_path.relative_to(ROOT)),
        "sha256": sha256_file(bsd_path),
        "size": str(bsd_path.stat().st_size),
    }
    kept = table[table["key"] != "bsd"]
    out = pd.concat([kept, pd.DataFrame([row])], ignore_index=True)
    protected = out[out["key"].isin(["rppa", "procan"])]
    if set(protected["key"]) != {"rppa", "procan"}:
        raise RuntimeError("rppa or procan row missing from frozen hashes")
    write_table(hash_path, out)
    verify_hash_table(hash_path)
    return row


def finite_mask(frame, columns):
    mask = np.ones(len(frame), dtype=bool)
    for col in columns:
        mask &= np.isfinite(frame[col].to_numpy(dtype=np.float64))
    return mask


def new_rng():
    return np.random.default_rng(BOOT_SEED)


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


def mean_median_ci(values, rng):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    empty = {
        "mean": np.nan,
        "median": np.nan,
        "mean_ci_low": np.nan,
        "mean_ci_high": np.nan,
        "median_ci_low": np.nan,
        "median_ci_high": np.nan,
    }
    if values.size == 0:
        return empty
    empty["mean"] = float(values.mean())
    empty["median"] = float(np.median(values))
    drawn = rng.choice(values, size=(BOOT_N, values.size), replace=True)
    means = drawn.mean(axis=1)
    medians = np.median(drawn, axis=1)
    empty["mean_ci_low"] = float(np.quantile(means, 0.025))
    empty["mean_ci_high"] = float(np.quantile(means, 0.975))
    empty["median_ci_low"] = float(np.quantile(medians, 0.025))
    empty["median_ci_high"] = float(np.quantile(medians, 0.975))
    return empty


def unit_means(pairs, unit_col):
    rows = []
    for unit, sub in pairs.groupby(unit_col, sort=True):
        d = sub["delta"].to_numpy(dtype=np.float64)
        d = d[np.isfinite(d)]
        if d.size == 0:
            continue
        rows.append({unit_col: unit, "gene_symbol": sub["gene_symbol"].iloc[0], "n_contexts": int(d.size), "mean_delta": float(d.mean())})
    return pd.DataFrame(rows)


def blank_hypothesis(name):
    return {
        "hypothesis": name,
        "n_units": np.nan,
        "n_pairs": np.nan,
        "n_zero": np.nan,
        "statistic": np.nan,
        "pvalue": np.nan,
        "mean": np.nan,
        "mean_ci_low": np.nan,
        "mean_ci_high": np.nan,
        "median": np.nan,
        "median_ci_low": np.nan,
        "median_ci_high": np.nan,
        "S_F": np.nan,
        "S_F_ci_low": np.nan,
        "S_F_ci_high": np.nan,
        "S_B": np.nan,
        "S_B_ci_low": np.nan,
        "S_B_ci_high": np.nan,
        "diff": np.nan,
        "diff_ci_low": np.nan,
        "diff_ci_high": np.nan,
        "n_diff_le_0": np.nan,
        "n_nonfinite_draws": np.nan,
        "error": "",
    }


def h_location(name, values):
    rec = blank_hypothesis(name)
    tested = wilcoxon_greater(values)
    effect = mean_median_ci(values, new_rng())
    rec.update(
        {
            "n_units": tested["n"],
            "n_zero": tested["n_zero"],
            "statistic": tested["statistic"],
            "pvalue": tested["pvalue"],
            "error": tested["error"],
        }
    )
    rec.update(effect)
    return rec


def paired_spearman_test(frame, unit_col):
    rec = blank_hypothesis("")
    use = frame.loc[finite_mask(frame, ["actual", "f1s_lgbm", "b_sd"])].copy()
    rec["n_pairs"] = int(len(use))
    if len(use) == 0:
        rec["error"] = "n=0"
        return rec
    genes = use[unit_col].to_numpy()
    uniq, codes = np.unique(genes, return_inverse=True)
    order = np.argsort(codes, kind="mergesort")
    codes_s = codes[order]
    pred = use["f1s_lgbm"].to_numpy(dtype=np.float64)[order]
    base = use["b_sd"].to_numpy(dtype=np.float64)[order]
    actual = use["actual"].to_numpy(dtype=np.float64)[order]
    breaks = np.flatnonzero(np.diff(codes_s)) + 1
    spans = np.split(np.arange(len(codes_s)), breaks)
    rec["n_units"] = int(len(uniq))
    rec["S_F"] = C09.finite_spearman(pred, actual, CONTEXT_MIN)[0]
    rec["S_B"] = C09.finite_spearman(base, actual, CONTEXT_MIN)[0]
    if np.isfinite(rec["S_F"]) and np.isfinite(rec["S_B"]):
        rec["diff"] = rec["S_F"] - rec["S_B"]
    rng = new_rng()
    draws_f = np.empty(BOOT_N)
    draws_b = np.empty(BOOT_N)
    n_units = len(uniq)
    for i in range(BOOT_N):
        chosen = rng.choice(n_units, size=n_units, replace=True)
        idx = np.concatenate([spans[j] for j in chosen])
        draws_f[i] = C09.finite_spearman(pred[idx], actual[idx], CONTEXT_MIN)[0]
        draws_b[i] = C09.finite_spearman(base[idx], actual[idx], CONTEXT_MIN)[0]
    diff = draws_f - draws_b
    ok = np.isfinite(diff)
    rec["n_nonfinite_draws"] = int((~ok).sum())
    rec["n_diff_le_0"] = int(np.sum(diff[ok] <= 0))
    rec["pvalue"] = (1.0 + rec["n_diff_le_0"]) / (BOOT_N + 1.0)
    def ci(v):
        v = v[np.isfinite(v)]
        if v.size == 0:
            return np.nan, np.nan
        return float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))
    rec["S_F_ci_low"], rec["S_F_ci_high"] = ci(draws_f)
    rec["S_B_ci_low"], rec["S_B_ci_high"] = ci(draws_b)
    rec["diff_ci_low"], rec["diff_ci_high"] = ci(diff)
    return rec


def holm_table(hypotheses):
    ranked = sorted(hypotheses, key=lambda row: (np.inf if not np.isfinite(row["pvalue"]) else row["pvalue"], HYP_ORDER.index(row["hypothesis"])))
    stopped = False
    rows = []
    for k, row in enumerate(ranked, start=1):
        threshold = ALPHA / (5 - k)
        exceeds = (not np.isfinite(row["pvalue"])) or row["pvalue"] > threshold
        if stopped or exceeds:
            stopped = True
            call = "not_rejected"
        else:
            call = "rejected"
        rows.append(
            {
                "hypothesis": row["hypothesis"],
                "pvalue": row["pvalue"],
                "rank": k,
                "threshold": threshold,
                "call": call,
            }
        )
    return pd.DataFrame(rows)


def rho_s0_from_vectors(raw, prot):
    z = C09.within_group_z(raw)
    z = np.where(np.isfinite(z), z, 0.0)
    if raw.size and not np.isfinite(raw).any():
        z = np.zeros(raw.shape, dtype=np.float64)
    return C09.finite_spearman(z, prot, MIN_N)


def rppa_used_and_protein(bundle):
    ab_index = {a: i for i, a in enumerate(bundle["rppa_antibodies"])}
    sample_index = {}
    for j, sample in enumerate(bundle["rppa_samples"]):
        sample_index.setdefault(sample, j)
    used = {}
    for cancer, ids in bundle["rppa_diseases"].items():
        used[cancer] = [s for s in ids if s in sample_index]
    return ab_index, sample_index, used


def attach_rppa_s0(pairs, bundle):
    ab_index, sample_index, used = rppa_used_and_protein(bundle)
    rho = []
    n_s0 = []
    for rec in pairs.itertuples(index=False):
        ids = used.get(rec.cancer, [])
        if rec.xena_antibody not in ab_index or not ids:
            rho.append(np.nan)
            n_s0.append(0)
            continue
        cols = [sample_index[s] for s in ids]
        prot = bundle["rppa"][ab_index[rec.xena_antibody], cols]
        raw = np.array([bundle["rppa_raw"].get((s, rec.gene_symbol), np.nan) for s in ids], dtype=np.float64)
        value, n = rho_s0_from_vectors(raw, prot)
        rho.append(value)
        n_s0.append(n)
    out = pairs.copy()
    out["rho_s0"] = rho
    out["n_rho_s0"] = n_s0
    out["delta"] = [a - s if np.isfinite(a) and np.isfinite(s) else np.nan for a, s in zip(out["ar"], out["rho_s0"])]
    return out


def procan_used(bundle):
    proteins = bundle["procan_proteins"]
    protein_index = {p: i for i, p in enumerate(proteins["protein_id"].tolist())}
    sample_index = {s: i for i, s in enumerate(bundle["procan_samples"])}
    mapping = bundle["procan_mapping"]
    project_of = dict(zip(mapping["SIDM"], mapping["Project_Identifier"]))
    lineage_of = dict(zip(mapping["SIDM"], mapping["Tissue_type"]))
    models = {}
    for sidm, lin in lineage_of.items():
        if lin not in bundle["lineages"]:
            continue
        project = project_of.get(sidm)
        if project not in sample_index:
            continue
        models.setdefault(lin, []).append(sidm)
    raw_sidms = {sidm for sidm, _gene in bundle["procan_raw"]}
    models = {lin: [s for s in sids if s in raw_sidms] for lin, sids in models.items()}
    return protein_index, sample_index, project_of, models


def attach_procan_s0(pairs, bundle):
    protein_index, sample_index, project_of, models = procan_used(bundle)
    rho = []
    n_s0 = []
    for rec in pairs.itertuples(index=False):
        sidms = models.get(rec.lineage, [])
        if rec.protein_id not in protein_index or not sidms:
            rho.append(np.nan)
            n_s0.append(0)
            continue
        col = protein_index[rec.protein_id]
        prot = np.array(
            [bundle["procan"][sample_index[project_of[s]], col] if project_of.get(s) in sample_index else np.nan for s in sidms],
            dtype=np.float64,
        )
        raw = np.array([bundle["procan_raw"].get((s, rec.gene_symbol), np.nan) for s in sidms], dtype=np.float64)
        value, n = rho_s0_from_vectors(raw, prot)
        rho.append(value)
        n_s0.append(n)
    out = pairs.copy()
    out["rho_s0"] = rho
    out["n_rho_s0"] = n_s0
    out["delta"] = [a - s if np.isfinite(a) and np.isfinite(s) else np.nan for a, s in zip(out["ar"], out["rho_s0"])]
    return out


def load_bsd():
    path = ROOT / PLAN["b_sd"]["output"]
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame["b_sd"] = pd.to_numeric(frame["b_sd"], errors="coerce")
    return frame


def merge_bsd(pairs, bsd, cohort, id_col, context_col):
    sub = bsd[bsd["cohort"] == cohort][[id_col, context_col, "b_sd"]]
    merged = pairs.merge(sub, on=[id_col, context_col], how="left")
    return merged


def pair_distribution(deltas):
    d = np.asarray(deltas, dtype=np.float64)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return {"n": 0, "median": np.nan, "q1": np.nan, "q3": np.nan, "frac_gt_0.05": np.nan, "frac_lt_m0.05": np.nan}
    return {
        "n": int(d.size),
        "median": float(np.median(d)),
        "q1": float(np.quantile(d, 0.25)),
        "q3": float(np.quantile(d, 0.75)),
        "frac_gt_0.05": float(np.mean(d > 0.05)),
        "frac_lt_m0.05": float(np.mean(d < -0.05)),
    }


def context_medians(pairs, context_col, cohort):
    rows = []
    for context, sub in pairs.groupby(context_col, sort=True):
        d = sub["delta"].to_numpy(dtype=np.float64)
        d = d[np.isfinite(d)]
        rows.append({"cohort": cohort, "context": context, "n": int(d.size), "median_delta": float(np.median(d)) if d.size else np.nan})
    return pd.DataFrame(rows)


def auroc(actual, pred):
    label = actual >= 0.5
    if np.unique(label).size < 2:
        return np.nan
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(label, pred))


def spearman_only(frame, pred_col):
    m = finite_mask(frame, ["actual", pred_col])
    a = frame.loc[m, "actual"].to_numpy(dtype=np.float64)
    p = frame.loc[m, pred_col].to_numpy(dtype=np.float64)
    rho, _ = C09.finite_spearman(p, a, CONTEXT_MIN)
    return {"n": int(a.size), "spearman": rho, "auroc": auroc(a, p) if a.size else np.nan}


def within_antibody(frame, pred_col):
    rhos = []
    n_eligible = 0
    for _, sub in frame.groupby("xena_antibody", sort=True):
        m = finite_mask(sub, ["actual", pred_col])
        if int(m.sum()) < WITHIN_MIN:
            continue
        n_eligible += 1
        rho, _ = C09.finite_spearman(sub.loc[m, pred_col], sub.loc[m, "actual"], WITHIN_MIN)
        if np.isfinite(rho):
            rhos.append(rho)
    arr = np.asarray(rhos, dtype=np.float64)
    return {
        "n_antibodies_ge_min": n_eligible,
        "n_defined": int(arr.size),
        "median": float(np.median(arr)) if arr.size else np.nan,
        "frac_gt_0": float(np.mean(arr > 0)) if arr.size else np.nan,
    }


def predictor_rows(pairs, cohort, with_antibody):
    context_col = "cancer" if cohort == "rppa" else "lineage"
    rows = []
    within_rows = []
    bias_rows = []
    context_frames = []
    names = list(PLAN["descriptive"]["predictors"])
    for name in names:
        pred = "b_sd" if name == "b_sd" else name
        if pred not in pairs.columns:
            continue
        block = spearman_only(pairs, pred)
        block.update({"cohort": cohort, "model": name, "item": "pooled"})
        if name != "b_sd":
            half80 = f"{name}_half_80"
            half90 = f"{name}_half_90"
            extra = C09.reliability_block(pairs, pred, half80, half90)
            block["pearson"] = extra["pearson"]
            block["mae"] = extra["mae"]
            block["coverage_80"] = extra["coverage_80"]
            block["width_80"] = extra["width_80"]
            block["coverage_90"] = extra["coverage_90"]
            block["width_90"] = extra["width_90"]
            block["n_with_intervals"] = extra["n"]
        rows.append(block)
        by, summary = C09.within_context(pairs, context_col, pred)
        summary.update({"cohort": cohort, "model": name, "item": "within_context"})
        within_rows.append(summary)
        by.insert(0, "model", name)
        by.insert(0, "cohort", cohort)
        context_frames.append(by)
        if name != "b_sd":
            bias = C09.bias_table(pairs, context_col, pred)
            bias.insert(0, "model", name)
            bias.insert(0, "cohort", cohort)
            bias_rows.append(bias)
        if with_antibody:
            antibody = within_antibody(pairs, pred)
            antibody.update({"cohort": cohort, "model": name, "item": "within_antibody"})
            within_rows.append(antibody)
    return rows, within_rows, bias_rows, context_frames


def cptac_compare(pairs, cohort):
    m = finite_mask(pairs, ["actual", "f0"])
    rho, n = C09.finite_spearman(pairs.loc[m, "f0"], pairs.loc[m, "actual"], CONTEXT_MIN)
    return {"cohort": cohort, "n": int(m.sum()), "spearman_f0_actual": rho, "n_spearman": n}


def ar_sets():
    antibodies = set(pd.read_parquet(RES / "frozen_ar_tcga.parquet", columns=["xena_antibody"])["xena_antibody"])
    proteins = set(pd.read_parquet(RES / "frozen_ar_procan.parquet", columns=["protein_id"])["protein_id"])
    return antibodies, proteins


def run_metrics(bundle, dest, mode):
    print("pairs", flush=True)
    rppa = C09.build_rppa_pairs(bundle)
    procan = C09.build_procan_pairs(bundle)
    rppa = attach_rppa_s0(rppa, bundle)
    procan = attach_procan_s0(procan, bundle)
    bsd = load_bsd()
    rppa = merge_bsd(rppa, bsd, "rppa", "xena_antibody", "cancer")
    procan = merge_bsd(procan, bsd, "procan", "protein_id", "lineage")
    antibodies, proteins = ar_sets()
    h1_pairs = rppa[rppa["xena_antibody"].isin(antibodies)].copy()
    h2_pairs = procan[procan["protein_id"].isin(proteins)].copy()
    h1_units = unit_means(h1_pairs, "xena_antibody")
    h2_units = unit_means(h2_pairs, "protein_id")
    h1 = h_location("H1", h1_units["mean_delta"] if len(h1_units) else np.array([]))
    h2 = h_location("H2", h2_units["mean_delta"] if len(h2_units) else np.array([]))
    h3 = paired_spearman_test(rppa, "xena_antibody")
    h3["hypothesis"] = "H3"
    h4 = paired_spearman_test(procan, "protein_id")
    h4["hypothesis"] = "H4"
    hypotheses = [h1, h2, h3, h4]
    holm = holm_table(hypotheses)
    h1_pairs_test = h_location("H1_pairs", h1_pairs["delta"])
    h2_pairs_test = h_location("H2_pairs", h2_pairs["delta"])
    h1_pairs_test["n_pairs"] = int(np.isfinite(h1_pairs["delta"]).sum())
    h2_pairs_test["n_pairs"] = int(np.isfinite(h2_pairs["delta"]).sum())
    rel_rows = []
    within_rows = []
    bias_rows = []
    context_frames = []
    for frame, cohort, with_antibody in ((rppa, "rppa", True), (procan, "procan", False)):
        rows, within, bias, contexts = predictor_rows(frame, cohort, with_antibody)
        rel_rows.extend(rows)
        within_rows.extend(within)
        bias_rows.extend(bias)
        context_frames.extend(contexts)
    adc_rppa = rppa[rppa["gene_symbol"].isin(ADC_GENES)].copy()
    adc_procan = procan[C09.boolish(procan["is_adc"])].copy()
    adc_rel = []
    adc_within = []
    adc_bias = []
    adc_contexts = []
    for frame, cohort, with_antibody in ((adc_rppa, "rppa", True), (adc_procan, "procan", False)):
        rows, within, bias, contexts = predictor_rows(frame, cohort, with_antibody)
        adc_rel.extend(rows)
        adc_within.extend(within)
        adc_bias.extend(bias)
        adc_contexts.extend(contexts)
    by_context = pd.concat(
        [
            context_medians(h1_pairs, "cancer", "rppa"),
            context_medians(h2_pairs, "lineage", "procan"),
        ],
        ignore_index=True,
    )
    dist_rows = [
        {"cohort": "rppa", **pair_distribution(h1_pairs["delta"])},
        {"cohort": "procan", **pair_distribution(h2_pairs["delta"])},
    ]
    compare = pd.DataFrame([cptac_compare(rppa, "rppa"), cptac_compare(procan, "procan")])
    counts = pd.DataFrame(
        [
            {"item": "mode", "value": mode},
            {"item": "n_rppa_pairs", "value": int(len(rppa))},
            {"item": "n_procan_pairs", "value": int(len(procan))},
            {"item": "n_h1_antibodies_in_ar", "value": int(h1_pairs["xena_antibody"].nunique())},
            {"item": "n_h1_units", "value": int(len(h1_units))},
            {"item": "n_h1_pairs", "value": int(len(h1_pairs))},
            {"item": "n_h1_delta_na", "value": int((~np.isfinite(h1_pairs["delta"])).sum())},
            {"item": "n_h2_proteins_in_ar", "value": int(h2_pairs["protein_id"].nunique())},
            {"item": "n_h2_units", "value": int(len(h2_units))},
            {"item": "n_h2_pairs", "value": int(len(h2_pairs))},
            {"item": "n_h2_delta_na", "value": int((~np.isfinite(h2_pairs["delta"])).sum())},
            {"item": "n_h3_dropped", "value": int(len(rppa) - h3["n_pairs"])},
            {"item": "n_h4_dropped", "value": int(len(procan) - h4["n_pairs"])},
            {"item": "n_adc_rppa_pairs", "value": int(len(adc_rppa))},
            {"item": "n_adc_procan_pairs", "value": int(len(adc_procan))},
            {"item": "n_adc_procan_genes", "value": int(adc_procan["gene_symbol"].nunique()) if len(adc_procan) else 0},
        ]
    )
    write_table(dest / "hypotheses.tsv", pd.DataFrame(hypotheses))
    write_table(dest / "holm.tsv", holm)
    write_table(dest / "h1_units.tsv", h1_units)
    write_table(dest / "h2_units.tsv", h2_units)
    write_table(dest / "rppa_pairs.tsv", rppa)
    write_table(dest / "procan_pairs.tsv", procan)
    write_table(dest / "ar_by_context.tsv", by_context)
    write_table(dest / "ar_pair_distribution.tsv", pd.DataFrame(dist_rows))
    write_table(dest / "ar_pair_tests.tsv", pd.DataFrame([h1_pairs_test, h2_pairs_test]))
    write_table(dest / "reliability.tsv", pd.DataFrame(rel_rows))
    write_table(dest / "within_summary.tsv", pd.DataFrame(within_rows))
    if bias_rows:
        write_table(dest / "bias.tsv", pd.concat(bias_rows, ignore_index=True))
    if context_frames:
        write_table(dest / "within_context.tsv", pd.concat(context_frames, ignore_index=True))
    write_table(dest / "adc_rppa_pairs.tsv", adc_rppa)
    write_table(dest / "adc_procan_pairs.tsv", adc_procan)
    write_table(dest / "adc_reliability.tsv", pd.DataFrame(adc_rel))
    write_table(dest / "adc_within.tsv", pd.DataFrame(adc_within))
    if adc_bias:
        write_table(dest / "adc_bias.tsv", pd.concat(adc_bias, ignore_index=True))
    if adc_contexts:
        write_table(dest / "adc_within_context.tsv", pd.concat(adc_contexts, ignore_index=True))
    write_table(dest / "cptac_compare.tsv", compare)
    write_table(dest / "counts.tsv", counts)
    print("metrics written", dest, flush=True)


def g_rna():
    return [g for g in (RES / "G_RNA_genes.txt").read_text().split() if g]


def sd_percentile(sd_by_gene, genes):
    values = {gene: sd_by_gene.get(gene, np.nan) for gene in genes}
    mapped = percentile_map(values)
    n_finite = int(sum(np.isfinite(v) for v in values.values()))
    return {gene: pair[1] for gene, pair in mapped.items()}, n_finite


def freeze_bsd():
    stage3 = C09.load_stage3()
    genes = g_rna()
    ar = pd.read_parquet(RES / "frozen_ar_tcga.parquet", columns=["cancer", "sample_id"])
    samples = {cancer: sub["sample_id"].drop_duplicates().tolist() for cancer, sub in ar.groupby("cancer", sort=True)}
    disease, _path = stage3.disease_samples()
    sym = 0
    for cancer in set(samples) | set(disease):
        sym += len(set(samples.get(cancer, [])) ^ set(disease.get(cancer, [])))
    n_tcga = len({s for ids in samples.values() for s in ids})
    print(f"tcga samples {n_tcga} symdiff {sym}", flush=True)
    mapped = pd.read_csv(RES / "tcga_rna_id_map_b3.tsv", sep="\t", dtype=str)
    mat = stage3.load_tcga_matrix(mapped, samples)
    rppa_forecast = pd.read_csv(RES / "frozen_forecasts_rppa.tsv", sep="\t", usecols=["xena_antibody", "gene_symbol", "cancer"])
    seen = set(pd.read_csv(ROOT / "config" / "rppa_excluded_seen.tsv", sep="\t", dtype=str)["antibody"])
    rppa_forecast = rppa_forecast[~rppa_forecast["xena_antibody"].isin(seen)].copy()
    rppa_rows = []
    for cancer, ids in samples.items():
        _mu, sd = stage3.gene_mean_sd(mat, ids)
        pct, n_finite = sd_percentile(sd, genes)
        part = rppa_forecast[rppa_forecast["cancer"] == cancer]
        for rec in part.itertuples(index=False):
            rppa_rows.append(
                {
                    "cohort": "rppa",
                    "xena_antibody": rec.xena_antibody,
                    "protein_id": "",
                    "gene_symbol": rec.gene_symbol,
                    "cancer": rec.cancer,
                    "lineage": "",
                    "b_sd": pct.get(rec.gene_symbol, np.nan),
                    "n_samples": len(ids),
                    "n_finite_genes": n_finite,
                }
            )
    print("rppa bsd rows", len(rppa_rows), flush=True)
    pro = pd.read_parquet(RES / "frozen_ar_procan.parquet", columns=["lineage", "model_id"])
    by_lin = {lin: sub["model_id"].drop_duplicates().tolist() for lin, sub in pro.groupby("lineage", sort=True)}
    model_ids = []
    for ids in by_lin.values():
        model_ids.extend(ids)
    model_ids = list(dict.fromkeys(model_ids))
    print("procan models", len(model_ids), {k: len(v) for k, v in by_lin.items()}, flush=True)
    mindex = {m: i for i, m in enumerate(model_ids)}
    gindex = {g: i for i, g in enumerate(genes)}
    arr = np.full((len(model_ids), len(genes)), np.nan, dtype=np.float32)
    filled = np.zeros((len(model_ids), len(genes)), dtype=bool)
    import csv

    path = C09.RAW / "cmp" / "rnaseq_merged" / "rnaseq_merged_20260323.csv"
    seen_models = set()
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
            if model_id not in mindex:
                continue
            gene = parts[i_gene]
            if gene not in gindex:
                continue
            r = mindex[model_id]
            c = gindex[gene]
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
            seen_models.add(model_id)
    pro_forecast = pd.read_csv(RES / "frozen_forecasts_procan.tsv", sep="\t", usecols=["protein_id", "gene_symbol", "lineage"])
    pro_rows = []
    for lin, ids in by_lin.items():
        rows = [mindex[m] for m in ids if m in mindex]
        block = np.asarray(arr[rows], dtype=np.float64)
        with np.errstate(all="ignore"):
            sd = np.nanstd(block, axis=0, ddof=0)
        sd_map = {genes[i]: float(sd[i]) if np.isfinite(sd[i]) else np.nan for i in range(len(genes))}
        pct, n_finite = sd_percentile(sd_map, genes)
        part = pro_forecast[pro_forecast["lineage"] == lin]
        for rec in part.itertuples(index=False):
            pro_rows.append(
                {
                    "cohort": "procan",
                    "xena_antibody": "",
                    "protein_id": rec.protein_id,
                    "gene_symbol": rec.gene_symbol,
                    "cancer": "",
                    "lineage": rec.lineage,
                    "b_sd": pct.get(rec.gene_symbol, np.nan),
                    "n_samples": len(ids),
                    "n_finite_genes": n_finite,
                }
            )
    print("procan bsd rows", len(pro_rows), "models_with_rna", len(seen_models), flush=True)
    frame = pd.DataFrame(rppa_rows + pro_rows)
    frame = frame.sort_values(["cohort", "cancer", "lineage", "xena_antibody", "protein_id"], kind="mergesort")
    bsd_path = ROOT / PLAN["b_sd"]["output"]
    write_table(bsd_path, frame)
    expected_tcga = int(PLAN["b_sd"]["expected_tcga_samples"])
    expected_pro = int(PLAN["b_sd"]["expected_procan_models"])
    check = pd.DataFrame(
        [
            {"item": "n_tcga_samples", "value": n_tcga},
            {"item": "n_tcga_symdiff_vs_stage3", "value": sym},
            {"item": "n_rppa_pairs", "value": len(rppa_rows)},
            {"item": "n_rppa_bsd_na", "value": int((~np.isfinite(pd.Series([r["b_sd"] for r in rppa_rows]))).sum())},
            {"item": "n_procan_models", "value": len(model_ids)},
            {"item": "n_procan_models_with_rna", "value": len(seen_models)},
            {"item": "n_procan_pairs", "value": len(pro_rows)},
            {"item": "n_procan_bsd_na", "value": int((~np.isfinite(pd.Series([r["b_sd"] for r in pro_rows]))).sum())},
            {"item": "n_haematopoietic", "value": len(by_lin.get("Haematopoietic and Lymphoid", []))},
            {"item": "n_large_intestine", "value": len(by_lin.get("Large Intestine", []))},
            {"item": "n_lung", "value": len(by_lin.get("Lung", []))},
        ]
    )
    write_table(RES / "frozen_baseline_bsd_check.tsv", check)
    if n_tcga != expected_tcga or sym != 0 or len(model_ids) != expected_pro:
        raise RuntimeError(f"sample set differs tcga {n_tcga} sym {sym} procan {len(model_ids)}")
    if len(rppa_rows) != int(PLAN["definitions"]["h3_n_pairs_expected"]) or len(pro_rows) != int(PLAN["definitions"]["h4_n_pairs_expected"]):
        raise RuntimeError("pair count differs from the preregistered counts")
    row = append_bsd_hash(bsd_path)
    print("bsd hash", row["sha256"], row["size"], flush=True)


def load_real_bundle():
    """Locked load. Row ids must match the structure files.

    pandas renames the second copy of a duplicate RPPA sample id. Those 10 ids
    are already outside the analysis sample set. Column order is unchanged, so
    the raw header names are restored by position.
    """
    stage3 = C09.load_stage3()
    gene_map = pd.read_csv(C09.RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    samples_by, raw_rppa, _z = C09.rppa_rna_z(stage3, sorted(set(gene_map["gene_symbol"])))
    antibodies = pd.read_csv(C09.RES / "rppa_antibodies.tsv", sep="\t", dtype=str)["antibody"].tolist()
    header_samples = C09.rppa_header_samples()
    rppa_frame = C09.load_locked_matrix("tcga_rppa", unlock=True)
    if list(rppa_frame.index) != antibodies or len(rppa_frame.columns) != len(header_samples):
        raise RuntimeError("locked RPPA ids do not match the saved structure files")
    rppa_frame.columns = header_samples
    proteins = pd.read_csv(C09.RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    samples = pd.read_csv(C09.RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"].tolist()
    procan_frame = C09.load_locked_matrix("procan_protein", unlock=True)
    if list(procan_frame.index) != samples or list(procan_frame.columns) != proteins["protein_id"].tolist():
        raise RuntimeError("locked ProCan ids do not match the saved structure files")
    lineages = C09.procan_lineages()
    symbols = [s for s in proteins["symbol"].tolist() if s not in ("", "NA")]
    mapping, raw_procan, _zpro = C09.procan_rna_log2(symbols, lineages)
    return {
        "rppa": rppa_frame.to_numpy(dtype=np.float64),
        "rppa_antibodies": antibodies,
        "rppa_samples": header_samples,
        "rppa_raw": raw_rppa,
        "rppa_diseases": samples_by,
        "procan": procan_frame.to_numpy(dtype=np.float64),
        "procan_samples": samples,
        "procan_proteins": proteins,
        "procan_raw": raw_procan,
        "procan_mapping": mapping,
        "lineages": lineages,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze-bsd", action="store_true")
    parser.add_argument("--fake", action="store_true")
    parser.add_argument("--unlock", action="store_true")
    parser.add_argument("--fix1", action="store_true")
    args = parser.parse_args()
    chosen = [args.freeze_bsd, args.fake, args.unlock, args.fix1]
    if sum(chosen) != 1:
        raise SystemExit("pass exactly one of --freeze-bsd, --fake, --unlock, --fix1")
    t0 = time.perf_counter()
    if args.freeze_bsd:
        freeze_bsd()
        mode = "freeze-bsd"
        dest = None
    else:
        verify_hash_table(ROOT / "config" / "frozen_hashes.tsv")
        if args.fake:
            bundle = C09.make_fake()
            mode = "fake"
            dest = out_dir(mode)
        elif args.fix1:
            bundle = load_real_bundle()
            mode = "fix1"
            dest = ROOT / "results" / "stage4" / "confirm_fix1"
            dest.mkdir(parents=True, exist_ok=True)
        else:
            bundle = load_real_bundle()
            mode = "unlock"
            dest = out_dir(mode)
        run_metrics(bundle, dest, mode)
    elapsed = time.perf_counter() - t0
    run_path = (dest or RES) / ("run.tsv" if dest is not None else "frozen_baseline_bsd_run.tsv")
    write_table(
        run_path,
        pd.DataFrame(
            [
                {"item": "elapsed_sec", "value": elapsed},
                {"item": "n_errors", "value": 0},
                {"item": "mode", "value": mode},
            ]
        ),
    )
    print(f"done {mode} {elapsed:.1f}s", flush=True)


if __name__ == "__main__":
    main()
