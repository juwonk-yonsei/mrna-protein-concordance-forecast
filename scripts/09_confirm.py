"""Confirmation metrics.

Refuses to run unless --unlock or --fake is passed. This stage runs --fake only.
The fake path reads the RPPA header line and does not read locked values or
locked missingness. --unlock is implemented for the next stage and is not run now.
"""
import argparse
import gzip
import os
import sys
import time
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from scripts.b2_lib import DDOF, RAW, RES, ROOT
from scripts.b3_lib import PLAN
from scripts.loader import load_locked_matrix

MIN_N = 30
BOOT_N = int(PLAN["confirm"]["bootstrap_draws"])
BOOT_SEED = int(PLAN["confirm"]["bootstrap_seed"])
WITHIN_MIN = int(PLAN["confirm"]["within_unit_min_contexts"])
CENTER_MIN = int(PLAN["confirm"]["centered_min_contexts"])
RPPA_ADC = ["CDH3", "EGFR", "ERBB2", "ERBB3"]


def rppa_header_samples():
    """Sample ids from the header only. Value rows are not read."""
    path = RAW / "tcga" / "TCGA-RPPA-pancan-clean.xena.gz"
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if not header or header[0] == "":
        raise RuntimeError("RPPA header was empty")
    return header[1:]


def boolish(series):
    if series.dtype == object:
        return series.astype(str).str.lower().isin(["true", "1"])
    return series.astype(bool)


def finite_spearman(a, b, min_n):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < min_n:
        return np.nan, n
    x = a[m]
    y = b[m]
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, n
    rx = pd.Series(x).rank(method="average").to_numpy()
    ry = pd.Series(y).rank(method="average").to_numpy()
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if den == 0:
        return np.nan, n
    return float((rx * ry).sum() / den), n


def within_group_z(values):
    values = np.asarray(values, dtype=np.float64)
    out = np.full(values.shape, np.nan, dtype=np.float64)
    m = np.isfinite(values)
    if int(m.sum()) < 2:
        return out
    mu = float(values[m].mean())
    sd = float(values[m].std(ddof=DDOF))
    if not np.isfinite(sd) or sd == 0:
        return out
    out[m] = (values[m] - mu) / sd
    return out


def load_stage3():
    import importlib.util

    spec = importlib.util.spec_from_file_location("b3stage", ROOT / "scripts" / "08_stage3.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rppa_rna_z(stage3, genes):
    mapped = pd.read_csv(RES / "tcga_rna_id_map_b3.tsv", sep="\t", dtype=str)
    samples, _ = stage3.disease_samples()
    mat = stage3.load_tcga_matrix(mapped, samples)
    z_lookup = {}
    raw_lookup = {}
    for disease, ids in samples.items():
        for gene in genes:
            if gene not in mat.index:
                continue
            raw = pd.to_numeric(mat.loc[gene, ids], errors="coerce").to_numpy(dtype=np.float64)
            z = within_group_z(raw)
            for sample, rv, zv in zip(ids, raw, z):
                if np.isfinite(rv):
                    raw_lookup[(sample, gene)] = float(rv)
                if np.isfinite(zv):
                    z_lookup[(sample, gene)] = float(zv)
    return samples, raw_lookup, z_lookup


def procan_lineages():
    table = pd.read_csv(RES / "procan_sanger_by_lineage.tsv", sep="\t")
    keep = table.loc[table["n_ge_30"].astype(str).isin(["True", "true", "1"]), "Tissue_type"].tolist()
    return keep


def procan_rna_log2(symbols, lineages):
    import csv

    mapping = pd.read_csv(RAW / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t", dtype=str)
    sample_ids = set(pd.read_csv(RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    mapping = mapping[mapping["Project_Identifier"].isin(sample_ids)].copy()
    wanted = mapping[mapping["Tissue_type"].isin(lineages)]
    cands = wanted["SIDM"].drop_duplicates().tolist()
    mindex = {m: i for i, m in enumerate(cands)}
    gindex = {g: i for i, g in enumerate(symbols)}
    arr = np.full((len(cands), len(symbols)), np.nan, dtype=np.float32)
    filled = np.zeros((len(cands), len(symbols)), dtype=bool)
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
    sidm_lineage = dict(zip(wanted["SIDM"], wanted["Tissue_type"]))
    z_lookup = {}
    raw_lookup = {}
    for lin in lineages:
        rows = [i for i, sidm in enumerate(cands) if sidm_lineage.get(sidm) == lin]
        if not rows:
            continue
        block = np.asarray(arr[rows], dtype=np.float64)
        for j, gene in enumerate(symbols):
            raw = block[:, j]
            z = within_group_z(raw)
            for local, row_i in enumerate(rows):
                sidm = cands[row_i]
                if np.isfinite(raw[local]):
                    raw_lookup[(sidm, gene)] = float(raw[local])
                if np.isfinite(z[local]):
                    z_lookup[(sidm, gene)] = float(z[local])
    return mapping, raw_lookup, z_lookup


def make_fake():
    """One Generator(1): RPPA noise, ProCan noise, then the ProCan 20% mask."""
    stage3 = load_stage3()
    gene_map = pd.read_csv(RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    rppa_genes = sorted(set(gene_map["gene_symbol"]))
    samples_by, raw_rppa, z_rppa = rppa_rna_z(stage3, rppa_genes)
    antibodies = pd.read_csv(RES / "rppa_antibodies.tsv", sep="\t", dtype=str)["antibody"].tolist()
    header_samples = rppa_header_samples()
    gene_of = dict(zip(gene_map["xena_antibody"], gene_map["gene_symbol"]))
    proteins = pd.read_csv(RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    samples = pd.read_csv(RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"].tolist()
    symbols = [s for s in proteins["symbol"].tolist() if s not in ("", "NA")]
    lineages = procan_lineages()
    print("fake RNA scan", flush=True)
    mapping, raw_procan, z_procan = procan_rna_log2(symbols, lineages)
    rng = np.random.default_rng(int(PLAN["confirm"]["fake_seed"]))
    rppa = rng.standard_normal((len(antibodies), len(header_samples)))
    procan = rng.standard_normal((len(samples), len(proteins)))
    mask = rng.random((len(samples), len(proteins))) < 0.2
    columns_of = {}
    for j, sample in enumerate(header_samples):
        columns_of.setdefault(sample, []).append(j)
    for i, antibody in enumerate(antibodies):
        gene = gene_of.get(antibody)
        if gene is None:
            continue
        for sample, cols in columns_of.items():
            z = z_rppa.get((sample, gene))
            if z is None:
                continue
            for j in cols:
                rppa[i, j] = 0.5 * z + rppa[i, j]
    sidm_of = dict(zip(mapping["Project_Identifier"], mapping["SIDM"]))
    symbol_of = proteins["symbol"].tolist()
    for row_i, sample in enumerate(samples):
        sidm = sidm_of.get(sample)
        if sidm is None:
            continue
        for col_i, symbol in enumerate(symbol_of):
            if symbol in ("", "NA"):
                continue
            z = z_procan.get((sidm, symbol))
            if z is None:
                continue
            procan[row_i, col_i] = 0.5 * z + procan[row_i, col_i]
    procan[mask] = np.nan
    bundle = {
        "rppa": rppa,
        "rppa_antibodies": antibodies,
        "rppa_samples": header_samples,
        "rppa_raw": raw_rppa,
        "rppa_diseases": samples_by,
        "procan": procan,
        "procan_samples": samples,
        "procan_proteins": proteins,
        "procan_raw": raw_procan,
        "procan_mapping": mapping,
        "lineages": lineages,
    }
    print(
        f"fake shapes rppa {rppa.shape} procan {procan.shape} procan_missing {int(mask.sum())}",
        flush=True,
    )
    return bundle


def load_real_bundle():
    """Locked numeric load. Not called unless --unlock is passed."""
    stage3 = load_stage3()
    gene_map = pd.read_csv(RES / "rppa_total_protein_map.tsv", sep="\t", dtype=str)
    samples_by, raw_rppa, _z = rppa_rna_z(stage3, sorted(set(gene_map["gene_symbol"])))
    antibodies = pd.read_csv(RES / "rppa_antibodies.tsv", sep="\t", dtype=str)["antibody"].tolist()
    header_samples = rppa_header_samples()
    rppa_frame = load_locked_matrix("tcga_rppa", unlock=True)
    if list(rppa_frame.index) != antibodies or list(rppa_frame.columns) != header_samples:
        raise RuntimeError("locked RPPA ids do not match the saved structure files")
    proteins = pd.read_csv(RES / "procan_protein_ids.tsv", sep="\t", dtype=str, keep_default_na=False)
    samples = pd.read_csv(RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"].tolist()
    procan_frame = load_locked_matrix("procan_protein", unlock=True)
    if list(procan_frame.index) != samples or list(procan_frame.columns) != proteins["protein_id"].tolist():
        raise RuntimeError("locked ProCan ids do not match the saved structure files")
    lineages = procan_lineages()
    symbols = [s for s in proteins["symbol"].tolist() if s not in ("", "NA")]
    mapping, raw_procan, _zpro = procan_rna_log2(symbols, lineages)
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


def test_record(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    rec = {
        "n_units": int(values.size),
        "mean": float(values.mean()) if values.size else np.nan,
        "statistic": np.nan,
        "pvalue": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "error": "",
    }
    if values.size == 0:
        rec["error"] = "n=0"
        return rec
    try:
        result = wilcoxon(values, alternative="greater", zero_method="wilcox")
        rec["statistic"] = float(result.statistic)
        rec["pvalue"] = float(result.pvalue)
    except ValueError as exc:
        rec["error"] = f"ValueError: {exc}"
    rng = np.random.default_rng(BOOT_SEED)
    draws = rng.choice(values, size=(BOOT_N, values.size), replace=True).mean(axis=1)
    rec["ci_low"] = float(np.quantile(draws, 0.025))
    rec["ci_high"] = float(np.quantile(draws, 0.975))
    return rec


def context_summary(pairs, context_col):
    rows = []
    for context, sub in pairs.groupby(context_col, sort=True):
        d = sub["delta"].to_numpy(dtype=np.float64)
        d = d[np.isfinite(d)]
        rows.append(
            {
                context_col: context,
                "n": int(d.size),
                "median_delta": float(np.median(d)) if d.size else np.nan,
                "frac_lt_m0.05": float(np.mean(d < -0.05)) if d.size else np.nan,
            }
        )
    return pd.DataFrame(rows)


def unit_mean_table(pairs, unit_col):
    rows = []
    for unit, sub in pairs.groupby(unit_col, sort=True):
        d = sub["delta"].to_numpy(dtype=np.float64)
        d = d[np.isfinite(d)]
        if d.size == 0:
            continue
        rows.append({unit_col: unit, "gene_symbol": sub["gene_symbol"].iloc[0], "n_contexts": int(d.size), "mean_delta": float(d.mean())})
    return pd.DataFrame(rows)


def reliability_block(frame, pred_col, half80, half90):
    actual = frame["actual"].to_numpy(dtype=np.float64)
    pred = frame[pred_col].to_numpy(dtype=np.float64)
    h80 = frame[half80].to_numpy(dtype=np.float64)
    h90 = frame[half90].to_numpy(dtype=np.float64)
    m = np.isfinite(actual) & np.isfinite(pred) & np.isfinite(h80) & np.isfinite(h90)
    a = actual[m]
    p = pred[m]
    if a.size == 0:
        return {"n": 0, "spearman": np.nan, "pearson": np.nan, "mae": np.nan, "coverage_80": np.nan, "width_80": np.nan, "coverage_90": np.nan, "width_90": np.nan, "auroc": np.nan}
    spear, _ = finite_spearman(p, a, 3)
    pear, _ = finite_spearman(p, a, 3)
    # Pearson, not Spearman. Recompute with corrcoef.
    if np.unique(p).size >= 2 and np.unique(a).size >= 2:
        pear = float(np.corrcoef(p, a)[0, 1])
    else:
        pear = np.nan
    label = a >= 0.5
    auroc = np.nan
    if np.unique(label).size == 2:
        from sklearn.metrics import roc_auc_score

        auroc = float(roc_auc_score(label, p))
    return {
        "n": int(a.size),
        "spearman": spear,
        "pearson": pear,
        "mae": float(np.mean(np.abs(p - a))),
        "coverage_80": float(np.mean(np.abs(p - a) <= h80[m])),
        "width_80": float(np.mean(2 * h80[m])),
        "coverage_90": float(np.mean(np.abs(p - a) <= h90[m])),
        "width_90": float(np.mean(2 * h90[m])),
        "auroc": auroc,
    }


def within_context(frame, context_col, pred_col):
    rows = []
    for context, sub in frame.groupby(context_col, sort=True):
        rho, n = finite_spearman(sub[pred_col], sub["actual"], 3)
        rows.append({context_col: context, "n": n, "spearman": rho})
    out = pd.DataFrame(rows)
    defined = out["spearman"].dropna()
    summary = {
        "n_contexts": int(defined.size),
        "median_spearman": float(defined.median()) if defined.size else np.nan,
    }
    return out, summary


def within_unit(frame, unit_col, pred_col):
    rows = []
    for unit, sub in frame.groupby(unit_col, sort=True):
        rho, n = finite_spearman(sub[pred_col], sub["actual"], WITHIN_MIN)
        rows.append({unit_col: unit, "n_contexts": n, "spearman": rho})
    units = pd.DataFrame(rows)
    defined = units["spearman"].dropna()
    summary = {
        "n_units": int(defined.size),
        "median": float(defined.median()) if defined.size else np.nan,
        "q1": float(defined.quantile(0.25)) if defined.size else np.nan,
        "q3": float(defined.quantile(0.75)) if defined.size else np.nan,
    }
    if pred_col == "f0":
        summary.update({"statistic": np.nan, "pvalue": np.nan, "error": "F0 is constant across contexts"})
        return units, summary
    tested = test_record(defined.to_numpy())
    summary.update({k: tested[k] for k in ("statistic", "pvalue", "error") if k in tested or True})
    summary["statistic"] = tested["statistic"]
    summary["pvalue"] = tested["pvalue"]
    summary["error"] = tested["error"]
    return units, summary


def centered(frame, unit_col, pred_col):
    preds = []
    actuals = []
    n_units = 0
    for _, sub in frame.groupby(unit_col, sort=True):
        m = np.isfinite(sub[pred_col].to_numpy()) & np.isfinite(sub["actual"].to_numpy())
        if int(m.sum()) < CENTER_MIN:
            continue
        p = sub[pred_col].to_numpy()[m]
        a = sub["actual"].to_numpy()[m]
        preds.append(p - p.mean())
        actuals.append(a - a.mean())
        n_units += 1
    if not preds:
        return {"n_units": 0, "n_pairs": 0, "spearman": np.nan}
    p = np.concatenate(preds)
    a = np.concatenate(actuals)
    rho, _ = finite_spearman(p, a, 3)
    return {"n_units": n_units, "n_pairs": int(p.size), "spearman": rho}


def bias_table(frame, context_col, pred_col):
    rows = []
    for context, sub in frame.groupby(context_col, sort=True):
        m = np.isfinite(sub["actual"]) & np.isfinite(sub[pred_col])
        diff = sub.loc[m, "actual"].to_numpy() - sub.loc[m, pred_col].to_numpy()
        rows.append({context_col: context, "n": int(diff.size), "mean_actual_minus_pred": float(diff.mean()) if diff.size else np.nan})
    return pd.DataFrame(rows)


def build_rppa_pairs(bundle):
    seen = set(pd.read_csv(ROOT / "config" / "rppa_excluded_seen.tsv", sep="\t", dtype=str)["antibody"])
    f0 = pd.read_csv(RES / "frozen_forecasts_rppa.tsv", sep="\t")
    f1s = pd.read_csv(RES / "frozen_forecasts_rppa_F1s.tsv", sep="\t")
    ar = pd.read_parquet(RES / "frozen_ar_tcga.parquet")
    f0 = f0[~f0["xena_antibody"].isin(seen)].copy()
    f1s = f1s[~f1s["xena_antibody"].isin(seen)].copy()
    ab_index = {a: i for i, a in enumerate(bundle["rppa_antibodies"])}
    sample_index = {}
    for j, sample in enumerate(bundle["rppa_samples"]):
        sample_index.setdefault(sample, j)
    rows = []
    ar_groups = {(r.xena_antibody, r.cancer): r for r in []}
    ar = ar[~ar["xena_antibody"].isin(seen)]
    grouped = {key: sub for key, sub in ar.groupby(["xena_antibody", "cancer"], sort=False)}
    f1s_key = f1s.set_index(["xena_antibody", "cancer"])
    for rec in f0.itertuples(index=False):
        ids = bundle["rppa_diseases"].get(rec.cancer, [])
        cols = [sample_index[s] for s in ids if s in sample_index]
        if rec.xena_antibody not in ab_index or not cols:
            prot = np.array([])
            used = []
        else:
            prot = bundle["rppa"][ab_index[rec.xena_antibody], cols]
            used = [s for s in ids if s in sample_index]
        rna = np.array([bundle["rppa_raw"].get((s, rec.gene_symbol), np.nan) for s in used], dtype=np.float64)
        actual, n_s0 = finite_spearman(rna, prot, MIN_N)
        sub = grouped.get((rec.xena_antibody, rec.cancer))
        if sub is None or not used:
            ar_rho, n_ar = np.nan, 0
        else:
            pred_map = dict(zip(sub["sample_id"], sub["ar_pred"]))
            ar_vec = np.array([pred_map.get(s, np.nan) for s in used], dtype=np.float64)
            ar_rho, n_ar = finite_spearman(ar_vec, prot, MIN_N)
        key = (rec.xena_antibody, rec.cancer)
        extra = f1s_key.loc[key] if key in f1s_key.index else None
        row = {
            "xena_antibody": rec.xena_antibody,
            "gene_symbol": rec.gene_symbol,
            "cancer": rec.cancer,
            "n_s0": n_s0,
            "n_ar": n_ar,
            "actual": actual,
            "s0": actual,
            "ar": ar_rho,
            "delta": ar_rho - actual if np.isfinite(ar_rho) and np.isfinite(actual) else np.nan,
            "is_adc": rec.gene_symbol in RPPA_ADC,
            "f0": rec.f0,
            "f0_half_80": rec.f0_half_80,
            "f0_half_90": rec.f0_half_90,
            "f1_ridge": rec.f1_ridge,
            "f1_ridge_half_80": rec.f1_ridge_half_80,
            "f1_ridge_half_90": rec.f1_ridge_half_90,
            "f1_lgbm": rec.f1_lgbm,
            "f1_lgbm_half_80": rec.f1_lgbm_half_80,
            "f1_lgbm_half_90": rec.f1_lgbm_half_90,
        }
        if extra is not None:
            row.update(
                {
                    "f1s_ridge": extra.f1s_ridge,
                    "f1s_ridge_half_80": extra.f1s_ridge_half_80,
                    "f1s_ridge_half_90": extra.f1s_ridge_half_90,
                    "f1s_lgbm": extra.f1s_lgbm,
                    "f1s_lgbm_half_80": extra.f1s_lgbm_half_80,
                    "f1s_lgbm_half_90": extra.f1s_lgbm_half_90,
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def build_procan_pairs(bundle):
    f0 = pd.read_csv(RES / "frozen_forecasts_procan.tsv", sep="\t")
    f1s = pd.read_csv(RES / "frozen_forecasts_procan_F1s.tsv", sep="\t")
    ar = pd.read_parquet(RES / "frozen_ar_procan.parquet")
    f0["is_adc"] = boolish(f0["is_adc"])
    proteins = bundle["procan_proteins"]
    protein_index = {p: i for i, p in enumerate(proteins["protein_id"].tolist())}
    sample_index = {s: i for i, s in enumerate(bundle["procan_samples"])}
    mapping = bundle["procan_mapping"]
    sidm_of = dict(zip(mapping["Project_Identifier"], mapping["SIDM"]))
    project_of = dict(zip(mapping["SIDM"], mapping["Project_Identifier"]))
    lineage_of_sidm = dict(zip(mapping["SIDM"], mapping["Tissue_type"]))
    models = {}
    for sidm, lin in lineage_of_sidm.items():
        if lin not in bundle["lineages"]:
            continue
        project = project_of.get(sidm)
        if project not in sample_index:
            continue
        models.setdefault(lin, []).append(sidm)
    raw_sidms = {sidm for sidm, _gene in bundle["procan_raw"]}
    models = {lin: [s for s in sids if s in raw_sidms] for lin, sids in models.items()}
    f1s_key = f1s.set_index(["protein_id", "lineage"])
    grouped = {key: sub for key, sub in ar.groupby(["protein_id", "lineage"], sort=False)}
    rows = []
    for rec in f0.itertuples(index=False):
        sidms = models.get(rec.lineage, [])
        projects = [project_of[s] for s in sidms if s in project_of and project_of[s] in sample_index]
        if rec.protein_id not in protein_index:
            prot = np.array([])
            used = []
        else:
            col = protein_index[rec.protein_id]
            used = sidms
            prot = np.array(
                [bundle["procan"][sample_index[project_of[s]], col] if project_of.get(s) in sample_index else np.nan for s in used],
                dtype=np.float64,
            )
        rna = np.array([bundle["procan_raw"].get((s, rec.gene_symbol), np.nan) for s in used], dtype=np.float64)
        actual, n_s0 = finite_spearman(rna, prot, MIN_N)
        sub = grouped.get((rec.protein_id, rec.lineage))
        if sub is None or not used:
            ar_rho, n_ar = np.nan, 0
        else:
            pred_map = dict(zip(sub["model_id"], sub["ar_pred"]))
            ar_vec = np.array([pred_map.get(s, np.nan) for s in used], dtype=np.float64)
            ar_rho, n_ar = finite_spearman(ar_vec, prot, MIN_N)
        key = (rec.protein_id, rec.lineage)
        extra = f1s_key.loc[key] if key in f1s_key.index else None
        row = {
            "protein_id": rec.protein_id,
            "gene_symbol": rec.gene_symbol,
            "lineage": rec.lineage,
            "n_s0": n_s0,
            "n_ar": n_ar,
            "actual": actual,
            "s0": actual,
            "ar": ar_rho,
            "delta": ar_rho - actual if np.isfinite(ar_rho) and np.isfinite(actual) else np.nan,
            "is_adc": bool(rec.is_adc),
            "f0": rec.f0,
            "f0_half_80": rec.f0_half_80,
            "f0_half_90": rec.f0_half_90,
            "f1_ridge": rec.f1_ridge,
            "f1_ridge_half_80": rec.f1_ridge_half_80,
            "f1_ridge_half_90": rec.f1_ridge_half_90,
            "f1_lgbm": rec.f1_lgbm,
            "f1_lgbm_half_80": rec.f1_lgbm_half_80,
            "f1_lgbm_half_90": rec.f1_lgbm_half_90,
        }
        if extra is not None:
            row.update(
                {
                    "f1s_ridge": extra.f1s_ridge,
                    "f1s_ridge_half_80": extra.f1s_ridge_half_80,
                    "f1s_ridge_half_90": extra.f1s_ridge_half_90,
                    "f1s_lgbm": extra.f1s_lgbm,
                    "f1s_lgbm_half_80": extra.f1s_lgbm_half_80,
                    "f1s_lgbm_half_90": extra.f1s_lgbm_half_90,
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


PREDICTORS = (
    ("f0", "f0", "f0_half_80", "f0_half_90"),
    ("f1_ridge", "f1_ridge", "f1_ridge_half_80", "f1_ridge_half_90"),
    ("f1_lgbm", "f1_lgbm", "f1_lgbm_half_80", "f1_lgbm_half_90"),
    ("f1s_ridge", "f1s_ridge", "f1s_ridge_half_80", "f1s_ridge_half_90"),
    ("f1s_lgbm", "f1s_lgbm", "f1s_lgbm_half_80", "f1s_lgbm_half_90"),
)


def write_e1(pairs, tag):
    units = unit_mean_table(pairs, "xena_antibody")
    units.to_csv(RES / f"confirm_fake_{tag}_e1_units.tsv", sep="\t", index=False)
    summary = test_record(units["mean_delta"] if len(units) else np.array([]))
    summary["metric"] = "E1"
    summary["subset"] = tag
    by = context_summary(pairs, "cancer")
    by.to_csv(RES / f"confirm_fake_{tag}_e1_by_cancer.tsv", sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(RES / f"confirm_fake_{tag}_e1_summary.tsv", sep="\t", index=False)
    return units, summary, by


def write_e2(pairs, tag):
    units = unit_mean_table(pairs, "protein_id")
    units.to_csv(RES / f"confirm_fake_{tag}_e2_units.tsv", sep="\t", index=False)
    summary = test_record(units["mean_delta"] if len(units) else np.array([]))
    summary["metric"] = "E2"
    summary["subset"] = tag
    by = context_summary(pairs, "lineage")
    by.to_csv(RES / f"confirm_fake_{tag}_e2_by_lineage.tsv", sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(RES / f"confirm_fake_{tag}_e2_summary.tsv", sep="\t", index=False)
    return units, summary, by


def write_reliability(pairs, tag, cohort, with_within):
    context_col = "cancer" if cohort == "rppa" else "lineage"
    unit_col = "xena_antibody" if cohort == "rppa" else "protein_id"
    overall = []
    within_sum = []
    centered_rows = []
    bias_rows = []
    for name, pred, h80, h90 in PREDICTORS:
        if pred not in pairs.columns:
            continue
        block = reliability_block(pairs, pred, h80, h90)
        block.update({"subset": tag, "model": name, "item": "pooled"})
        overall.append(block)
        by, summary = within_context(pairs, context_col, pred)
        by.to_csv(RES / f"confirm_fake_{tag}_{name}_by_context.tsv", sep="\t", index=False)
        summary.update({"subset": tag, "model": name, "item": "within_context"})
        within_sum.append(summary)
        bias = bias_table(pairs, context_col, pred)
        bias["model"] = name
        bias["subset"] = tag
        bias_rows.append(bias)
        if with_within:
            _units, usum = within_unit(pairs, unit_col, pred)
            usum.update({"subset": tag, "model": name, "item": "within_unit"})
            within_sum.append(usum)
            cen = centered(pairs, unit_col, pred)
            cen.update({"subset": tag, "model": name, "item": "centered"})
            centered_rows.append(cen)
    pd.DataFrame(overall).to_csv(RES / f"confirm_fake_{tag}_reliability.tsv", sep="\t", index=False)
    pd.DataFrame(within_sum).to_csv(RES / f"confirm_fake_{tag}_within.tsv", sep="\t", index=False)
    if centered_rows:
        pd.DataFrame(centered_rows).to_csv(RES / f"confirm_fake_{tag}_centered.tsv", sep="\t", index=False)
    if bias_rows:
        pd.concat(bias_rows, ignore_index=True).to_csv(RES / f"confirm_fake_{tag}_bias.tsv", sep="\t", index=False)
    return overall


def write_examples(named_frames):
    chunks = []
    for name, frame in named_frames:
        if frame is None or len(frame) == 0:
            continue
        head = frame.head(3).copy()
        head.insert(0, "table", name)
        chunks.append(head)
    if chunks:
        pd.concat(chunks, ignore_index=True).to_csv(RES / "confirm_fake_examples.tsv", sep="\t", index=False)


def run(bundle):
    print("E metrics", flush=True)
    rppa = build_rppa_pairs(bundle)
    procan = build_procan_pairs(bundle)
    rppa.to_csv(RES / "confirm_fake_rppa_pairs.tsv", sep="\t", index=False)
    procan.to_csv(RES / "confirm_fake_procan_pairs.tsv", sep="\t", index=False)
    e1_units, e1_summary, e1_by = write_e1(rppa, "all")
    e2_units, e2_summary, e2_by = write_e2(procan, "all")
    e3 = write_reliability(rppa, "e3", "rppa", True)
    e4 = write_reliability(procan, "e4", "procan", False)
    e5_rppa = rppa[rppa["gene_symbol"].isin(RPPA_ADC)].copy()
    e5_procan = procan[boolish(procan["is_adc"])].copy()
    write_e1(e5_rppa, "e5")
    write_e2(e5_procan, "e5")
    write_reliability(e5_rppa, "e5_rppa", "rppa", True)
    write_reliability(e5_procan, "e5_procan", "procan", False)
    write_examples(
        [
            ("e1_summary", pd.DataFrame([e1_summary])),
            ("e1_units", e1_units),
            ("e1_by_cancer", e1_by),
            ("e2_summary", pd.DataFrame([e2_summary])),
            ("e2_units", e2_units),
            ("e2_by_lineage", e2_by),
            ("e3_reliability", pd.DataFrame(e3)),
            ("e4_reliability", pd.DataFrame(e4)),
            ("rppa_pairs", rppa),
            ("procan_pairs", procan),
        ]
    )
    pd.DataFrame(
        [
            {"item": "n_rppa_pairs", "value": int(len(rppa))},
            {"item": "n_procan_pairs", "value": int(len(procan))},
            {"item": "n_e5_rppa_pairs", "value": int(len(e5_rppa))},
            {"item": "n_e5_procan_pairs", "value": int(len(e5_procan))},
            {"item": "n_e5_procan_genes", "value": int(e5_procan["gene_symbol"].nunique()) if len(e5_procan) else 0},
        ]
    ).to_csv(RES / "confirm_fake_counts.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unlock", action="store_true")
    parser.add_argument("--fake", action="store_true")
    args = parser.parse_args()
    if args.unlock and args.fake:
        raise SystemExit("pass only one of --unlock or --fake")
    if not args.unlock and not args.fake:
        raise SystemExit("refusing to run without --unlock or --fake")
    t0 = time.perf_counter()
    errors = []
    try:
        if args.unlock:
            bundle = load_real_bundle()
        else:
            bundle = make_fake()
        run(bundle)
    except Exception as exc:
        errors.append(f"{type(exc).__name__}: {exc}")
        pd.DataFrame([{"item": "error", "value": errors[-1]}]).to_csv(RES / "confirm_fake_run.tsv", sep="\t", index=False)
        raise
    elapsed = time.perf_counter() - t0
    pd.DataFrame(
        [
            {"item": "elapsed_sec", "value": elapsed},
            {"item": "n_errors", "value": len(errors)},
            {"item": "mode", "value": "fake" if args.fake else "unlock"},
        ]
    ).to_csv(RES / "confirm_fake_run.tsv", sep="\t", index=False)
    print(f"confirm done {elapsed:.1f}s errors {len(errors)}", flush=True)


if __name__ == "__main__":
    main()
