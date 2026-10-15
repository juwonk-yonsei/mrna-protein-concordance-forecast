"""Stage-5 post-hoc sensitivity, H2 decomposition, and training-context distance.

Reads the unlocked RPPA matrix for the duplicate-sample check. Does not write
into results/stage4/confirm or the frozen files.
"""
import gzip
import importlib.util
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from scripts.b2_lib import RAW, RES, ROOT
from scripts.b3_lib import percentile_map, predict_ar_matrix

LABEL = "사후"
OUT = RES / "stage5" / "posthoc"
ATLAS = RES / "stage5" / "atlas"
MIN_N = 30


def load_09():
    spec = importlib.util.spec_from_file_location("confirm_b3", ROOT / "scripts" / "09_confirm.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write(path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = frame.copy()
    frame.insert(0, "label", LABEL)
    frame.to_csv(path, sep="\t", index=False)
    print("wrote", path, len(frame), flush=True)


def finite_rho(c09, a, b):
    value, _n = c09.finite_spearman(a, b, MIN_N)
    return value


def rho_s0(c09, raw, prot):
    z = c09.within_group_z(raw)
    z = np.where(np.isfinite(z), z, 0.0)
    if raw.size and not np.isfinite(raw).any():
        z = np.zeros(raw.shape, dtype=np.float64)
    return finite_rho(c09, z, prot)


def profile_spearman(c09, left, right, genes):
    a = np.array([left.get(g, np.nan) for g in genes], dtype=np.float64)
    b = np.array([right.get(g, np.nan) for g in genes], dtype=np.float64)
    value, _n = c09.finite_spearman(a, b, 3)
    return value


def is_blood(name):
    text = str(name).lower()
    return any(token in text for token in ("leukemia", "lymphoma", "myeloma", "haematopoietic", "hematopoietic"))


def sample_type(sample):
    parts = str(sample).split("-")
    if len(parts) < 4:
        return ""
    return parts[3][:2]


def wilcox_mean(values):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    rec = {"n": int(values.size), "mean": np.nan, "pvalue": np.nan, "statistic": np.nan, "error": ""}
    if values.size == 0:
        rec["error"] = "n=0"
        return rec
    rec["mean"] = float(values.mean())
    try:
        tested = wilcoxon(values, alternative="greater", zero_method="wilcox")
        rec["statistic"] = float(tested.statistic)
        rec["pvalue"] = float(tested.pvalue)
    except ValueError as exc:
        rec["error"] = str(exc)
    return rec


def h2_tables():
    pairs = pd.read_csv(RES / "stage4" / "confirm" / "procan_pairs.tsv", sep="\t")
    model = pd.read_parquet(RES / "ar_final_models.parquet", columns=["gene_symbol", "lambda"])
    pset = set(g for g in (RES / "P_genes.txt").read_text().split() if g)
    use = pairs[np.isfinite(pairs["delta"])].copy()
    use = use.merge(model, on="gene_symbol", how="left")
    use["in_P"] = use["gene_symbol"].isin(pset)
    rows = []
    for (lineage, lam), sub in use.groupby(["lineage", "lambda"], dropna=False, sort=True):
        rows.append({"table": "lineage_x_lambda", "lineage": lineage, "lambda": lam, "in_P": "", "n": int(len(sub)), "median_delta": float(sub["delta"].median())})
    for (lineage, flag), sub in use.groupby(["lineage", "in_P"], sort=True):
        rows.append({"table": "lineage_x_in_P", "lineage": lineage, "lambda": np.nan, "in_P": bool(flag), "n": int(len(sub)), "median_delta": float(sub["delta"].median())})
    write(OUT / "posthoc_22_h2_decomposition.tsv", pd.DataFrame(rows))
    return use


def median_profile(mat, genes, ids):
    block = mat.reindex(index=genes, columns=ids).to_numpy(dtype=np.float64)
    with np.errstate(all="ignore"):
        med = np.nanmedian(block, axis=1)
        sd = np.nanstd(block, axis=1, ddof=0)
    med_map = percentile_map({g: float(v) if np.isfinite(v) else np.nan for g, v in zip(genes, med)})
    sd_map = percentile_map({g: float(v) if np.isfinite(v) else np.nan for g, v in zip(genes, sd)})
    return {g: med_map[g][1] for g in genes}, {g: sd_map[g][1] for g in genes}, {g: med_map[g][1] for g in genes}


def sd_percentile(stage3, mat, genes, ids):
    _mu, sd = stage3.gene_mean_sd(mat, ids)
    mapped = percentile_map({g: sd.get(g, np.nan) for g in genes})
    return {g: mapped[g][1] for g in genes}


def predict_ar(stage3, mat, ids, genes, gindex, keep_pos, self_index, feat_index, b, lams, coefs):
    needed = []
    pos_of = {}
    for j in keep_pos:
        for g_i in [int(self_index[j]), *[int(v) for v in feat_index[j]]]:
            gene = genes[g_i]
            if gene not in pos_of:
                pos_of[gene] = len(needed)
                needed.append(gene)
    z_need = stage3.z_for_genes(mat, ids, needed)
    z_full = np.zeros((len(ids), len(genes)), dtype=np.float64)
    for gene, col in pos_of.items():
        z_full[:, gindex[gene]] = z_need[:, col]
    z_full = np.nan_to_num(z_full, nan=0.0)
    pred = predict_ar_matrix(
        z_full,
        self_index[keep_pos],
        feat_index[keep_pos].T,
        b[keep_pos],
        lams[keep_pos],
        np.asarray(coefs[keep_pos], dtype=np.float64).T,
    )
    return {genes[int(self_index[j])]: pred[:, t] for t, j in enumerate(keep_pos)}


def rppa_values():
    from scripts.loader import load_locked_matrix

    frame = load_locked_matrix("tcga_rppa", unlock=True)
    with gzip.open(RAW / "tcga" / "TCGA-RPPA-pancan-clean.xena.gz", "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")[1:]
    if len(frame.columns) != len(header):
        raise RuntimeError("RPPA width does not match the header")
    frame.columns = header
    columns_of = {}
    for j, sample in enumerate(header):
        columns_of.setdefault(sample, []).append(j)
    return frame.index.astype(str).tolist(), columns_of, frame.to_numpy(dtype=np.float64)


def primary_by_disease(header_counts):
    ph = pd.read_csv(RAW / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", dtype=str, usecols=["sample", "_primary_disease"])
    by = {}
    for sample, disease in zip(ph["sample"], ph["_primary_disease"]):
        if header_counts.get(sample, 0) != 1:
            continue
        if sample_type(sample) not in {"01", "03"}:
            continue
        if disease is None or disease == "" or disease == "nan":
            continue
        by.setdefault(disease, [])
        if sample not in by[disease]:
            by[disease].append(sample)
    return by


def sensitivity(c09, stage3, mat, genes, pairs, h1):
    antibodies, columns_of, values = rppa_values()
    ab_index = {a: i for i, a in enumerate(antibodies)}
    base, _path = stage3.disease_samples()
    ph = pd.read_csv(RAW / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", dtype=str, usecols=["sample", "_primary_disease"])
    disease_of = dict(zip(ph["sample"], ph["_primary_disease"]))
    dups = pd.read_csv(RES / "rppa_duplicate_sample_ids.tsv", sep="\t", dtype=str)["sample_id"].tolist()
    with gzip.open(RAW / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz", "rt", encoding="utf-8", errors="replace") as handle:
        rna_header = handle.readline().rstrip("\n").split("\t")[1:]
    from collections import Counter

    rna_counts = Counter(rna_header)
    audit = []
    extra = []
    for sample in dups:
        disease = disease_of.get(sample, "")
        in_list = sample in set(base.get(disease, []))
        add = disease in base and rna_counts[sample] == 1 and not in_list
        audit.append(
            {
                "sample_id": sample,
                "n_rppa_columns": len(columns_of.get(sample, [])),
                "in_rna_once": int(rna_counts[sample] == 1),
                "phenotype_disease": disease,
                "disease_in_31": int(disease in base),
                "already_in_preregistered_samples": int(in_list),
                "added_in_first_and_mean": int(add),
            }
        )
        if add:
            extra.append((disease, sample))
    write(OUT / "posthoc_21_duplicate_ids.tsv", pd.DataFrame(audit))
    rule = pd.DataFrame(
        [
            {
                "item": "preregistered_sample_rule",
                "value": "disease_samples drops the 10 duplicate ids before intersecting RNA. attach_rppa_s0 requests only those ids. setdefault would keep the first RPPA column if an id were requested, but these ids are not requested, so neither column enters actual rho, rho_S0, AR, or B_sd.",
            },
            {"item": "n_duplicate_ids_in_preregistered_lists", "value": str(int(sum(r["already_in_preregistered_samples"] for r in audit)))},
            {"item": "n_ids_added_in_first_and_mean", "value": str(int(sum(r["added_in_first_and_mean"] for r in audit)))},
            {"item": "first", "value": "one row per id, protein value from the first header column"},
            {"item": "mean", "value": "one row per id, protein value is the mean of the finite values in every header column with that id"},
            {"item": "drop", "value": "same sample lists as the preregistered confirmation"},
        ]
    )
    write(OUT / "posthoc_21_duplicate_rule.tsv", rule)
    first = {k: list(v) for k, v in base.items()}
    for disease, sample in extra:
        first[disease].append(sample)
    scenarios = {"drop": base, "first": first, "mean": first}
    model, self_index, feat_index, coefs = stage3.models_from_parquet(RES / "ar_final_models.parquet", genes)
    b = model["b"].to_numpy(dtype=np.float64)
    lams = model["lambda"].to_numpy(dtype=np.float64)
    gindex = {g: i for i, g in enumerate(genes)}
    symbol_at = {i: g for i, g in enumerate(model["gene_symbol"])}
    keep_pos = [i for i, g in symbol_at.items() if g in set(h1["gene_symbol"])]
    frozen = pd.read_parquet(RES / "frozen_ar_tcga.parquet", columns=["cancer", "sample_id", "gene_symbol", "ar_pred"])
    rows = []
    ref = pairs.set_index(["xena_antibody", "cancer"])
    h1_genes = dict(zip(h1["xena_antibody"], h1["gene_symbol"]))
    for mode, lists in (("drop", base), ("first", first), ("mean", first)):
        actual_gap = []
        ar_gap = []
        bsd_gap = []
        built = []
        for cancer, ids in lists.items():
            pred = predict_ar(stage3, mat, ids, genes, gindex, keep_pos, self_index, feat_index, b, lams, coefs)
            bsd = sd_percentile(stage3, mat, genes, ids)
            if mode == "drop":
                frozen_c = frozen[frozen["cancer"] == cancer]
                for gene, vec in pred.items():
                    sub = frozen_c[frozen_c["gene_symbol"] == gene].set_index("sample_id")
                    for sample, value in zip(ids, vec):
                        if sample in sub.index:
                            ar_gap.append(abs(float(value) - float(sub.loc[sample, "ar_pred"] if np.ndim(sub.loc[sample, "ar_pred"]) == 0 else sub.loc[sample, "ar_pred"].iloc[0])))
            part = pairs[pairs["cancer"] == cancer]
            for rec in part.itertuples(index=False):
                if rec.xena_antibody not in ab_index:
                    continue
                cols = []
                for sample in ids:
                    pos = columns_of[sample]
                    if mode == "mean" and len(pos) > 1:
                        cols.append(float(np.nanmean(values[ab_index[rec.xena_antibody], pos])))
                    else:
                        cols.append(float(values[ab_index[rec.xena_antibody], pos[0]]))
                prot = np.asarray(cols, dtype=np.float64)
                if rec.gene_symbol in mat.index:
                    raw = pd.to_numeric(mat.loc[rec.gene_symbol, ids], errors="coerce").to_numpy(dtype=np.float64)
                else:
                    raw = np.full(len(ids), np.nan)
                actual = finite_rho(c09, raw, prot)
                s0 = rho_s0(c09, raw, prot)
                ar_vec = pred.get(rec.gene_symbol)
                if ar_vec is None:
                    ar = np.nan
                else:
                    ar = finite_rho(c09, ar_vec, prot)
                delta = ar - s0 if np.isfinite(ar) and np.isfinite(s0) else np.nan
                key = (rec.xena_antibody, cancer)
                old = ref.loc[key]
                if mode == "drop":
                    if np.isfinite(actual) and np.isfinite(old.actual):
                        actual_gap.append(abs(actual - float(old.actual)))
                    if np.isfinite(ar) and np.isfinite(old.ar):
                        bsd_gap.append(0.0)
                    if np.isfinite(bsd.get(rec.gene_symbol, np.nan)) and np.isfinite(old.b_sd):
                        bsd_gap[-1] = abs(bsd[rec.gene_symbol] - float(old.b_sd)) if bsd_gap else abs(bsd[rec.gene_symbol] - float(old.b_sd))
                        if not (np.isfinite(ar) and np.isfinite(old.ar)):
                            pass
                built.append(
                    {
                        "xena_antibody": rec.xena_antibody,
                        "gene_symbol": rec.gene_symbol,
                        "cancer": cancer,
                        "actual": actual,
                        "rho_s0": s0,
                        "ar": ar,
                        "delta": delta,
                        "f1s_lgbm": float(old.f1s_lgbm),
                        "b_sd": bsd.get(rec.gene_symbol, np.nan),
                        "in_h1": rec.xena_antibody in h1_genes,
                    }
                )
            print(mode, cancer, len(ids), flush=True)
        built = pd.DataFrame(built)
        if mode == "drop":
            # Recompute the B_sd gap without the placeholder above.
            bsd_gap = []
            ar_rho_gap = []
            for rec in built.itertuples(index=False):
                old = ref.loc[(rec.xena_antibody, rec.cancer)]
                if np.isfinite(rec.b_sd) and np.isfinite(old.b_sd):
                    bsd_gap.append(abs(float(rec.b_sd) - float(old.b_sd)))
                if np.isfinite(rec.ar) and np.isfinite(old.ar):
                    ar_rho_gap.append(abs(float(rec.ar) - float(old.ar)))
            print(
                "drop match",
                "actual", max(actual_gap) if actual_gap else None,
                "ar_rho", max(ar_rho_gap) if ar_rho_gap else None,
                "bsd", max(bsd_gap) if bsd_gap else None,
                "ar_pred", max(ar_gap) if ar_gap else None,
                flush=True,
            )
            if not actual_gap or max(actual_gap) > 1e-6:
                raise RuntimeError("drop scenario actual rho does not match the preregistered pairs")
            if not ar_gap or max(ar_gap) > 1e-2:
                raise RuntimeError("drop scenario AR predictions do not match the frozen file")
        h1_means = []
        for antibody, sub in built[built["in_h1"]].groupby("xena_antibody", sort=True):
            d = sub["delta"].to_numpy(dtype=np.float64)
            d = d[np.isfinite(d)]
            if d.size:
                h1_means.append(float(d.mean()))
        tested = wilcox_mean(h1_means)
        use = built[np.isfinite(built["actual"]) & np.isfinite(built["f1s_lgbm"]) & np.isfinite(built["b_sd"])]
        s_f, _ = c09.finite_spearman(use["f1s_lgbm"], use["actual"], 3)
        s_b, _ = c09.finite_spearman(use["b_sd"], use["actual"], 3)
        rows.append(
            {
                "method": mode,
                "n_samples": int(sum(len(v) for v in lists.values())),
                "h1_n": tested["n"],
                "h1_mean": tested["mean"],
                "h1_statistic": tested["statistic"],
                "h1_pvalue": tested["pvalue"],
                "h3_n_pairs": int(len(use)),
                "S_F": s_f,
                "S_B": s_b,
                "diff": s_f - s_b if np.isfinite(s_f) and np.isfinite(s_b) else np.nan,
            }
        )
        del built
    hyp = pd.read_csv(RES / "stage4" / "confirm" / "hypotheses.tsv", sep="\t")
    h1row = hyp[hyp["hypothesis"] == "H1"].iloc[0]
    h3row = hyp[hyp["hypothesis"] == "H3"].iloc[0]
    rows.append(
        {
            "method": "preregistered",
            "n_samples": int(sum(len(v) for v in base.values())),
            "h1_n": int(h1row.n_units),
            "h1_mean": float(h1row["mean"]),
            "h1_statistic": float(h1row.statistic),
            "h1_pvalue": float(h1row.pvalue),
            "h3_n_pairs": int(h3row.n_pairs),
            "S_F": float(h3row.S_F),
            "S_B": float(h3row.S_B),
            "diff": float(h3row["diff"]),
        }
    )
    write(OUT / "posthoc_21_sensitivity.tsv", pd.DataFrame(rows))
    del values


def residual_rows(c09, extras, sd_pct, deltas, cohort, context):
    fs = []
    ds = []
    for gene, delta in deltas:
        extra = extras.get(gene)
        if extra is None or not np.isfinite(delta):
            continue
        flags = []
        for name in extra:
            pct = sd_pct.get(name, np.nan)
            flags.append(bool(np.isfinite(pct) and pct < 10))
        if len(flags) != 50:
            continue
        fs.append(float(np.mean(flags)))
        ds.append(float(delta))
    return {
        "cohort": cohort,
        "context": context,
        "n": len(fs),
        "spearman_f_delta": profile_spearman(c09, {}, {}, []) if False else (c09.finite_spearman(np.asarray(fs), np.asarray(ds), 3)[0] if fs else np.nan),
        "median_f": float(np.median(fs)) if fs else np.nan,
    }


def cptac_profiles(genes):
    paired = pd.read_csv(RES / "cptac_paired_samples.tsv", sep="\t", dtype=str)
    paired = paired[paired["tissue"] == "tumor"]
    from scripts.b2_lib import PROC

    profiles = {}
    for cancer in ["BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC"]:
        rna = pd.read_parquet(PROC / "cptac" / f"{cancer}__tumor__rna.parquet")
        sids = [s for s in paired.loc[paired["cancer"] == cancer, "sample_id"] if s in rna.columns]
        block = rna.reindex(index=genes, columns=sids).to_numpy(dtype=np.float64)
        with np.errstate(all="ignore"):
            med = np.nanmedian(block, axis=1)
        mapped = percentile_map({g: float(v) if np.isfinite(v) else np.nan for g, v in zip(genes, med)})
        profiles[cancer] = {g: mapped[g][1] for g in genes}
        print("cptac profile", cancer, len(sids), flush=True)
        del rna, block
    return profiles


def max_similarity(c09, profile, profiles, genes, skip=None):
    vals = []
    for name, other in profiles.items():
        if skip is not None and name == skip:
            continue
        vals.append(profile_spearman(c09, profile, other, genes))
    finite = [v for v in vals if np.isfinite(v)]
    return float(np.max(finite)) if finite else np.nan


def choose_tau(c09, cptac, genes, ar_pairs):
    s_of = {cancer: max_similarity(c09, cptac[cancer], cptac, genes, skip=cancer) for cancer in cptac}
    ordered = sorted(v for v in s_of.values() if np.isfinite(v))
    mids = [(ordered[i] + ordered[i + 1]) / 2 for i in range(len(ordered) - 1)]

    def score(tau):
        meds = []
        n_ar = 0
        for cancer, s in s_of.items():
            sub = ar_pairs[ar_pairs["heldout_cancer"] == cancer]
            both = sub[np.isfinite(sub["s0"]) & np.isfinite(sub["ar"])]
            use_ar = np.isfinite(s) and s >= tau
            meds.append(float((both["ar"] if use_ar else both["s0"]).median()))
            n_ar += int(use_ar)
        return float(np.mean(meds)), n_ar

    always_ar, _ = score(-1e9)
    always_s0, _ = score(1e9)
    rows = [{"tau": -1e9, "role": "always_AR", "score": always_ar, "n_contexts_AR": 10}]
    best = None
    for tau in mids:
        value, n_ar = score(tau)
        rows.append({"tau": tau, "role": "candidate", "score": value, "n_contexts_AR": n_ar})
        if 0 < n_ar < 10 and value > always_ar and value > always_s0:
            if best is None or value > best[0] or (value == best[0] and tau < best[1]):
                best = (value, tau, n_ar)
    rows.append({"tau": 1e9, "role": "always_S0", "score": always_s0, "n_contexts_AR": 0})
    decision = pd.DataFrame(rows)
    decision["identified"] = best is not None
    decision["chosen_tau"] = best[1] if best else np.nan
    write(OUT / "posthoc_22_tau_development.tsv", decision)
    pd.DataFrame([{"context": k, "s_vs_other_nine": v} for k, v in s_of.items()]).assign(label=LABEL).to_csv(
        OUT / "posthoc_22_tau_cptac_s.tsv", sep="\t", index=False
    )
    return best[1] if best else None, s_of


def apply_tau(c09, tau, context_s, rppa_pairs, pro_pairs):
    if tau is None:
        write(OUT / "posthoc_22_tau_applied.tsv", pd.DataFrame([{"status": "not_identified"}]))
        return
    rows = []
    for platform, frame, unit, context in (
        ("rppa", rppa_pairs, "xena_antibody", "cancer"),
        ("procan", pro_pairs, "protein_id", "lineage"),
    ):
        use = frame[np.isfinite(frame["delta"])].copy()
        means = []
        for _name, sub in use.groupby(unit, sort=True):
            vals = []
            for rec in sub.itertuples(index=False):
                s = context_s.get(getattr(rec, context), np.nan)
                vals.append(float(rec.delta) if np.isfinite(s) and s >= tau else 0.0)
            if vals:
                means.append(float(np.mean(vals)))
        tested = wilcox_mean(means)
        tested.update({"platform": platform, "tau": tau, "status": "applied_from_cptac_only"})
        rows.append(tested)
    write(OUT / "posthoc_22_tau_applied.tsv", pd.DataFrame(rows))


def plot_similarity(frame):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for cohort, marker in (("TCGA", "o"), ("ProCan", "s")):
        sub = frame[frame["cohort"] == cohort]
        ax.scatter(sub["s"], sub["median_delta"], s=28, marker=marker, label=cohort)
    blood = frame[frame["is_blood"] == 1]
    if len(blood):
        ax.scatter(blood["s"], blood["median_delta"], s=64, facecolors="none", edgecolors="black", linewidths=1.0, label="blood")
    ax.set_xlabel("s")
    ax.set_ylabel("median delta")
    ax.legend(frameon=False)
    fig.tight_layout()
    path = OUT / "posthoc_22_similarity.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path, flush=True)


def procan_blocks(genes):
    import csv

    ar = pd.read_parquet(RES / "frozen_ar_procan.parquet", columns=["lineage", "model_id"])
    by_lin = {lin: sub["model_id"].drop_duplicates().tolist() for lin, sub in ar.groupby("lineage", sort=True)}
    model_ids = list(dict.fromkeys(s for ids in by_lin.values() for s in ids))
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
            if parts[i_model] not in mindex or parts[i_gene] not in gindex:
                continue
            r = mindex[parts[i_model]]
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
    print("procan rna", len(model_ids), flush=True)
    return by_lin, mindex, np.asarray(arr, dtype=np.float64)


def main():
    print("thread env", {k: os.environ.get(k, "") for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, flush=True)
    h2_tables()
    c09 = load_09()
    stage3 = c09.load_stage3()
    genes = [g for g in (RES / "G_RNA_genes.txt").read_text().split() if g]
    pairs = pd.read_csv(RES / "stage4" / "confirm" / "rppa_pairs.tsv", sep="\t")
    h1 = pd.read_csv(RES / "stage4" / "confirm" / "h1_units.tsv", sep="\t")
    base, _path = stage3.disease_samples()
    with gzip.open(RAW / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz", "rt", encoding="utf-8", errors="replace") as handle:
        rna_header = handle.readline().rstrip("\n").split("\t")[1:]
    from collections import Counter

    primary = primary_by_disease(Counter(rna_header))
    dups = set(pd.read_csv(RES / "rppa_duplicate_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    load_lists = {k: list(v) for k, v in primary.items()}
    for disease, ids in base.items():
        load_lists.setdefault(disease, [])
        for sample in ids:
            if sample not in load_lists[disease]:
                load_lists[disease].append(sample)
    # Duplicate ids are not in `base`. Add them when they are in the RNA header so the sensitivity can see them.
    ph = pd.read_csv(RAW / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", dtype=str, usecols=["sample", "_primary_disease"])
    disease_of = dict(zip(ph["sample"], ph["_primary_disease"]))
    rna_once = {s for s, n in Counter(rna_header).items() if n == 1}
    for sample in dups:
        disease = disease_of.get(sample, "")
        if disease in load_lists and sample in rna_once and sample not in load_lists[disease]:
            load_lists[disease].append(sample)
    mapped = pd.read_csv(RES / "tcga_rna_id_map_b3.tsv", sep="\t", dtype=str)
    print("reading TCGA RNA", sum(len(v) for v in load_lists.values()), flush=True)
    mat = stage3.load_tcga_matrix(mapped, load_lists)
    sensitivity(c09, stage3, mat, genes, pairs, h1)
    model = pd.read_parquet(RES / "ar_final_models.parquet", columns=["gene_symbol", "extra_genes"])
    extras = {r.gene_symbol: list(r.extra_genes) for r in model.itertuples(index=False)}
    confirm_profiles = {}
    atlas_rows = []
    residual = []
    for disease, ids in base.items():
        med_pct, sd_pct, _same = median_profile(mat, genes, ids)
        confirm_profiles[disease] = med_pct
        deltas = [(r.gene_symbol, r.delta) for r in pairs[pairs["cancer"] == disease].itertuples(index=False)]
        residual.append(residual_rows(c09, extras, sd_percentile(stage3, mat, genes, ids), deltas, "TCGA", disease))
        print("tcga context", disease, flush=True)
    for disease, ids in primary.items():
        med_pct, sd_pct, _same = median_profile(mat, genes, ids)
        for gene in genes:
            atlas_rows.append(
                {
                    "disease": disease,
                    "gene_symbol": gene,
                    "n_samples": len(ids),
                    "median_percentile": med_pct[gene],
                    "sd_percentile": sd_pct[gene],
                    "is_blood": int(is_blood(disease)),
                    "sample_set": "primary_01_or_03",
                }
            )
        print("atlas rna", disease, len(ids), flush=True)
    del mat
    write(ATLAS / "tcga_rna_percentiles.tsv", pd.DataFrame(atlas_rows))
    cptac = cptac_profiles(genes)
    tau, _cptac_s = choose_tau(c09, cptac, genes, pd.read_csv(RES / "ar_pairs.tsv", sep="\t"))
    by_lin, mindex, arr = procan_blocks(genes)
    pro_pairs = pd.read_csv(RES / "stage4" / "confirm" / "procan_pairs.tsv", sep="\t")
    pro_profiles = {}
    gindex = {g: i for i, g in enumerate(genes)}
    for lin, ids in by_lin.items():
        rows = [mindex[m] for m in ids]
        block = arr[rows]
        with np.errstate(all="ignore"):
            med = np.nanmedian(block, axis=0)
            sd = np.nanstd(block, axis=0, ddof=0)
        med_map = percentile_map({g: float(med[gindex[g]]) if np.isfinite(med[gindex[g]]) else np.nan for g in genes})
        sd_map = {g: float(sd[gindex[g]]) if np.isfinite(sd[gindex[g]]) else np.nan for g in genes}
        sd_pct = percentile_map(sd_map)
        pro_profiles[lin] = {g: med_map[g][1] for g in genes}
        deltas = [(r.gene_symbol, r.delta) for r in pro_pairs[pro_pairs["lineage"] == lin].itertuples(index=False)]
        residual.append(
            residual_rows(
                c09,
                extras,
                {g: sd_pct[g][1] for g in genes},
                deltas,
                "ProCan",
                lin,
            )
        )
        print("procan context", lin, len(ids), flush=True)
    write(OUT / "posthoc_22_residual_f.tsv", pd.DataFrame(residual))
    sim_rows = []
    for disease, profile in confirm_profiles.items():
        sub = pairs[pairs["cancer"] == disease]
        d = sub["delta"].to_numpy(dtype=np.float64)
        d = d[np.isfinite(d)]
        sim_rows.append(
            {
                "cohort": "TCGA",
                "context": disease,
                "s": max_similarity(c09, profile, cptac, genes),
                "median_delta": float(np.median(d)) if d.size else np.nan,
                "n_pairs": int(d.size),
                "is_blood": int(is_blood(disease)),
                "in_confirmation_34": 1,
                "sample_set": "rppa_rna_overlap",
            }
        )
    for lin, profile in pro_profiles.items():
        sub = pro_pairs[pro_pairs["lineage"] == lin]
        d = sub["delta"].to_numpy(dtype=np.float64)
        d = d[np.isfinite(d)]
        sim_rows.append(
            {
                "cohort": "ProCan",
                "context": lin,
                "s": max_similarity(c09, profile, cptac, genes),
                "median_delta": float(np.median(d)) if d.size else np.nan,
                "n_pairs": int(d.size),
                "is_blood": int(is_blood(lin)),
                "in_confirmation_34": 1,
                "sample_set": "ar_models",
            }
        )
    for disease, ids in primary.items():
        if disease in confirm_profiles:
            continue
        # s for cancers outside the confirmation set uses the same primary-tumor profile as the atlas.
        med_pct = {}
        sub = [r for r in atlas_rows if r["disease"] == disease]
        med_pct = {r["gene_symbol"]: r["median_percentile"] for r in sub}
        sim_rows.append(
            {
                "cohort": "TCGA",
                "context": disease,
                "s": max_similarity(c09, med_pct, cptac, genes),
                "median_delta": np.nan,
                "n_pairs": 0,
                "is_blood": int(is_blood(disease)),
                "in_confirmation_34": 0,
                "sample_set": "primary_01_or_03",
            }
        )
    sim = pd.DataFrame(sim_rows)
    write(OUT / "posthoc_22_similarity.tsv", sim)
    confirm = sim[sim["in_confirmation_34"] == 1]
    both = confirm[np.isfinite(confirm["s"]) & np.isfinite(confirm["median_delta"])]
    rho, n = c09.finite_spearman(both["s"], both["median_delta"], 3)
    write(
        OUT / "posthoc_22_similarity_summary.tsv",
        pd.DataFrame([{"n": int(len(both)), "spearman_s_median_delta": rho, "n_finite_pairs_in_spearman": n}]),
    )
    plot_similarity(confirm)
    context_s = dict(zip(confirm["context"], confirm["s"]))
    apply_tau(c09, tau, context_s, pairs, pro_pairs)
    print("scope done", flush=True)


if __name__ == "__main__":
    main()
