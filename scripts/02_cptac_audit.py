import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.common import load_paths, load_plan

PATHS = load_paths()
PLAN = load_plan()
CANCERS = PLAN["cancers"]
LOCUS = PLAN["gene_id"]["protein_coding_locus_group"]
RAW = PATHS["raw_dir"] / "cptac"
PROC = PATHS["processed_dir"] / "cptac"
RES = PATHS["results_dir"]
PROC.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)


def read_header_samples(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    dup = sorted(s for s, n in Counter(samples).items() if n > 1)
    return samples, dup


def read_matrix(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def load_hgnc():
    h = pd.read_csv(PATHS["raw_dir"] / "hgnc" / "hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
    h["ensembl_gene_id"] = h["ensembl_gene_id"].fillna("")
    by = {}
    for rec in h.to_dict("records"):
        eid = rec.get("ensembl_gene_id") or ""
        if not eid:
            continue
        by.setdefault(eid, []).append(rec)
    symbols = set(h.loc[h["locus_group"] == LOCUS, "symbol"].dropna())
    return by, symbols, h


def map_row(gene_id, by_ensembl):
    gid = str(gene_id)
    if gid.startswith("ENSG"):
        base = gid.split(".", 1)[0]
        version = gid.split(".", 1)[1] if "." in gid else ""
        version_n = int(version) if version.isdigit() else -1
        hits = by_ensembl.get(base, [])
        pc = [r for r in hits if r.get("locus_group") == LOCUS]
        if len(pc) == 1:
            return pc[0]["symbol"], "ensembl_protein_coding", version_n
        if len(pc) > 1:
            return None, "ambiguous_ensembl", version_n
        if hits:
            return None, "not_protein_coding", version_n
        return None, "unmapped", version_n
    return None, "not_ensembl", -1


def collapse(df, by_ensembl):
    mapped = [map_row(g, by_ensembl) for g in df.index]
    status = Counter(s for _, s, _ in mapped)
    keep_symbol = []
    keep_version = []
    keep_pos = []
    for i, (sym, st, ver) in enumerate(mapped):
        if sym is None:
            continue
        keep_symbol.append(sym)
        keep_version.append(ver)
        keep_pos.append(i)
    sub = df.iloc[keep_pos].copy()
    sub.insert(0, "_symbol", keep_symbol)
    sub.insert(1, "_version", keep_version)
    sub = sub.sort_values(["_symbol", "_version"], ascending=[True, False], kind="mergesort")
    n_dup_rows = int(sub["_symbol"].duplicated().sum())
    sub = sub.drop_duplicates("_symbol", keep="first")
    symbols = sub["_symbol"].tolist()
    sub = sub.drop(columns=["_symbol", "_version"])
    sub.index = symbols
    sub.index.name = "gene_symbol"
    return sub, status, n_dup_rows


def value_stats(df, cancer, tissue, omics):
    arr = df.to_numpy(dtype=np.float64, copy=False)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        vmin = vmed = vmax = np.nan
    else:
        vmin = float(np.min(finite))
        vmed = float(np.median(finite))
        vmax = float(np.max(finite))
    with np.errstate(all="ignore"):
        sample_med = np.nanmedian(arr, axis=0)
    sample_med = sample_med[np.isfinite(sample_med)]
    def q(a, p):
        if a.size == 0:
            return np.nan
        return float(np.quantile(a, p))
    return {
        "cancer": cancer,
        "tissue": tissue,
        "omics": omics,
        "n_rows": int(df.shape[0]),
        "n_samples": int(df.shape[1]),
        "n_finite": int(finite.size),
        "min": vmin,
        "median": vmed,
        "max": vmax,
        "sample_median_min": q(sample_med, 0),
        "sample_median_q1": q(sample_med, 0.25),
        "sample_median_median": q(sample_med, 0.5),
        "sample_median_q3": q(sample_med, 0.75),
        "sample_median_max": q(sample_med, 1),
        "frac_missing": float(1 - finite.size / arr.size) if arr.size else np.nan,
    }


def matrix_path(cancer, tissue, omics):
    if omics == "rna":
        kind = f"RNAseq_gene_RSEM_coding_UQ_1500_log2_{tissue}"
    else:
        kind = f"proteomics_gene_abundance_log2_reference_intensity_normalized_{tissue}"
    return RAW / cancer / f"{cancer}_{kind}.txt"


def drop_duplicate_samples(df, dup):
    if not dup:
        return df
    drop_cols = []
    dupset = set(dup)
    for c in df.columns:
        base = str(c)
        if base in dupset:
            drop_cols.append(c)
            continue
        # pandas may rename duplicates as name.1
        stem = base.rsplit(".", 1)[0]
        if stem in dupset and base[len(stem) + 1 :].isdigit():
            drop_cols.append(c)
    return df.drop(columns=drop_cols)


def main():
    by_ensembl, protein_coding_symbols, _hgnc = load_hgnc()
    targets = pd.read_csv(PATHS["root"] / "config" / "adc_targets.tsv", sep="\t")
    target_genes = targets["gene_symbol"].tolist()

    dist_rows = []
    count_rows = []
    map_rows = []
    dup_rows = []
    rna_symbols = {}
    protein_rates = []
    adc_rows = []
    paired_rows = []

    for cancer in CANCERS:
        presence = {}
        collapsed = {}
        for tissue in ("Tumor", "Normal"):
            for omics in ("rna", "protein"):
                path = matrix_path(cancer, tissue, omics)
                if not path.exists():
                    presence[(tissue, omics)] = set()
                    dup_rows.append(
                        {"cancer": cancer, "tissue": tissue, "omics": omics, "n_duplicate_sample_ids": 0, "file": "ABSENT"}
                    )
                    continue
                samples, dup = read_header_samples(path)
                dup_rows.append(
                    {
                        "cancer": cancer,
                        "tissue": tissue,
                        "omics": omics,
                        "n_duplicate_sample_ids": len(dup),
                        "duplicate_ids": ",".join(dup),
                        "file": path.name,
                    }
                )
                df = read_matrix(path)
                df = drop_duplicate_samples(df, dup)
                # Exact string sample ids.
                df.columns = [str(c) for c in df.columns]
                dist_rows.append(value_stats(df, cancer, tissue, omics))
                sub, status, n_dup_sym = collapse(df, by_ensembl)
                for k, v in status.items():
                    map_rows.append(
                        {"cancer": cancer, "tissue": tissue, "omics": omics, "status": k, "n": v}
                    )
                map_rows.append(
                    {
                        "cancer": cancer,
                        "tissue": tissue,
                        "omics": omics,
                        "status": "collapsed_duplicate_symbol_rows",
                        "n": n_dup_sym,
                    }
                )
                out = PROC / f"{cancer}__{tissue.lower()}__{omics}.parquet"
                sub.astype(np.float32).to_parquet(out)
                presence[(tissue, omics)] = set(sub.columns)
                collapsed[(tissue, omics)] = sub
                print(f"{cancer} {tissue} {omics} {sub.shape}", flush=True)

        for tissue in ("Tumor", "Normal"):
            rna = presence.get((tissue, "rna"), set())
            prot = presence.get((tissue, "protein"), set())
            both = rna & prot
            count_rows.append(
                {
                    "cancer": cancer,
                    "tissue": tissue.lower(),
                    "n_rna": len(rna),
                    "n_protein": len(prot),
                    "n_rna_only": len(rna - prot),
                    "n_protein_only": len(prot - rna),
                    "n_both": len(both),
                }
            )
            for sid in sorted(both):
                paired_rows.append({"cancer": cancer, "tissue": tissue.lower(), "sample_id": sid})

        if ("Tumor", "rna") in collapsed:
            rna_symbols[cancer] = set(collapsed[("Tumor", "rna")].index)
        else:
            rna_symbols[cancer] = set()

        if ("Tumor", "protein") in collapsed:
            prot = collapsed[("Tumor", "protein")]
            n_s = prot.shape[1]
            measured = np.isfinite(prot.to_numpy(dtype=np.float32, copy=False)).sum(axis=1)
            for sym, n_m in zip(prot.index, measured):
                protein_rates.append(
                    {
                        "cancer": cancer,
                        "gene_symbol": sym,
                        "n_samples": n_s,
                        "n_measured": int(n_m),
                        "rate": float(n_m) / n_s if n_s else np.nan,
                    }
                )

        for tissue in ("Tumor", "Normal"):
            rna = collapsed.get((tissue, "rna"))
            prot = collapsed.get((tissue, "protein"))
            for gene in target_genes:
                row = {"cancer": cancer, "tissue": tissue.lower(), "gene_symbol": gene}
                if prot is not None and gene in prot.index:
                    v = prot.loc[gene].to_numpy(dtype=np.float64)
                    m = np.isfinite(v)
                    row["n_protein_samples"] = int(prot.shape[1])
                    row["n_protein_measured"] = int(m.sum())
                    row["protein_rate"] = float(m.mean()) if v.size else np.nan
                else:
                    row["n_protein_samples"] = 0 if prot is None else int(prot.shape[1])
                    row["n_protein_measured"] = 0 if prot is not None else np.nan
                    row["protein_rate"] = np.nan if prot is None or gene not in (prot.index if prot is not None else []) else np.nan
                    if prot is None:
                        row["n_protein_measured"] = np.nan
                        row["protein_rate"] = np.nan
                        row["n_protein_samples"] = np.nan
                if rna is not None and gene in rna.index:
                    v = rna.loc[gene].to_numpy(dtype=np.float64)
                    finite = v[np.isfinite(v)]
                    row["n_rna_samples"] = int(rna.shape[1])
                    row["rna_median"] = float(np.median(finite)) if finite.size else np.nan
                else:
                    row["n_rna_samples"] = np.nan if rna is None else int(rna.shape[1])
                    row["rna_median"] = np.nan
                adc_rows.append(row)

    # G_RNA and coverage summaries
    if rna_symbols:
        gset = set.intersection(*(rna_symbols[c] for c in CANCERS)) if all(rna_symbols[c] for c in CANCERS) else set()
    else:
        gset = set()
    gset &= protein_coding_symbols
    g_sorted = sorted(gset)
    (RES / "G_RNA_genes.txt").write_text("\n".join(g_sorted) + "\n")

    rates = pd.DataFrame(protein_rates)
    rates.to_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t", index=False)
    cov_rows = []
    for cancer in CANCERS:
        sub = rates[rates["cancer"] == cancer] if len(rates) else rates
        if sub.empty:
            cov_rows.append(
                {"cancer": cancer, "n_proteins": 0, "n_rate_100": 0, "n_rate_ge_70": 0, "n_rate_ge_30": 0}
            )
            continue
        cov_rows.append(
            {
                "cancer": cancer,
                "n_proteins": int(sub.shape[0]),
                "n_rate_100": int((sub["rate"] >= 1 - 1e-12).sum()),
                "n_rate_ge_70": int((sub["rate"] >= 0.70).sum()),
                "n_rate_ge_30": int((sub["rate"] >= 0.30).sum()),
            }
        )
    # P
    min_c = PLAN["p_protein"]["min_cancers_with_measurement_rate"]
    thr = PLAN["p_protein"]["measurement_rate_threshold"]
    if len(rates):
        ok = rates[rates["rate"] >= thr]
        n_c = ok.groupby("gene_symbol")["cancer"].nunique()
        p_genes = sorted(set(n_c[n_c >= min_c].index) & gset)
    else:
        p_genes = []
    (RES / "P_genes.txt").write_text("\n".join(p_genes) + "\n")

    pd.DataFrame(count_rows).to_csv(RES / "cptac_sample_counts.tsv", sep="\t", index=False)
    pd.DataFrame(dist_rows).to_csv(RES / "cptac_value_distribution.tsv", sep="\t", index=False)
    pd.DataFrame(cov_rows).to_csv(RES / "cptac_protein_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(map_rows).to_csv(RES / "cptac_gene_map_status.tsv", sep="\t", index=False)
    pd.DataFrame(dup_rows).to_csv(RES / "cptac_duplicate_samples.tsv", sep="\t", index=False)
    pd.DataFrame(adc_rows).to_csv(RES / "cptac_adc_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(paired_rows).to_csv(RES / "cptac_paired_samples.tsv", sep="\t", index=False)
    summary = pd.DataFrame(
        [
            {"item": "n_G_RNA", "value": len(g_sorted)},
            {"item": "n_P", "value": len(p_genes)},
            {"item": "p_min_cancers", "value": min_c},
            {"item": "p_rate_threshold", "value": thr},
        ]
    )
    summary.to_csv(RES / "cptac_set_sizes.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
