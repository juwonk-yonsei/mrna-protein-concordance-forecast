"""APOLLO LUAD RNA processing and protein-header unit bounds.

Does not read APOLLO protein values. Does not import scripts/04_models.py or
scripts/02_cptac_audit.py. Does not read config/paths.yaml.
"""
import io
import os
import tarfile

os.environ["B_CPU_CAP"] = "none"

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.common import ROOT, sha256_file
from scripts.loader import structure_locked

RES = ROOT / "results"
RAW = ROOT / "data" / "raw" / "apollo_luad"
OUT = ROOT / "data" / "processed" / "apollo_luad"
PREP = RES / "stage6" / "prep"
TAR = RAW / "RNA_quant_files_20181226.tar.gz"


def read_quant(tar, member):
    handle = tar.extractfile(member)
    if handle is None:
        raise RuntimeError(f"missing {member.name}")
    text = io.TextIOWrapper(handle, encoding="utf-8")
    frame = pd.read_csv(text, sep="\t", usecols=["Gene", "coding", "TPM"])
    text.detach()
    return frame


def middle_id(protein_column):
    parts = str(protein_column).split(".")
    if len(parts) < 2:
        return None
    return f"{parts[0]}-{parts[1]}"


def apply_row_rules(row_names):
    multi = []
    single = []
    for name in row_names:
        text = "" if name is None else str(name)
        if any(sep in text for sep in (";", "|", ",")):
            multi.append(text)
        else:
            single.append(text)
    counts = Counter(single)
    duplicated = {symbol for symbol, n in counts.items() if n > 1}
    kept = [symbol for symbol in single if symbol not in duplicated]
    return {
        "n_rows": len(row_names),
        "n_multi_symbol_rows": len(multi),
        "n_duplicate_symbol_rows": int(sum(counts[s] for s in duplicated)),
        "n_duplicate_symbols": len(duplicated),
        "n_rows_kept": len(kept),
        "symbols": sorted(set(kept)),
    }


def main():
    print("thread env", {k: os.environ.get(k, "") for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}, flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    PREP.mkdir(parents=True, exist_ok=True)
    blocks = []
    n_files_with_dup = 0
    with tarfile.open(TAR, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".counts.quant.txt")]
        print("n_rna_files", len(members), flush=True)
        for member in members:
            frame = read_quant(tar, member)
            stem = Path(member.name).name.replace(".counts.quant.txt", "")
            if frame["Gene"].duplicated().any():
                n_files_with_dup += 1
            collapsed = frame.groupby("Gene", sort=False)["TPM"].sum(min_count=1)
            blocks.append(collapsed.rename(stem))
    # Re-read coding from the first file only for the format count.
    with tarfile.open(TAR, "r:gz") as tar:
        first = next(m for m in tar.getmembers() if m.name.endswith(".counts.quant.txt"))
        raw = read_quant(tar, first)
        first_header_line = tar.extractfile(first).readline().decode("utf-8").rstrip("\n")
    n_coding = int((raw["coding"] == 1).sum())
    n_noncoding = int((raw["coding"] == 0).sum())
    n_ensembl = int(raw["Gene"].astype(str).str.startswith("ENSG").sum())
    n_unique_gene = int(raw["Gene"].nunique())
    n_dup_gene_rows = int(raw["Gene"].duplicated().sum())

    tpm = pd.concat(blocks, axis=1, join="outer")
    tpm = tpm.apply(pd.to_numeric, errors="coerce")
    log2 = np.log2(tpm.to_numpy(dtype=np.float64) + 1.0)
    mat = pd.DataFrame(log2, index=tpm.index.astype(str), columns=tpm.columns.astype(str))
    mat.index.name = "gene_symbol"

    info = structure_locked("apollo_luad_protein")
    rules = apply_row_rules(info["row_names"])
    protein_cols = list(info["col_names"])
    rows = []
    for col in protein_cols:
        mid = middle_id(col)
        rna_name = f"{mid}.counts.quant.txt" if mid else ""
        rows.append(
            {
                "rna_file": f"RNA_quant_files/{rna_name}",
                "protein_column": col,
                "middle_id": mid if mid else "",
                "rna_in_matrix": int(mid in mat.columns) if mid else 0,
            }
        )
    sample_map = pd.DataFrame(rows)
    if sample_map["middle_id"].duplicated().any() or (sample_map["rna_in_matrix"] != 1).any():
        raise RuntimeError("protein columns do not match RNA files one to one")
    mat = mat.loc[:, sample_map["middle_id"].tolist()]
    parquet = OUT / "rna_log2tpm.parquet"
    mat.to_parquet(parquet)
    sample_path = ROOT / "config" / "apollo_sample_map.tsv"
    sample_map.drop(columns=["rna_in_matrix"]).to_csv(sample_path, sep="\t", index=False)

    g_rna = [line.strip() for line in (RES / "G_RNA_genes.txt").read_text().splitlines() if line.strip()]
    ar_set = set(line.strip() for line in (RES / "ar_set_genes.txt").read_text().splitlines() if line.strip())
    p_set = set(line.strip() for line in (RES / "P_genes.txt").read_text().splitlines() if line.strip())
    rna_genes = set(mat.index)
    apollo = set(rules["symbols"])
    absent = [g for g in g_rna if g not in rna_genes]
    u_p1 = sorted(ar_set & apollo)
    u_p2 = sorted(set(g_rna) & apollo)
    u_q = sorted(ar_set & p_set & apollo)

    qc_rows = []
    tpm_use = tpm.loc[mat.index, mat.columns]
    for sample in mat.columns:
        values = tpm_use[sample].to_numpy(dtype=np.float64)
        finite = values[np.isfinite(values)]
        qc_rows.append(
            {
                "middle_id": sample,
                "n_genes_finite": int(finite.size),
                "n_genes_tpm_gt_0": int(np.sum(finite > 0)),
                "tpm_sum": float(np.nansum(values)),
            }
        )
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(PREP / "rna_sample_qc.tsv", sep="\t", index=False)

    hashes = pd.DataFrame(
        [
            {"file": str(parquet.relative_to(ROOT)), "sha256": sha256_file(parquet), "size_bytes": parquet.stat().st_size},
            {"file": str(sample_path.relative_to(ROOT)), "sha256": sha256_file(sample_path), "size_bytes": sample_path.stat().st_size},
        ]
    )
    hashes.to_csv(PREP / "rna_hashes.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {"item": "n_rna_files", "value": sample_map.shape[0]},
            {"item": "header", "value": first_header_line},
            {"item": "n_gene_rows_one_file", "value": len(raw)},
            {"item": "n_unique_gene_ids_one_file", "value": n_unique_gene},
            {"item": "n_duplicate_gene_rows_one_file", "value": n_dup_gene_rows},
            {"item": "n_files_with_duplicate_genes", "value": n_files_with_dup},
            {"item": "n_coding_1", "value": n_coding},
            {"item": "n_coding_0", "value": n_noncoding},
            {"item": "n_ensembl_ids", "value": n_ensembl},
            {"item": "n_genes_in_matrix", "value": mat.shape[0]},
            {"item": "n_samples", "value": mat.shape[1]},
            {"item": "n_G_RNA", "value": len(g_rna)},
            {"item": "n_G_RNA_absent", "value": len(absent)},
            {"item": "n_G_RNA_present", "value": len(g_rna) - len(absent)},
            {"item": "protein_n_rows", "value": rules["n_rows"]},
            {"item": "protein_n_cols", "value": len(protein_cols)},
            {"item": "protein_multi_symbol_rows", "value": rules["n_multi_symbol_rows"]},
            {"item": "protein_duplicate_symbol_rows", "value": rules["n_duplicate_symbol_rows"]},
            {"item": "protein_duplicate_symbols", "value": rules["n_duplicate_symbols"]},
            {"item": "protein_rows_kept", "value": rules["n_rows_kept"]},
            {"item": "protein_symbols_kept", "value": len(rules["symbols"])},
            {"item": "U_P1", "value": len(u_p1)},
            {"item": "U_P2", "value": len(u_p2)},
            {"item": "U_Q", "value": len(u_q)},
            {"item": "tpm_sum_min", "value": float(qc["tpm_sum"].min())},
            {"item": "tpm_sum_max", "value": float(qc["tpm_sum"].max())},
            {"item": "n_detected_min", "value": int(qc["n_genes_tpm_gt_0"].min())},
            {"item": "n_detected_max", "value": int(qc["n_genes_tpm_gt_0"].max())},
        ]
    )
    summary.to_csv(PREP / "rna_summary.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": u_p1}).to_csv(PREP / "U_P1_genes.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": u_p2}).to_csv(PREP / "U_P2_genes.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_symbol": u_q}).to_csv(PREP / "U_Q_genes.tsv", sep="\t", index=False)
    dup_rows = []
    counts = Counter(s for s in info["row_names"] if not any(sep in str(s) for sep in (";", "|", ",")))
    for symbol, n in sorted(counts.items()):
        if n > 1:
            dup_rows.append({"gene_symbol": symbol, "n_rows": n})
    pd.DataFrame(dup_rows).to_csv(PREP / "protein_duplicate_symbols.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)
    print(hashes.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
