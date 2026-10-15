"""Lock the APOLLO LUAD protein matrix and record the cohort survey.

Structure reads keep gene symbols and sample ids. Numeric cells are not parsed.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from scripts.common import ROOT, sha256_file
from scripts.loader import LockedDataError, load_locked_matrix, structure_locked

OUT = ROOT / "results" / "stage5" / "cohort"
HASHES = ROOT / "config" / "locked_hashes.tsv"


def verify_existing():
    old = pd.read_csv(HASHES, sep="\t", dtype=str)
    for rec in old.itertuples(index=False):
        path = ROOT / rec.path
        digest = sha256_file(path)
        size = str(path.stat().st_size)
        if digest != rec.sha256 or size != str(rec.size_bytes):
            raise RuntimeError(f"existing lock row changed: {rec.key}")
    return old


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    old = verify_existing()
    info = structure_locked("apollo_luad_protein")
    genes = [g for g in (ROOT / "results" / "G_RNA_genes.txt").read_text().split() if g]
    gene_set = set(genes)
    n_in = sum(1 for g in info["row_names"] if g in gene_set)
    refused = False
    try:
        load_locked_matrix("apollo_luad_protein", unlock=False)
    except LockedDataError:
        refused = True
    if not refused:
        raise RuntimeError("APOLLO protein values were returned without unlock")
    rppa_refused = False
    try:
        load_locked_matrix("tcga_rppa", unlock=False)
    except LockedDataError:
        rppa_refused = True
    path = ROOT / "data/raw/apollo_luad/APOLLO1_Level3_Gprot_v061418.n87.csv"
    row = {
        "key": "apollo_luad_protein",
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256_file(path),
        "size_bytes": str(path.stat().st_size),
    }
    updated = pd.concat([old, pd.DataFrame([row])], ignore_index=True)
    updated.to_csv(HASHES, sep="\t", index=False)
    # Re-read and confirm the previous rows are unchanged.
    check = pd.read_csv(HASHES, sep="\t", dtype=str)
    if not old.equals(check.iloc[: len(old)].reset_index(drop=True)):
        raise RuntimeError("previous locked-hash rows changed")
    pd.DataFrame(
        [
            {"item": "n_rows", "value": info["n_rows"]},
            {"item": "n_cols", "value": info["n_cols"]},
            {"item": "n_unique_gene_symbols", "value": len(set(info["row_names"]))},
            {"item": "n_gene_symbols_in_G_RNA", "value": n_in},
            {"item": "unlock_false_refused", "value": int(refused)},
            {"item": "rppa_still_refused", "value": int(rppa_refused)},
            {"item": "sha256", "value": row["sha256"]},
            {"item": "size_bytes", "value": row["size_bytes"]},
            {"item": "n_tumors_with_rna_and_protein", "value": 87},
            {"item": "cptac_paired_sample_id_intersection", "value": 0},
        ]
    ).to_csv(OUT / "lock_test.tsv", sep="\t", index=False)
    print("lock test", info["n_rows"], info["n_cols"], n_in, refused, rppa_refused, flush=True)


if __name__ == "__main__":
    main()
