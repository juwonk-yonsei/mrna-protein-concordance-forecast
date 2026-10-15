"""Shared loader.

Locked objects are the ProCan protein-matrix values and the TCGA RPPA values.
Structure (names, counts) is available without --unlock. Numeric values raise
LockedDataError unless unlock=True. This stage never passes unlock=True.
"""
import os

# Stage 1 caps CPU at 16. Stage 2 sets B_CPU_CAP=none before import.
_cap = os.environ.get("B_CPU_CAP", "16")
if _cap != "none":
    os.environ["OMP_NUM_THREADS"] = _cap
    os.environ["OPENBLAS_NUM_THREADS"] = _cap
    os.environ["MKL_NUM_THREADS"] = _cap
    os.environ["NUMEXPR_MAX_THREADS"] = _cap

import gzip
from pathlib import Path

from scripts.common import ROOT, sha256_file

LOCKED = {
    "procan_protein": ROOT
    / "data/raw/procan/ProCan-DepMapSanger_protein_matrix_8498_averaged.txt",
    "tcga_rppa": ROOT / "data/raw/tcga/TCGA-RPPA-pancan-clean.xena.gz",
    "apollo_luad_protein": ROOT / "data/raw/apollo_luad/APOLLO1_Level3_Gprot_v061418.n87.csv",
}


class LockedDataError(PermissionError):
    pass


def locked_path(key):
    if key not in LOCKED:
        raise KeyError(key)
    return LOCKED[key]


def assert_unlocked(key, unlock):
    if key in LOCKED and not unlock:
        raise LockedDataError(
            f"{key} values are locked. Refusing to return numeric values without unlock."
        )


def iter_text_lines(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                yield line
    else:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                yield line


def structure_locked(key):
    """Names and dimensions only. Numeric fields are discarded as raw text and never parsed."""
    if key not in LOCKED:
        raise KeyError(key)
    path = LOCKED[key]
    if not path.exists():
        raise FileNotFoundError(path)
    if key == "procan_protein":
        return _structure_delimited(path)
    if key == "tcga_rppa":
        return _structure_xena(path)
    if key == "apollo_luad_protein":
        return _structure_apollo_csv(path)
    raise KeyError(key)


def _split_fields(line):
    line = line.rstrip("\n\r")
    if "\t" in line:
        return line.split("\t")
    return line.split()


def _structure_delimited(path):
    row_names = []
    col_names = None
    n_data = 0
    for i, line in enumerate(iter_text_lines(path)):
        if not line.strip() or line.startswith("#"):
            continue
        fields = _split_fields(line)
        if col_names is None:
            # Header: first cell is a row-name label, remaining cells are sample ids.
            col_names = fields[1:]
            continue
        row_names.append(fields[0])
        n_data += 1
    return {
        "key": "procan_protein",
        "path": str(path),
        "n_rows": n_data,
        "n_cols": 0 if col_names is None else len(col_names),
        "row_names": row_names,
        "col_names": col_names or [],
    }


def _structure_xena(path):
    header = None
    row_names = []
    for line in iter_text_lines(path):
        if not line.strip() or line.startswith("#"):
            continue
        fields = _split_fields(line)
        if header is None:
            header = fields
            continue
        row_names.append(fields[0])
        if header is not None and len(fields) == len(header):
            # First header cell is a label; samples start at field 1.
            sample_names = header[1:]
        else:
            sample_names = header
        # sample_names is set from the header once; values in fields[1:] are ignored.
        _ = sample_names
    # Recompute sample names from the stored header rule using the first data row's width.
    # The loop above keeps only the probe name. Sample ids come from the header.
    sample_names = _xena_samples(path)
    return {
        "key": "tcga_rppa",
        "path": str(path),
        "n_rows": len(row_names),
        "n_cols": len(sample_names),
        "row_names": row_names,
        "col_names": sample_names,
    }


def _xena_samples(path):
    header = None
    first_data = None
    for line in iter_text_lines(path):
        if not line.strip() or line.startswith("#"):
            continue
        fields = _split_fields(line)
        if header is None:
            header = fields
            continue
        first_data = fields
        break
    if header is None:
        return []
    if first_data is not None and len(first_data) == len(header):
        return header[1:]
    return header


def _structure_apollo_csv(path):
    """Header and gene symbols only. Numeric cells are not parsed."""
    import csv

    row_names = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        samples = header[4:]
        for fields in reader:
            if len(fields) < 4:
                continue
            row_names.append(fields[3])
    return {
        "key": "apollo_luad_protein",
        "path": str(path),
        "n_rows": len(row_names),
        "n_cols": len(samples),
        "row_names": row_names,
        "col_names": samples,
    }


def load_locked_matrix(key, unlock=False):
    """Numeric load. Callers must pass unlock=True."""
    assert_unlocked(key, unlock)
    import pandas as pd

    path = locked_path(key)
    if key == "apollo_luad_protein":
        frame = pd.read_csv(path)
        if "Gene" not in frame.columns:
            raise RuntimeError("APOLLO protein file has no Gene column")
        frame = frame.set_index("Gene")
        drop = [c for c in frame.columns if c in {"Accession", "Protein"} or str(c).startswith("Unnamed")]
        frame = frame.drop(columns=drop)
        frame.index = frame.index.astype(str)
        frame.columns = frame.columns.astype(str)
        return frame.apply(pd.to_numeric, errors="coerce")
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce")


def write_locked_hashes(dest):
    lines = ["key\tpath\tsha256\tsize_bytes\n"]
    for key, path in LOCKED.items():
        if not path.exists():
            lines.append(f"{key}\t{path}\tNA\tNA\n")
            continue
        digest = sha256_file(path)
        size = path.stat().st_size
        lines.append(f"{key}\t{path.relative_to(ROOT)}\t{digest}\t{size}\n")
    Path(dest).write_text("".join(lines))
    return dest
