import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import re
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.common import load_paths, load_plan
from scripts.loader import structure_locked, write_locked_hashes

PATHS = load_paths()
PLAN = load_plan()
RES = PATHS["results_dir"]
RAW = PATHS["raw_dir"]
PROC = PATHS["processed_dir"]
RES.mkdir(parents=True, exist_ok=True)
LOCUS = PLAN["gene_id"]["protein_coding_locus_group"]


def split_multi(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    text = str(val).strip()
    if not text or text == "nan":
        return []
    return [p.strip() for p in text.split("|") if p.strip()]


def hgnc_maps():
    h = pd.read_csv(RAW / "hgnc" / "hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False)
    pc = h[h["locus_group"] == LOCUS]
    approved = set(pc["symbol"].dropna())
    prev = {}
    alias = {}
    for rec in pc.to_dict("records"):
        sym = rec["symbol"]
        for p in split_multi(rec.get("prev_symbol")):
            prev.setdefault(p, set()).add(sym)
        for a in split_multi(rec.get("alias_symbol")):
            alias.setdefault(a, set()).add(sym)
    return approved, prev, alias


def map_symbol(token, approved, prev, alias):
    if token in approved:
        return token, "approved"
    if token in prev and len(prev[token]) == 1:
        return next(iter(prev[token])), "prev_symbol"
    if token in prev and len(prev[token]) > 1:
        return None, "ambiguous_prev"
    if token in alias and len(alias[token]) == 1:
        return next(iter(alias[token])), "alias"
    if token in alias and len(alias[token]) > 1:
        return None, "ambiguous_alias"
    return None, "unmapped"


def parse_depmap_col(col):
    text = str(col)
    if text.endswith(")") and " (" in text:
        return text.rsplit(" (", 1)[0]
    return text


def ccle_audit(approved, prev, alias):
    prot_path = RAW / "ccle" / "protein_quant_current_normalized.csv.gz"
    model_path = RAW / "depmap" / "Model.csv"
    rna_path = RAW / "depmap" / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    prot = pd.read_csv(prot_path, low_memory=False)
    meta_cols = [
        c
        for c in prot.columns
        if c
        in {
            "Protein_Id",
            "Gene_Symbol",
            "Description",
            "Group_ID",
            "Uniprot",
            "Uniprot_Acc",
            "TenPx",
        }
    ]
    if "Gene_Symbol" not in prot.columns:
        raise SystemExit(f"CCLE protein has no Gene_Symbol column: {list(prot.columns)[:30]}")
    sample_cols = [
        c
        for c in prot.columns
        if c not in meta_cols and not str(c).endswith("_Peptides")
    ]
    models = pd.read_csv(model_path, dtype=str, low_memory=False)
    need = ["ModelID", "OncotreeLineage", "CCLEName", "StrippedCellLineName"]
    missing = [c for c in need if c not in models.columns]
    if missing:
        (RES / "ccle_model_columns.txt").write_text("\n".join(models.columns) + "\n")
        raise SystemExit(f"Model.csv missing {missing}")
    lookup = {field: {} for field in PLAN["ccle"]["id_match_priority"]}
    for rec in models.to_dict("records"):
        mid = rec["ModelID"]
        for field in PLAN["ccle"]["id_match_priority"]:
            val = rec.get(field)
            if val and val != "nan":
                lookup[field].setdefault(val, mid)
    def match_key(col):
        # Column names are {CCLEName}_TenPx{plex}. Try the full string, then
        # the string with one trailing _TenPx<digits> removed.
        keys = [col]
        stripped = re.sub(r"_TenPx\d+$", "", col)
        if stripped != col:
            keys.append(stripped)
        return keys

    matched = []
    for col in sample_cols:
        hit = None
        how = None
        for key in match_key(col):
            for field in PLAN["ccle"]["id_match_priority"]:
                if key in lookup[field]:
                    hit = lookup[field][key]
                    how = field if key == col else field + "+strip_TenPx"
                    break
            if hit:
                break
        matched.append({"protein_column": col, "model_id": hit if hit else "NA", "matched_on": how if how else "NA"})
    match_df = pd.DataFrame(matched)
    # One model id must not be claimed by two protein columns.
    claimed = match_df[match_df["model_id"] != "NA"]
    dup_models = claimed[claimed["model_id"].duplicated(keep=False)]["model_id"].unique().tolist()
    if dup_models:
        match_df.loc[match_df["model_id"].isin(dup_models), "matched_on"] = "duplicate_model_excluded"
        match_df.loc[match_df["model_id"].isin(dup_models), "model_id"] = "NA"

    mapped = [map_symbol(str(s), approved, prev, alias) for s in prot["Gene_Symbol"].fillna("")]
    prot = prot.copy()
    prot["_symbol"] = [m[0] for m in mapped]
    prot["_map"] = [m[1] for m in mapped]
    map_counts = prot["_map"].value_counts().to_dict()
    use = prot[prot["_symbol"].notna()].copy()
    value_cols = sample_cols
    n_measured = np.isfinite(use[value_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)).sum(axis=1)
    use["_n"] = n_measured
    use = use.sort_values(["_symbol", "_n"], ascending=[True, False], kind="mergesort")
    n_collapsed = int(use["_symbol"].duplicated().sum())
    use = use.drop_duplicates("_symbol", keep="first")
    num = use[value_cols].apply(pd.to_numeric, errors="coerce")
    num.index = use["_symbol"].tolist()
    # Keep matched columns and rename to model id.
    keep_cols = match_df[match_df["model_id"] != "NA"]
    num = num[keep_cols["protein_column"].tolist()]
    num.columns = keep_cols["model_id"].tolist()
    num.index.name = "gene_symbol"
    outdir = PROC / "ccle"
    outdir.mkdir(parents=True, exist_ok=True)
    num.astype(np.float32).to_parquet(outdir / "protein.parquet")

    rna = pd.read_csv(rna_path, index_col=0)
    rna.index = rna.index.astype(str)
    gene_tokens = [parse_depmap_col(c) for c in rna.columns]
    gene_map = [map_symbol(t, approved, prev, alias) for t in gene_tokens]
    rna_status = pd.Series([m[1] for m in gene_map]).value_counts().to_dict()
    rna.columns = [m[0] if m[0] else f"UNMAPPED::{i}" for i, m in enumerate(gene_map)]
    rna = rna.loc[:, ~rna.columns.astype(str).str.startswith("UNMAPPED::")]
    rna = rna.T
    rna.index.name = "gene_symbol"
    n_measured_rna = np.isfinite(rna.to_numpy(dtype=np.float32)).sum(axis=1)
    rna = rna.assign(_n=n_measured_rna, _ord=np.arange(rna.shape[0])).reset_index()
    rna = rna.sort_values(["gene_symbol", "_n", "_ord"], ascending=[True, False, True], kind="mergesort")
    n_rna_dup = int(rna["gene_symbol"].duplicated().sum())
    rna = rna.drop_duplicates("gene_symbol", keep="first").drop(columns=["_n", "_ord"]).set_index("gene_symbol")
    # columns are model ids
    model_ids = [c for c in num.columns if c in set(rna.columns)]
    rna = rna[model_ids]
    num = num[model_ids]
    rna.astype(np.float32).to_parquet(outdir / "rna.parquet")
    meta = models[models["ModelID"].isin(model_ids)][
        ["ModelID", "OncotreeLineage", "CCLEName", "StrippedCellLineName"]
    ].drop_duplicates("ModelID")
    meta.to_csv(outdir / "models.tsv", sep="\t", index=False)

    lineage = meta["OncotreeLineage"].fillna("NA").value_counts()
    targets = pd.read_csv(PATHS["root"] / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    adc = []
    for gene in targets:
        if gene in num.index:
            v = num.loc[gene].to_numpy(dtype=np.float64)
            m = np.isfinite(v)
            adc.append(
                {
                    "gene_symbol": gene,
                    "n_cell_lines": int(v.size),
                    "n_measured": int(m.sum()),
                    "rate": float(m.mean()) if v.size else np.nan,
                }
            )
        else:
            adc.append({"gene_symbol": gene, "n_cell_lines": int(num.shape[1]), "n_measured": 0, "rate": np.nan})
    pd.DataFrame(adc).to_csv(RES / "ccle_adc_coverage.tsv", sep="\t", index=False)
    lineage.rename_axis("OncotreeLineage").reset_index(name="n_cell_lines").to_csv(
        RES / "ccle_lineage_counts.tsv", sep="\t", index=False
    )
    summary = pd.DataFrame(
        [
            {"item": "n_protein_rows_raw", "value": int(pd.read_csv(prot_path, usecols=["Gene_Symbol"]).shape[0])},
            {"item": "n_protein_sample_columns", "value": len(sample_cols)},
            {"item": "n_protein_columns_matched", "value": int((match_df["model_id"] != "NA").sum())},
            {"item": "n_duplicate_model_ids_excluded", "value": len(dup_models)},
            {"item": "n_protein_symbols", "value": int(num.shape[0])},
            {"item": "n_collapsed_duplicate_symbols", "value": n_collapsed},
            {"item": "n_paired_models", "value": len(model_ids)},
            {"item": "n_rna_models_in_file", "value": int(pd.read_csv(rna_path, usecols=[0]).shape[0])},
            {"item": "n_rna_genes_mapped", "value": int(rna.shape[0])},
            {"item": "n_rna_collapsed_duplicate_symbols", "value": n_rna_dup},
        ]
    )
    # The rna variable was transposed to genes x models, so shape[0] is genes.
    summary.to_csv(RES / "ccle_summary.tsv", sep="\t", index=False)
    pd.Series(map_counts).rename("n").rename_axis("protein_gene_map").reset_index().to_csv(
        RES / "ccle_protein_gene_map_status.tsv", sep="\t", index=False
    )
    pd.Series(rna_status).rename("n").rename_axis("rna_gene_map").reset_index().to_csv(
        RES / "ccle_rna_gene_map_status.tsv", sep="\t", index=False
    )
    match_df.to_csv(RES / "ccle_id_match.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)


def name_matches(antibody, symbol):
    if antibody == symbol:
        return "exact"
    if antibody.startswith(symbol + "_") or antibody.startswith(symbol + "-"):
        return "prefix"
    return ""


def locked_and_hpa(approved):
    write_locked_hashes(PATHS["root"] / "config" / "locked_hashes.tsv")
    targets = pd.read_csv(PATHS["root"] / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    pro = structure_locked("procan_protein")
    # The header label is Project_Identifier and the other header fields are
    # UniProt;SYMBOL_HUMAN. Rows are samples. Columns are proteins.
    def procan_symbol(name):
        if ";" in name and name.endswith("_HUMAN"):
            sym = name.split(";", 1)[1][: -len("_HUMAN")]
            return sym if sym else None
        return None

    protein_ids = pro["col_names"]
    sample_ids = pro["row_names"]
    parsed = [(pid, procan_symbol(pid)) for pid in protein_ids]
    symbols = [s for _, s in parsed if s]
    adc_rows = []
    for gene in targets:
        hits = [pid for pid, sym in parsed if sym == gene]
        adc_rows.append(
            {
                "gene_symbol": gene,
                "n_columns": len(hits),
                "column_ids": ";".join(hits) if hits else "NA",
            }
        )
    pd.DataFrame(adc_rows).to_csv(RES / "procan_adc_name_presence.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "n_samples": len(sample_ids),
                "n_protein_columns": len(protein_ids),
                "n_parsed_symbols": len(symbols),
                "n_unique_symbols": len(set(symbols)),
                "sample_id_example_1": sample_ids[0] if sample_ids else "NA",
                "protein_id_example_1": protein_ids[0] if protein_ids else "NA",
            }
        ]
    ).to_csv(RES / "procan_structure.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": sample_ids}).to_csv(RES / "procan_sample_ids.tsv", sep="\t", index=False)
    pd.DataFrame({"protein_id": [p for p, _ in parsed], "symbol": [s if s else "NA" for _, s in parsed]}).to_csv(
        RES / "procan_protein_ids.tsv", sep="\t", index=False
    )

    rppa = structure_locked("tcga_rppa")
    pd.DataFrame({"antibody": rppa["row_names"]}).to_csv(RES / "rppa_antibodies.tsv", sep="\t", index=False)
    hits = []
    for ab in rppa["row_names"]:
        for gene in targets:
            how = name_matches(ab, gene)
            if how:
                hits.append({"antibody": ab, "gene_symbol": gene, "name_rule": how})
    pd.DataFrame(hits).to_csv(RES / "rppa_adc_antibodies.tsv", sep="\t", index=False)
    ph_path = RAW / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
    ph = pd.read_csv(ph_path, sep="\t", dtype=str)
    ids = set(rppa["col_names"])
    sub = ph[ph["sample"].isin(ids)]
    counts = sub.groupby("_primary_disease").size().rename("n_samples").reset_index()
    counts.to_csv(RES / "rppa_samples_by_disease.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "n_antibodies": rppa["n_rows"],
                "n_samples": rppa["n_cols"],
                "n_samples_in_phenotype": int(sub.shape[0]),
                "n_samples_unmatched_phenotype": int(len(ids - set(ph["sample"]))),
                "antibody_example_1": rppa["row_names"][0] if rppa["row_names"] else "NA",
            }
        ]
    ).to_csv(RES / "rppa_structure.tsv", sep="\t", index=False)

    # CMP RNA structure: column names and whether a sequencing-source column exists.
    zpath = RAW / "cmp" / "rnaseq_merged_20260323.zip"
    zdest = RAW / "cmp" / "rnaseq_merged"
    if zpath.exists():
        with zipfile.ZipFile(zpath) as zf:
            names = zf.namelist()
            (RES / "cmp_rnaseq_zip_members.txt").write_text("\n".join(names) + "\n")
            zf.extractall(zdest)
    model = pd.read_csv(RAW / "cmp" / "model_list_20260921.csv", nrows=5, dtype=str)
    (RES / "cmp_model_list_columns.txt").write_text("\n".join(model.columns) + "\n")
    # Find a wide or annotation table and record headers only.
    headers = []
    for path in sorted(zdest.rglob("*")) if zdest.exists() else []:
        if path.suffix.lower() not in {".csv", ".tsv", ".txt"}:
            continue
        if path.stat().st_size > 5_000_000_000:
            continue
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        with open(path, errors="replace") as f:
            header = f.readline().rstrip("\n")
        headers.append(f"{path.name}\t{header[:2000]}")
    (RES / "cmp_rnaseq_headers.txt").write_text("\n".join(headers) + "\n")

    for label, zname in (
        ("hpa_normal", RAW / "hpa" / "normal_ihc_data.tsv.zip"),
        ("hpa_cancer", RAW / "hpa" / "cancer_data.tsv.zip"),
    ):
        if not zname.exists():
            continue
        with zipfile.ZipFile(zname) as zf:
            member = zf.namelist()[0]
            with zf.open(member) as f:
                header = f.readline().decode("utf-8", "replace").rstrip("\n")
                n = 1
                for _ in f:
                    n += 1
        (RES / f"{label}_header.txt").write_text(f"member\t{member}\nn_lines\t{n}\nheader\t{header}\n")
    print("aux structure written", flush=True)


def main():
    approved, prev, alias = hgnc_maps()
    ccle_audit(approved, prev, alias)
    locked_and_hpa(approved)


if __name__ == "__main__":
    main()
