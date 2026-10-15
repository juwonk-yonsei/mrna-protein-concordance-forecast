"""Section 2. Antibody names and sample IDs only. No RPPA or ProCan values."""
import os

os.environ["B_CPU_CAP"] = "none"

import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.common import ROOT, sha256_file
from scripts.loader import structure_locked

PLAN = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
RES = ROOT / "results"
RAW = ROOT / "data" / "raw"
MAP = PLAN["rppa_map"]


def xena_key(text):
    cleaned = re.sub(r"[^A-Za-z0-9]", "", text).upper()
    if cleaned[:1].isdigit():
        cleaned = "X" + cleaned
    return cleaned


def name_flagged(text):
    return any(token in text for token in MAP["drop_if_name_contains"])


def load_cbioportal_names():
    files = sorted((RAW / "rppa_map").glob("*.antibody_col1.txt"))
    panels = []
    for path in files:
        lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
        if not lines or lines[0] != "Composite.Element.REF":
            raise RuntimeError(f"unexpected header in {path}")
        panels.append(tuple(lines[1:]))
    if not panels:
        raise RuntimeError("no cBioPortal antibody columns")
    if len(set(panels)) != 1:
        raise RuntimeError("cBioPortal RPPA name columns differ across studies")
    return files, list(panels[0])


def map_antibodies():
    files, names = load_cbioportal_names()
    xena = pd.read_csv(RES / "rppa_antibodies.tsv", sep="\t", dtype=str)["antibody"].tolist()
    by_key = defaultdict(list)
    for antibody in xena:
        by_key[xena_key(antibody)].append(antibody)
    seen = set(pd.read_csv(ROOT / "config" / "rppa_excluded_seen.tsv", sep="\t", dtype=str)["antibody"])
    g_rna = set(Path(RES / "G_RNA_genes.txt").read_text().split())
    targets = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t", dtype=str)["gene_symbol"].tolist()
    rows = []
    for raw in names:
        genes_field, token = raw.split("|", 1)
        symbols = [g for g in genes_field.split() if g and g != "-"]
        key = xena_key(token)
        hits = by_key.get(key, [])
        xena_name = hits[0] if len(hits) == 1 else "NA"
        reason = []
        if len(hits) != 1:
            reason.append("no_unique_xena_key")
        if len(symbols) != 1:
            reason.append("not_one_gene")
        if name_flagged(token) or (xena_name != "NA" and name_flagged(xena_name)):
            reason.append("phospho_or_cleaved")
        if xena_name in seen:
            reason.append("seen_in_stage1")
        gene = symbols[0] if len(symbols) == 1 else "NA"
        if gene != "NA" and gene not in g_rna:
            reason.append("not_in_g_rna")
        kept = not reason
        rows.append(
            {
                "cbioportal_name": raw,
                "antibody_token": token,
                "xena_antibody": xena_name,
                "gene_symbol": gene,
                "n_genes": len(symbols),
                "gene_symbols": " ".join(symbols) if symbols else "NA",
                "in_g_rna": gene in g_rna,
                "kept_for_forecast": kept,
                "exclude_reason": "ok" if kept else ";".join(reason),
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(RES / "rppa_antibody_map.tsv", sep="\t", index=False)
    kept = table[table["kept_for_forecast"]].copy()
    kept.to_csv(RES / "rppa_total_protein_map.tsv", sep="\t", index=False)
    adc_rows = []
    for gene in targets:
        sub = kept[kept["gene_symbol"] == gene]
        adc_rows.append(
            {
                "gene_symbol": gene,
                "n_antibodies": int(sub.shape[0]),
                "xena_antibodies": ";".join(sub["xena_antibody"]) if len(sub) else "NA",
            }
        )
    pd.DataFrame(adc_rows).to_csv(RES / "rppa_adc_total_antibodies.tsv", sep="\t", index=False)
    summary = pd.DataFrame(
        [
            {"item": "n_cbioportal_files", "value": len(files)},
            {"item": "n_cbioportal_antibodies", "value": len(names)},
            {"item": "n_unique_panels", "value": 1},
            {"item": "n_xena_antibodies", "value": len(xena)},
            {"item": "n_unique_xena_key_match", "value": int((table["xena_antibody"] != "NA").sum())},
            {"item": "n_kept_total_protein", "value": int(kept.shape[0])},
            {"item": "n_kept_in_g_rna", "value": int(kept["in_g_rna"].sum())},
            {"item": "n_adc_genes_with_antibody", "value": int(sum(r["n_antibodies"] > 0 for r in adc_rows))},
        ]
    )
    summary.to_csv(RES / "rppa_map_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)
    return kept


def rppa_rna_overlap():
    rppa = structure_locked("tcga_rppa")
    dup = set(pd.read_csv(RES / "rppa_duplicate_sample_ids.tsv", sep="\t", dtype=str)["sample_id"])
    rppa_ids = [s for s in rppa["col_names"] if s not in dup]
    # Drop the duplicated ID entirely, including a single remaining copy: the
    # file has two columns, and decision 11 drops both.
    rppa_ids = set(rppa_ids)
    rna_path = RAW / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"
    import gzip

    with gzip.open(rna_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    rna_ids = header[1:]
    from collections import Counter

    counts = Counter(rna_ids)
    rna_unique = {s for s, n in counts.items() if n == 1}
    ph = pd.read_csv(RAW / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", dtype=str)
    both = rppa_ids & rna_unique
    sub = ph[ph["sample"].isin(both)]
    grouped = sub.groupby("_primary_disease").size().rename("n_samples").reset_index()
    grouped["n_ge_30"] = grouped["n_samples"] >= 30
    grouped.to_csv(RES / "tcga_rppa_rna_overlap_by_disease.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"item": "n_rppa_ids_after_dropping_duplicates", "value": len(rppa_ids)},
            {"item": "n_duplicate_ids_dropped", "value": len(dup)},
            {"item": "n_rna_ids", "value": len(rna_ids)},
            {"item": "n_rna_ids_unique", "value": len(rna_unique)},
            {"item": "n_intersection", "value": len(both)},
            {"item": "n_intersection_in_phenotype", "value": int(sub["sample"].nunique())},
            {"item": "n_intersection_missing_phenotype", "value": int(len(both - set(ph["sample"])))},
            {"item": "n_diseases_ge_30", "value": int(grouped["n_ge_30"].sum())},
        ]
    ).to_csv(RES / "tcga_rppa_rna_overlap_summary.tsv", sep="\t", index=False)
    print("overlap diseases", int(grouped["n_ge_30"].sum()), "intersection", len(both), flush=True)


def procan_sanger_counts():
    mapping = pd.read_csv(
        RAW / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt", sep="\t", dtype=str
    )
    sample_ids = pd.read_csv(RES / "procan_sample_ids.tsv", sep="\t", dtype=str)["sample_id"]
    mapping = mapping[mapping["Project_Identifier"].isin(set(sample_ids))].copy()
    sanger = set()
    n_rows = 0
    path = RAW / "cmp" / "rnaseq_merged" / "rnaseq_merged_20260323.csv"
    with open(path, errors="replace") as handle:
        header = handle.readline().rstrip("\n").split(",")
        i_model = header.index("model_id")
        i_source = header.index("data_source")
        for line in handle:
            n_rows += 1
            parts = line.rstrip("\n").split(",")
            if len(parts) <= i_source:
                continue
            if parts[i_source] == "Sanger":
                sanger.add(parts[i_model])
    mapping["has_sanger"] = mapping["SIDM"].isin(sanger)
    use = mapping[mapping["has_sanger"]]
    counts = use.groupby("Tissue_type").size().rename("n_cell_lines").reset_index()
    counts["n_ge_30"] = counts["n_cell_lines"] >= 30
    counts.to_csv(RES / "procan_sanger_by_lineage.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"item": "n_procan_rows_in_matrix", "value": int(sample_ids.shape[0])},
            {"item": "n_mapping_rows_matched", "value": int(mapping.shape[0])},
            {"item": "n_sanger_model_ids", "value": len(sanger)},
            {"item": "n_procan_with_sanger", "value": int(use.shape[0])},
            {"item": "n_long_file_rows_scanned", "value": n_rows},
            {"item": "n_lineages_ge_30", "value": int(counts["n_ge_30"].sum())},
        ]
    ).to_csv(RES / "procan_sanger_summary.tsv", sep="\t", index=False)
    print(counts.to_string(index=False), flush=True)


def log_downloads():
    log = ROOT / "logs" / "download_log.tsv"
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = []
    for path in sorted((RAW / "rppa_map").glob("*.antibody_col1.txt")):
        study = path.name.split(".antibody_col1")[0]
        url = MAP["cbioportal_url"].format(study=study)
        lines.append(
            f"cbioportal\t{path.relative_to(ROOT)}\t{url}\t{path.stat().st_size}\t{sha256_file(path)}\t{stamp}\tok\n"
        )
    lines.append(
        "cbioportal\tdata/raw/rppa_map/laml_tcga_pan_can_atlas_2018.antibody_col1.txt\t"
        + MAP["cbioportal_url"].format(study="laml_tcga_pan_can_atlas_2018")
        + f"\tNA\tNA\t{stamp}\tfailed_http_404\n"
    )
    with log.open("a") as handle:
        handle.writelines(lines)


def main():
    log_downloads()
    map_antibodies()
    rppa_rna_overlap()
    procan_sanger_counts()
    print("maps done", flush=True)


if __name__ == "__main__":
    main()
