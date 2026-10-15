import os

os.environ["OMP_NUM_THREADS"] = "16"
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["MKL_NUM_THREADS"] = "16"
os.environ["NUMEXPR_MAX_THREADS"] = "16"

import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.common import free_gb, load_paths, sha256_file

PATHS = load_paths()
BASE = "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized"
# Files listed on https://kb.linkedomics.org/download for data freeze v1.2.
# Phosphoproteome, mutation, methylation, and isoform files are omitted.
CPTAC_REL = """
BRCA/BRCA_meta.txt
BRCA/BRCA_Tumor_CaseList.txt
BRCA/BRCA_phenotype.txt
BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
BRCA/BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
CCRCC/CCRCC_meta.txt
CCRCC/CCRCC_Tumor_CaseList.txt
CCRCC/CCRCC_Normal_CaseList.txt
CCRCC/CCRCC_phenotype.txt
CCRCC/CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
CCRCC/CCRCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt
CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
CCRCC/CCRCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
COAD/COAD_meta.txt
COAD/COAD_Tumor_CaseList.txt
COAD/COAD_Normal_CaseList.txt
COAD/COAD_phenotype.txt
COAD/COAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
COAD/COAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
GBM/GBM_meta.txt
GBM/GBM_Tumor_CaseList.txt
GBM/GBM_phenotype.txt
GBM/GBM_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
GBM/GBM_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
HNSCC/HNSCC_meta.txt
HNSCC/HNSCC_Tumor_CaseList.txt
HNSCC/HNSCC_Normal_CaseList.txt
HNSCC/HNSCC_phenotype.txt
HNSCC/HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
HNSCC/HNSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt
HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
HNSCC/HNSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
LSCC/LSCC_meta.txt
LSCC/LSCC_Tumor_CaseList.txt
LSCC/LSCC_Normal_CaseList.txt
LSCC/LSCC_phenotype.txt
LSCC/LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
LSCC/LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt
LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
LUAD/LUAD_meta.txt
LUAD/LUAD_Tumor_CaseList.txt
LUAD/LUAD_Normal_CaseList.txt
LUAD/LUAD_phenotype.txt
LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
LUAD/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt
LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
OV/OV_meta.txt
OV/OV_Tumor_CaseList.txt
OV/OV_Normal_CaseList.txt
OV/OV_phenotype.txt
OV/OV_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
OV/OV_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
PDAC/PDAC_meta.txt
PDAC/PDAC_Tumor_CaseList.txt
PDAC/PDAC_Normal_CaseList.txt
PDAC/PDAC_phenotype.txt
PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
PDAC/PDAC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt
PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
PDAC/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
UCEC/UCEC_meta.txt
UCEC/UCEC_Tumor_CaseList.txt
UCEC/UCEC_Normal_CaseList.txt
UCEC/UCEC_phenotype.txt
UCEC/UCEC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
UCEC/UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt
""".strip().split()


def stamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log_row(source, file, url, size, digest, status):
    path = PATHS["logs_dir"] / "download_log.tsv"
    if not path.exists():
        path.write_text("source\tfile\turl\tsize_bytes\tsha256\tdatetime\tstatus\n")
    with path.open("a") as f:
        f.write(f"{source}\t{file}\t{url}\t{size}\t{digest}\t{stamp()}\t{status}\n")


def disk_ok(expected_bytes):
    free = free_gb()
    if expected_bytes and expected_bytes > 1_000_000_000:
        note = PATHS["logs_dir"] / "disk_before_large.txt"
        with note.open("a") as f:
            f.write(f"{stamp()}\tfree_gb={free:.1f}\texpected_bytes={expected_bytes}\n")
    remain = free - (expected_bytes or 0) / 1e9
    if remain < 100:
        return False
    return True


def content_length(url):
    try:
        out = subprocess.run(
            ["curl", "-sI", "-L", "-A", "Mozilla/5.0", "--max-time", "40", url],
            check=False,
            capture_output=True,
            text=True,
        ).stdout
    except OSError:
        return None
    length = None
    for line in out.splitlines():
        if line.lower().startswith("content-length:"):
            try:
                length = int(line.split(":", 1)[1].strip())
            except ValueError:
                length = None
    return length


def download_one(source, url, dest: Path, expected=None):
    dest.parent.mkdir(parents=True, exist_ok=True)
    rel = str(dest.relative_to(PATHS["root"]))
    if dest.exists() and dest.stat().st_size > 0:
        digest = sha256_file(dest)
        log_row(source, rel, url, dest.stat().st_size, digest, "exists")
        return True
    size_hint = expected if expected else content_length(url)
    if not disk_ok(size_hint or 0):
        log_row(source, rel, url, "NA", "NA", "skipped_disk")
        print(f"DISK skip {rel}", flush=True)
        return False
    last = "failed"
    tmp = dest.with_suffix(dest.suffix + ".partial")
    for attempt in range(1, 4):
        try:
            subprocess.run(
                [
                    "curl",
                    "-fL",
                    "-C",
                    "-",
                    "--retry",
                    "0",
                    "-A",
                    "Mozilla/5.0",
                    "--max-time",
                    "21600",
                    "-o",
                    str(tmp),
                    url,
                ],
                check=True,
            )
            if tmp.stat().st_size == 0:
                raise RuntimeError("empty file")
            tmp.replace(dest)
            digest = sha256_file(dest)
            log_row(source, rel, url, dest.stat().st_size, digest, "ok")
            print(f"OK {rel} {dest.stat().st_size}", flush=True)
            return True
        except Exception as exc:
            last = f"failed_attempt{attempt}:{type(exc).__name__}"
            print(f"FAIL {rel} attempt {attempt} {exc}", flush=True)
            if tmp.exists():
                tmp.unlink()
    log_row(source, rel, url, "NA", "NA", last)
    return False


def main():
    raw = PATHS["raw_dir"]
    jobs = []
    for rel in CPTAC_REL:
        jobs.append(("cptac", f"{BASE}/{rel}", raw / "cptac" / rel, None))

    jobs.extend(
        [
            (
                "hgnc",
                "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                raw / "hgnc" / "hgnc_complete_set.txt",
                16953545,
            ),
            (
                "ccle_proteome",
                "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
                raw / "ccle" / "protein_quant_current_normalized.csv.gz",
                30916232,
            ),
            (
                "depmap",
                "https://depmap.org/portal/api/no-captcha/download/files",
                raw / "depmap" / "download_files_nocaptcha.csv",
                None,
            ),
            (
                "depmap",
                "https://ndownloader.figshare.com/files/51065297",
                raw / "depmap" / "Model.csv",
                None,
            ),
            (
                "depmap",
                "https://ndownloader.figshare.com/files/51065489",
                raw / "depmap" / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                None,
            ),
            (
                "depmap",
                "https://ndownloader.figshare.com/files/51065795",
                raw / "depmap" / "README_24Q4.txt",
                None,
            ),
            (
                "procan",
                "https://ndownloader.figshare.com/files/34411172",
                raw / "procan" / "ProCan-DepMapSanger_protein_matrix_8498_averaged.txt",
                94222540,
            ),
            (
                "procan",
                "https://ndownloader.figshare.com/files/34411133",
                raw / "procan" / "ProCan-DepMapSanger_mapping_file_averaged.txt",
                82593,
            ),
            (
                "procan",
                "https://ndownloader.figshare.com/files/36155439",
                raw / "procan" / "figshare_README.md",
                3218,
            ),
            (
                "cmp",
                "https://cog.sanger.ac.uk/cmp/download/model_list_20260921.csv",
                raw / "cmp" / "model_list_20260921.csv",
                957555,
            ),
            (
                "cmp",
                "https://cog.sanger.ac.uk/cmp/download/rnaseq_merged_20260323.zip",
                raw / "cmp" / "rnaseq_merged_20260323.zip",
                1291769946,
            ),
            (
                "tcga",
                "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA-RPPA-pancan-clean.xena.gz",
                raw / "tcga" / "TCGA-RPPA-pancan-clean.xena.gz",
                5602638,
            ),
            (
                "tcga",
                "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
                raw / "tcga" / "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
                331000731,
            ),
            (
                "tcga",
                "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
                raw / "tcga" / "TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
                61165,
            ),
            (
                "hpa",
                "https://www.proteinatlas.org/download/normal_tissue.tsv.zip",
                raw / "hpa" / "normal_tissue.tsv.zip",
                None,
            ),
            (
                "hpa",
                "https://www.proteinatlas.org/download/pathology.tsv.zip",
                raw / "hpa" / "pathology.tsv.zip",
                None,
            ),
            (
                "hpa",
                "https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip",
                raw / "hpa" / "normal_ihc_data.tsv.zip",
                5732831,
            ),
            (
                "hpa",
                "https://www.proteinatlas.org/download/tsv/cancer_data.tsv.zip",
                raw / "hpa" / "cancer_data.tsv.zip",
                1724018,
            ),
            (
                "sanger_docs",
                "https://depmap.sanger.ac.uk/documentation/datasets/expression/",
                raw / "docs" / "sanger_expression.html",
                None,
            ),
            (
                "sanger_docs",
                "https://depmap.sanger.ac.uk/documentation/datasets/proteomics/",
                raw / "docs" / "sanger_proteomics.html",
                None,
            ),
        ]
    )
    # Non-critical probe-map attempts. Failure is logged and ignored.
    for url, name in [
        (
            "https://pancanatlas.xenahubs.net/download/probeMap/TCGA-RPPA-pancan-clean.probeMap",
            "TCGA-RPPA-pancan-clean.probeMap",
        ),
        (
            "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/probeMap/TCGA-RPPA-pancan-clean.probeMap",
            "TCGA-RPPA-pancan-clean.probeMap.s3",
        ),
    ]:
        jobs.append(("tcga_probemap", url, raw / "tcga" / name, None))

    ok = 0
    for source, url, dest, expected in jobs:
        if download_one(source, url, dest, expected):
            ok += 1
    print(f"done {ok}/{len(jobs)}", flush=True)


if __name__ == "__main__":
    main()
