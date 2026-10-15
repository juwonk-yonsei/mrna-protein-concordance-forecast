"""Assemble release v1.1.0 and review pack v2. v1.0.0 is left in place.

No new LightGBM file is added: the reproduction grid did not meet 1e-12.
"""

import os

os.environ["B_CPU_CAP"] = "none"

import hashlib
import importlib.util
import shutil
import subprocess
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location(
    "stage10_release", Path(__file__).resolve().parents[1] / "scripts" / "39_stage10_release.py"
)
rel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rel)

ROOT = rel.ROOT
DIST = ROOT / "dist" / "B_release_v1.1.0"
PACK = ROOT / "manuscript" / "review_pack_srep_v2"
rel.DIST = DIST
rel.PACK = PACK


def write_readme(dist):
    history = subprocess.check_output(
        ["git", "log", "--format=%H %ci %s", "--",
         "config/prereg_B4.md", "config/prereg_B6.md",
         "config/analysis_plan_B3.yaml", "config/analysis_plan_B4.yaml",
         "config/analysis_plan_B6.yaml", "config/analysis_plan_B10.yaml",
         "config/analysis_plan_B11.yaml", "config/analysis_plan_B12.yaml"],
        cwd=ROOT, text=True,
    ).strip()
    scripts = []
    for path in sorted((ROOT / "scripts").glob("*.py")):
        text = path.read_text(errors="replace")
        first = ""
        if text.startswith('"""'):
            first = text.split('"""', 2)[1].strip().splitlines()[0]
        scripts.append(f"- `scripts/{path.name}`: {first}")
    lines = [
        "# Project B release v1.1.0",
        "",
        "Local copy of the code, analysis plans, frozen outputs and Scientific Reports files for the mRNA–protein concordance forecast. v1.0.0 is unchanged and is not in this folder.",
        "",
        "Raw matrices are not in this folder. The manuscript Data availability section lists the public locations. This folder does not redistribute them.",
        "",
        "Also absent: sample identifiers, the planted-signal and null checks, and any file over 100 MB. The three large files are listed in `zenodo_only/files.tsv`.",
        "",
        "No license file is included. The author has not chosen one.",
        "Author names, the AI model name, funding, a GitHub URL and a Zenodo DOI are not filled in.",
        "",
        "## Forecast table v3",
        "",
        "`results/stage12/forecast_all_genes_v3.tsv.gz` adds `model_gene_features` to the v2 table. The other columns are the v2 columns.",
        "",
        "If model_gene_features is 0 and cptac_pooled_rho is available, the pooled CPTAC correlation was more accurate than the forecast in a post hoc comparison on APOLLO and should be preferred.",
        "",
        "## Forecast model",
        "",
        "One model file is released, `release/models/F1s_lightgbm.txt`. Reloading it reproduces the frozen APOLLO predictions. A post hoc grid over 1, 16 and 72 threads and over automatic, column-wise and row-wise histograms did not reproduce the frozen RPPA or ProCan predictions within 1e-12, so no second model file was written.",
        "",
        "## Commands",
        "",
        "Set `B_CPU_CAP=none` before a project script. Redraw the figures from `results/figures/source` and the stored result tables with:",
        "",
        "```",
        "B_CPU_CAP=none python scripts/36_stage10_figures.py",
        "```",
        "",
        "Scripts, in filename order:",
        "",
        *scripts,
        "",
        "`scripts_freeze/08_stage3.py` is the file from commit `224779714bffd58b02059f51de48d92d1324a24f`.",
        "",
        "## Plan commits",
        "",
        "```",
        history,
        "```",
        "",
        "Result tables keep their original column headers. `excluded.tsv` lists files that were not copied.",
        "",
    ]
    (dist / "README.md").write_text("\n".join(lines))


def before_after():
    rows = [
        ["R1 abstract", "at least as accurate as looking up observed CPTAC correlations, even from the matched cancer type", "at least as accurate as looking up observed CPTAC correlations, pooled or cancer-matched, for genes whose CPTAC data it had been given"],
        ["R2 results", "two paragraphs on the pooled lookup and genes outside P", "L1 pooled-lookup paragraph and the 600/6382 APOLLO split"],
        ["R3 ADC count", "Prescribing information opened for this revision confirmed eight approved ADCs acting on six of the 29 targets", "Eight approved ADCs act on six of the 29 targets"],
        ["R4 discussion", "the forecast was at least as accurate as the lookup even when the lookup came from the matched cancer type", "where the forecast had been given a gene's CPTAC data, it was at least as accurate as these lookups; where a lookup used correlations the features had not used, the lookup was better"],
        ["R5 methods", "interval half-width sentence ended the paragraph", "OPTION N: the released model reproduces APOLLO and not RPPA or ProCan, and the cause was not identified"],
        ["R6 pooled lookup", "mean of all finite CPTAC correlations, computed also where a matched correlation existed", "mean across CPTAC cancer types, including the matched cancer type"],
        ["R7 targets", "The 29 targets are the gene symbols listed in the project target table", "listed in the first analysis plan, before any data were analysed"],
        ["R8 DepMap", "DOI already ended in .v1", "unchanged"],
        ["R9 Table 3", "six rows", "ten rows, including the covered and not-covered APOLLO split"],
        ["R10 cover", "even from the matched cancer type", "pooled or matched, for genes whose CPTAC data it had been given; otherwise a lookup was better"],
        ["R11 note 2", "refit paragraph without the thread grid", "OPTION N, including refits at 1, 16 and 72 threads"],
        ["R12 note 3", "plans through the pooled lookup", "adds the B11 split and the B12 reproduction plan and output commits"],
        ["R13 note 6", "On the matched-cancer subset the counts were 600, 0 and 0", "600 APOLLO L1 pairs had default gene-level features; no RPPA or ProCan L1 pairs did"],
        ["R14 note 7 and S7", "labels opened on 14 October 2026; 23 targets not confirmed", "openFDA search on 15 October 2026; product and query table"],
    ]
    lines = ["| Block | Before | After |", "|---|---|---|"]
    for block, before, after in rows:
        lines.append(f"| {block} | {before} | {after} |")
    return "\n".join(lines) + "\n"


def review_pack():
    if PACK.exists():
        shutil.rmtree(PACK)
    PACK.mkdir(parents=True)
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep.md", PACK / "B_main_srep.md")
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep_claude.md", PACK / "B_main_srep_annotated.md")
    shutil.copy2(ROOT / "manuscript/srep/B_SI_srep.md", PACK / "B_SI_srep.md")
    shutil.copy2(ROOT / "manuscript/claude/B_revision_B12_claude.md", PACK / "B_revision_B12_claude.md")
    shutil.copy2(ROOT / "manuscript/srep/B_cover_letter_srep.md", PACK / "B_cover_letter_srep.md")
    for name in [
        "submission_srep/B_manuscript_srep.docx",
        "submission_srep/B_SI_srep.pdf",
        "submission_srep/B_cover_letter_srep.docx",
        "submission_srep/checklist_srep.md",
    ]:
        shutil.copy2(ROOT / name, PACK / Path(name).name)
    main_text = (ROOT / "manuscript/srep/B_main_srep.md").read_text()
    legends = main_text.split("## Figure legends", 1)[1].split("## Tables", 1)[0]
    tables = main_text.split("## Tables", 1)[1]
    (PACK / "figure_legends.md").write_text("## Figure legends\n" + legends)
    (PACK / "tables.md").write_text("## Tables\n" + tables)
    (PACK / "revision_before_after.md").write_text(before_after())
    ledger = pd.read_csv(ROOT / "manuscript/checks/ledger_srep.tsv", sep="\t")
    summary = ledger.groupby(["유형", "지위"]).size().reset_index(name="n")
    summary.to_csv(PACK / "ledger_summary.tsv", sep="\t", index=False)
    remaining = []
    for path in (ROOT / "manuscript/srep").glob("B_*srep.md"):
        if path.name.endswith("_claude.md"):
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if "AUTHOR TO FILL" in line or "AUTHOR TO CONFIRM" in line:
                remaining.append(f"{path.name}:{i}: {line.strip()[:240]}")
    (PACK / "remaining_placeholders.txt").write_text("\n".join(remaining) + "\n")
    pandoc = "/tmp/pandoc-3.6.4/bin/pandoc"
    soffice = "/tmp/lo/root/opt/libreoffice26.8/program/soffice"
    work = PACK / "_work"
    work.mkdir()
    subprocess.check_call([pandoc, str(PACK / "figure_legends.md"), "-o", str(work / "legends.docx")])
    subprocess.check_call([
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b12pack",
        "--convert-to", "pdf", "--outdir", str(work), str(work / "legends.docx"),
    ])
    subprocess.check_call([
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b12pack",
        "--convert-to", "pdf", "--outdir", str(work),
        str(ROOT / "submission_srep/B_manuscript_srep.docx"),
    ])
    figure_pdfs = [
        ROOT / "submission_srep/figures" / name
        for name in ["Fig1.pdf", "Fig2.pdf", "Fig3.pdf", "Fig4.pdf", "Fig5.pdf", "FigS1.pdf", "FigS2.pdf"]
    ]
    subprocess.check_call([
        "pdfunite", *[str(path) for path in figure_pdfs], str(work / "legends.pdf"),
        str(PACK / "figures_and_legends.pdf"),
    ])
    shutil.copy2(work / "B_manuscript_srep.pdf", PACK / "B_manuscript_srep.pdf")
    shutil.rmtree(work)
    fig_dir = PACK / "figures"
    fig_dir.mkdir()
    for path in (ROOT / "submission_srep" / "figures").glob("*.png"):
        shutil.copy2(path, fig_dir / path.name)


def main():
    if (ROOT / "dist" / "B_release_v1.0.0").exists() is False:
        raise SystemExit("v1.0.0 is missing")
    if DIST.exists():
        shutil.rmtree(DIST)
    DIST.mkdir(parents=True)
    copied, skipped = [], []
    for rel_path in [
        "scripts", "config",
        "results/stage4/confirm", "results/stage5/posthoc", "results/stage5/atlas",
        "results/stage6/confirm", "results/stage6/posthoc", "results/stage6/prep", "results/stage6/frozen",
        "results/stage7", "results/stage8", "results/stage9", "results/stage10",
        "results/stage11", "results/stage12",
        "results/figures/source",
        "results/f1s_predictions.tsv", "results/reliability_predictions.tsv", "results/f1s_model_fit.tsv",
        "results/cptac_sample_counts.tsv", "results/cptac_set_sizes.tsv",
        "results/rppa_structure.tsv", "results/procan_sanger_by_lineage.tsv",
        "release/models/F1s_lightgbm.txt", "release/models/F1s_lightgbm_features.json",
        "release/model_cards.md",
        "manuscript/srep", "manuscript/claude", "manuscript/checks",
        "submission_srep/figures",
        "submission_srep/B_manuscript_srep.docx", "submission_srep/B_SI_srep.pdf",
        "submission_srep/B_cover_letter_srep.docx", "submission_srep/checklist_srep.md",
    ]:
        rel.copy_tree(rel_path, copied, skipped)
    freeze_dir = DIST / "scripts_freeze"
    freeze_dir.mkdir()
    text = subprocess.check_output(["git", "show", f"{rel.FREEZE}:scripts/08_stage3.py"], cwd=ROOT)
    (freeze_dir / "08_stage3.py").write_bytes(text)
    digest = hashlib.sha256(text).hexdigest()
    if digest != "b14c554c86f2243c52c9dba339ae0205bcf773f001a69e964c40bfbe90fea747":
        raise SystemExit(f"freeze hash {digest}")
    (DIST / "requirements.txt").write_text(
        "numpy==1.26.4\npandas==2.2.2\nscipy==1.13.1\nscikit-learn==1.4.2\nlightgbm==4.7.0\nmatplotlib==3.8.4\n"
    )
    large = [
        ("release/models/r1r2_model_compressed.npz", 121062358, "4c1b0f6f7a75cbb91b1043f01005d154e01ec2ba74d28da0402918beccf6d4cc"),
        ("results/stage6/frozen/r1r2_model.npz", 121062358, "4c1b0f6f7a75cbb91b1043f01005d154e01ec2ba74d28da0402918beccf6d4cc"),
        ("results/reliability_training_rows.tsv", (ROOT / "results/reliability_training_rows.tsv").stat().st_size, rel.sha256(ROOT / "results/reliability_training_rows.tsv")),
    ]
    zdir = DIST / "zenodo_only"
    zdir.mkdir()
    (zdir / "files.tsv").write_text(
        "path\tbytes\tsha256\treason\n" + "".join(f"{path}\t{size}\t{digest}\tover 100 MB; not in this folder\n" for path, size, digest in large)
    )
    write_readme(DIST)
    spec = importlib.util.spec_from_file_location("figures", ROOT / "scripts" / "36_stage10_figures.py")
    figures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(figures)
    smoke = DIST / "_smoke"
    smoke.mkdir()
    figures.OUT = smoke
    for name in ("fig1", "fig2", "fig3", "fig4", "fig5", "fig_s1", "fig_s2"):
        getattr(figures, name)()
    smoke_lines = ["figure\tbytes"]
    for path in sorted(smoke.glob("Fig*.pdf")):
        if path.stat().st_size < 1000:
            raise SystemExit(f"smoke figure too small: {path}")
        smoke_lines.append(f"{path.name}\t{path.stat().st_size}")
    (DIST / "smoke_redraw.tsv").write_text("\n".join(smoke_lines) + "\n")
    shutil.rmtree(smoke)
    (DIST / "excluded.tsv").write_text(
        "path\tbytes\treason\n" + "".join(f"{path}\t{size}\t{reason}\n" for path, size, reason in sorted(skipped))
    )
    rows = ["path\tbytes\tsha256"]
    for path in sorted(p for p in DIST.rglob("*") if p.is_file() and p.name != "MANIFEST.tsv"):
        rows.append(f"{path.relative_to(DIST).as_posix()}\t{path.stat().st_size}\t{rel.sha256(path)}")
    (DIST / "MANIFEST.tsv").write_text("\n".join(rows) + "\n")
    review_pack()
    pack_bytes = sum(path.stat().st_size for path in PACK.rglob("*") if path.is_file())
    if pack_bytes > 25 * 1024 * 1024:
        raise SystemExit(f"review pack {pack_bytes} bytes")
    (DIST / ".gitignore").write_text("*\n!*/\n!.gitignore\n!README.md\n!MANIFEST.tsv\n!requirements.txt\n!excluded.tsv\n!smoke_redraw.tsv\n!zenodo_only/**\n!scripts_freeze/**\n")
    sentence = "If model_gene_features is 0 and cptac_pooled_rho is available, the pooled CPTAC correlation was more accurate than the forecast in a post hoc comparison on APOLLO and should be preferred."
    if sentence not in (DIST / "README.md").read_text():
        raise SystemExit("README sentence missing")
    if (DIST / "release/models/F1s_lightgbm_rppa_procan.txt").exists():
        raise SystemExit("OPTION N must not add a second model")
    print(f"copied {len(copied)} skipped {len(skipped)} pack {pack_bytes}", flush=True)


if __name__ == "__main__":
    main()
