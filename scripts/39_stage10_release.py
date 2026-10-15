"""Assemble the local release folder and the review pack.

The folder is a copy of code, plans, frozen outputs and the manuscript.
Raw matrices, sample identifiers, the planted and null checks, and files
over 100 MB stay out. The 121 MB anchor-residual archive is listed, not copied.
"""
import hashlib
import os
import re
import shutil
import subprocess

os.environ["B_CPU_CAP"] = "none"

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "dist" / "B_release_v1.0.0"
PACK = ROOT / "manuscript" / "review_pack_srep"
FREEZE = "224779714bffd58b02059f51de48d92d1324a24f"
SAMPLE = re.compile(
    r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}|C3[NL]-[A-Z0-9]+|AP-[A-Z0-9]{4}\b"
)
SKIP_NAMES = {
    "cptac_paired_samples.tsv",
    "apollo_sample_map.tsv",
    "r1r2_model_compressed.npz",
    "r1r2_model.npz",
    "reliability_training_rows.tsv",
}
SKIP_PARTS = {"fake_null", "fake_planted", "_work", "__pycache__", ".venv", "data"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def blocked(path):
    if any(part in SKIP_PARTS for part in path.parts):
        return "excluded directory"
    if path.name in SKIP_NAMES:
        return "excluded name"
    if path.stat().st_size > 100 * 1024 * 1024:
        return "over 100 MB"
    if path.suffix == ".npz" and path.stat().st_size < 80 * 1024 * 1024:
        import numpy as np

        archive = np.load(path, allow_pickle=True)
        for key in archive.files:
            values = archive[key]
            if getattr(values, "dtype", None) is not None and values.dtype.kind in "USO":
                text = " ".join(map(str, values.ravel()[:5000]))
                if SAMPLE.search(text):
                    return "sample identifier"
    if path.suffix in {".tsv", ".csv", ".txt"} and path.stat().st_size < 80 * 1024 * 1024:
        sample = path.read_text(errors="replace")[:200000]
        if SAMPLE.search(sample):
            return "sample identifier"
    return ""


def copy_tree(rel, copied, skipped):
    source = ROOT / rel
    if source.is_file():
        files = [source]
    else:
        files = [path for path in source.rglob("*") if path.is_file()]
    for path in files:
        reason = blocked(path)
        if reason:
            skipped.append((str(path.relative_to(ROOT)), path.stat().st_size, reason))
            continue
        dest = DIST / path.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
        copied.append(dest)


def write_readme(dist):
    history = subprocess.check_output(
        ["git", "log", "--format=%H %ci %s", "--", "config/prereg_B4.md",
         "config/prereg_B6.md", "config/analysis_plan_B3.yaml",
         "config/analysis_plan_B4.yaml", "config/analysis_plan_B6.yaml",
         "config/analysis_plan_B10.yaml"],
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
        "# Project B release v1.0.0",
        "",
        "Local copy of the code, analysis plans, frozen outputs and Scientific Reports files for the mRNA–protein concordance forecast.",
        "",
        "Raw matrices are not in this folder. The manuscript Data availability section lists the public locations (LinkedOmicsKB, TCGA PanCanAtlas Xena, figshare file 34411172, GDC APOLLO, Gygi CCLE, DepMap 24Q4, TCGA PanCanAtlas RNA). This folder does not redistribute them.",
        "",
        "Also absent: sample identifiers, the planted-signal and null checks, and any file over 100 MB. The three large files are listed in `zenodo_only/files.tsv`.",
        "",
        "No license file is included. The author has not chosen one.",
        "Author names, the AI model name, funding, a GitHub URL and a Zenodo DOI are not filled in.",
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
        "## Preregistration commits",
        "",
        "```",
        history,
        "```",
        "",
        "Result tables keep their original column headers. `excluded.tsv` lists files that were not copied.",
        "",
    ]
    (dist / "README.md").write_text("\n".join(lines))


def review_pack():
    if PACK.exists():
        shutil.rmtree(PACK)
    PACK.mkdir(parents=True)
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep.md", PACK / "B_main_srep.md")
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep_claude.md", PACK / "B_main_srep_annotated.md")
    shutil.copy2(ROOT / "manuscript/srep/B_SI_srep.md", PACK / "B_SI_srep.md")
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
    import pandas as pd

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
    subprocess.check_call([
        pandoc, str(PACK / "figure_legends.md"), "-o", str(work / "legends.docx"),
    ])
    subprocess.check_call([
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b10pack",
        "--convert-to", "pdf", "--outdir", str(work), str(work / "legends.docx"),
    ])
    subprocess.check_call([
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b10pack",
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
    if DIST.exists():
        shutil.rmtree(DIST)
    DIST.mkdir(parents=True)
    copied = []
    skipped = []
    for rel in [
        "scripts",
        "config",
        "results/stage4/confirm",
        "results/stage5/posthoc",
        "results/stage5/atlas",
        "results/stage6/confirm",
        "results/stage6/posthoc",
        "results/stage6/prep",
        "results/stage6/frozen",
        "results/stage7",
        "results/stage8",
        "results/stage9",
        "results/stage10",
        "results/figures/source",
        "results/f1s_predictions.tsv",
        "results/reliability_predictions.tsv",
        "results/f1s_model_fit.tsv",
        "results/cptac_sample_counts.tsv",
        "results/cptac_set_sizes.tsv",
        "results/rppa_structure.tsv",
        "results/procan_sanger_by_lineage.tsv",
        "release/models/F1s_lightgbm.txt",
        "release/models/F1s_lightgbm_features.json",
        "release/model_cards.md",
        "manuscript/srep",
        "manuscript/checks",
        "submission_srep/figures",
        "submission_srep/B_manuscript_srep.docx",
        "submission_srep/B_SI_srep.pdf",
        "submission_srep/B_cover_letter_srep.docx",
        "submission_srep/checklist_srep.md",
    ]:
        copy_tree(rel, copied, skipped)
    freeze_dir = DIST / "scripts_freeze"
    freeze_dir.mkdir()
    text = subprocess.check_output(
        ["git", "show", f"{FREEZE}:scripts/08_stage3.py"], cwd=ROOT
    )
    (freeze_dir / "08_stage3.py").write_bytes(text)
    digest = hashlib.sha256(text).hexdigest()
    if digest != "b14c554c86f2243c52c9dba339ae0205bcf773f001a69e964c40bfbe90fea747":
        raise SystemExit(f"freeze hash {digest}")
    (DIST / "requirements.txt").write_text(
        "\n".join(
            [
                "numpy==1.26.4",
                "pandas==2.2.2",
                "scipy==1.13.1",
                "scikit-learn==1.4.2",
                "lightgbm==4.7.0",
                "matplotlib==3.8.4",
                "",
            ]
        )
    )
    large = [
        ("release/models/r1r2_model_compressed.npz", 121062358,
         "4c1b0f6f7a75cbb91b1043f01005d154e01ec2ba74d28da0402918beccf6d4cc"),
        ("results/stage6/frozen/r1r2_model.npz", 121062358,
         "4c1b0f6f7a75cbb91b1043f01005d154e01ec2ba74d28da0402918beccf6d4cc"),
        ("results/reliability_training_rows.tsv",
         (ROOT / "results/reliability_training_rows.tsv").stat().st_size,
         sha256(ROOT / "results/reliability_training_rows.tsv")),
    ]
    zdir = DIST / "zenodo_only"
    zdir.mkdir()
    lines = ["path\tbytes\tsha256\treason"]
    for path, size, digest in large:
        lines.append(f"{path}\t{size}\t{digest}\tover 100 MB; not in this folder")
    (zdir / "files.tsv").write_text("\n".join(lines) + "\n")
    write_readme(DIST)
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "figures", ROOT / "scripts" / "36_stage10_figures.py"
    )
    figures = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(figures)
    smoke = DIST / "_smoke"
    smoke.mkdir()
    figures.OUT = smoke
    figures.fig1()
    figures.fig2()
    figures.fig3()
    figures.fig4()
    figures.fig5()
    figures.fig_s1()
    figures.fig_s2()
    smoke_lines = ["figure\tbytes"]
    for path in sorted(smoke.glob("Fig*.pdf")):
        if path.stat().st_size < 1000:
            raise SystemExit(f"smoke figure too small: {path}")
        smoke_lines.append(f"{path.name}\t{path.stat().st_size}")
    (DIST / "smoke_redraw.tsv").write_text("\n".join(smoke_lines) + "\n")
    shutil.rmtree(smoke)
    skip_rows = ["path\tbytes\treason"]
    for path, size, reason in sorted(skipped):
        skip_rows.append(f"{path}\t{size}\t{reason}")
    (DIST / "excluded.tsv").write_text("\n".join(skip_rows) + "\n")
    rows = ["path\tbytes\tsha256"]
    for path in sorted(path for path in DIST.rglob("*") if path.is_file() and path.name != "MANIFEST.tsv"):
        rows.append(
            f"{path.relative_to(DIST).as_posix()}\t{path.stat().st_size}\t{sha256(path)}"
        )
    (DIST / "MANIFEST.tsv").write_text("\n".join(rows) + "\n")
    review_pack()
    pack_bytes = sum(path.stat().st_size for path in PACK.rglob("*") if path.is_file())
    if pack_bytes > 25 * 1024 * 1024:
        raise SystemExit(f"review pack {pack_bytes} bytes")
    (DIST / ".gitignore").write_text(
        "\n".join(
            [
                "*",
                "!*/",
                "!.gitignore",
                "!README.md",
                "!MANIFEST.tsv",
                "!requirements.txt",
                "!excluded.tsv",
                "!smoke_redraw.tsv",
                "!zenodo_only/**",
                "!scripts_freeze/**",
                "",
            ]
        )
    )
    print(f"copied {len(copied)} files")
    print(f"skipped {len(skipped)}")
    print(f"review pack {pack_bytes} bytes")
    over = [path for path in DIST.rglob("*") if path.is_file() and path.stat().st_size > 100 * 1024 * 1024]
    print(f"over 100 MB in dist: {len(over)}")


if __name__ == "__main__":
    main()
