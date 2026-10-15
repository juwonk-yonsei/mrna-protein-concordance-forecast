"""Refresh v1.3.0 figures and the manifest, and build review pack v7.

v1.0.0, v1.1.0 and v1.2.0 stay in place. No second model file is written.
"""

import os

os.environ["B_CPU_CAP"] = "none"

import importlib.util
import shutil
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "dist" / "B_release_v1.3.0"
PACK = ROOT / "manuscript" / "review_pack_srep_v7"
spec = importlib.util.spec_from_file_location("stage10_release", ROOT / "scripts" / "39_stage10_release.py")
rel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rel)


def copy_over(rel_path):
    source = ROOT / rel_path
    dest = DIST / rel_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(source, dest)
    else:
        shutil.copy2(source, dest)


def side_by_side():
    before = ROOT / "results" / "figures_pre_c17" / "submission_srep_figures"
    after = ROOT / "submission_srep" / "figures"
    work = PACK / "_work"
    work.mkdir(parents=True)
    pages = []
    for name in ["Fig1", "Fig2", "Fig3", "Fig4", "Fig5", "FigS1", "FigS2"]:
        text = work / f"{name}.md"
        text.write_text(
            f"# {name}\n\n"
            f"Before\n\n![]({before / (name + '.png')})\n\n"
            f"After\n\n![]({after / (name + '.png')})\n"
        )
        docx = work / f"{name}.docx"
        subprocess.check_call(["/tmp/pandoc-3.6.4/bin/pandoc", str(text), "-o", str(docx)])
        subprocess.check_call([
            "/tmp/lo/root/opt/libreoffice26.8/program/soffice",
            "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b17pack",
            "--convert-to", "pdf", "--outdir", str(work), str(docx),
        ])
        pages.append(work / f"{name}.pdf")
    subprocess.check_call(["pdfunite", *[str(path) for path in pages], str(PACK / "figures_before_after.pdf")])


def main():
    if not DIST.exists():
        raise SystemExit("v1.3.0 is missing")
    for rel_path in [
        "scripts/36_stage10_figures.py",
        "scripts/37_stage10_ledger.py",
        "scripts/38_stage10_docx.py",
        "scripts/51_stage16_release.py",
        "scripts/52_stage17_refresh.py",
        "manuscript/srep/B_main_srep.md",
        "manuscript/srep/B_SI_srep.md",
        "manuscript/srep/B_cover_letter_srep.md",
        "manuscript/checks/ledger_srep.tsv",
        "manuscript/terms_of_use.md",
        "submission_srep/figures",
        "submission_srep/B_manuscript_srep.docx",
        "submission_srep/B_SI_srep.pdf",
        "submission_srep/B_cover_letter_srep.docx",
        "submission_srep/checklist_srep.md",
        "results/stage16/reload_gap.tsv",
        "results/stage17/figure_format.tsv",
    ]:
        copy_over(rel_path)
    spec_fig = importlib.util.spec_from_file_location("figures", ROOT / "scripts" / "36_stage10_figures.py")
    figures = importlib.util.module_from_spec(spec_fig)
    spec_fig.loader.exec_module(figures)
    smoke = DIST / "_smoke"
    smoke.mkdir()
    figures.OUT = smoke
    figures.FORMAT = []
    for name in ("fig1", "fig2", "fig3", "fig4", "fig5", "fig_s1", "fig_s2"):
        getattr(figures, name)()
    smoke_lines = ["figure\tbytes"]
    for path in sorted(smoke.glob("Fig*.pdf")):
        if path.stat().st_size < 1000:
            raise SystemExit(f"smoke figure too small: {path}")
        smoke_lines.append(f"{path.name}\t{path.stat().st_size}")
    (DIST / "smoke_redraw.tsv").write_text("\n".join(smoke_lines) + "\n")
    shutil.rmtree(smoke)
    if (DIST / "MANIFEST.tsv").exists():
        (DIST / "MANIFEST.tsv").unlink()
    rows = ["path\tbytes\tsha256"]
    for path in sorted(p for p in DIST.rglob("*") if p.is_file()):
        rows.append(f"{path.relative_to(DIST).as_posix()}\t{path.stat().st_size}\t{rel.sha256(path)}")
    manifest = DIST / "MANIFEST.tsv"
    manifest.write_text("\n".join(rows) + "\n")
    if PACK.exists():
        shutil.rmtree(PACK)
    PACK.mkdir(parents=True)
    for name in ("B_main_srep.md", "B_SI_srep.md", "B_cover_letter_srep.md"):
        shutil.copy2(ROOT / "manuscript/srep" / name, PACK / name)
    for name in [
        "submission_srep/B_manuscript_srep.docx",
        "submission_srep/B_SI_srep.pdf",
        "submission_srep/B_cover_letter_srep.docx",
        "submission_srep/checklist_srep.md",
    ]:
        shutil.copy2(ROOT / name, PACK / Path(name).name)
    main_text = (ROOT / "manuscript/srep/B_main_srep.md").read_text()
    legends = main_text.split("## Figure legends", 1)[1].split("## Tables", 1)[0]
    (PACK / "figure_legends.md").write_text("## Figure legends\n" + legends)
    si = (ROOT / "manuscript/srep/B_SI_srep.md").read_text().splitlines()
    si_legends = "\n".join(line for line in si if line.startswith("**Supplementary Figure"))
    (PACK / "supplementary_figure_legends.md").write_text(si_legends + "\n")
    hashes = ["file\tsha256_before\tsha256_after"]
    before = (ROOT / "results/stage17/source_sha256_before.tsv").read_text().splitlines()[1:]
    after = {line.split("\t", 1)[0]: line.split("\t", 1)[1] for line in (ROOT / "results/stage17/source_sha256_after.tsv").read_text().splitlines()[1:]}
    for line in before:
        path, digest = line.split("\t", 1)
        hashes.append(f"{path}\t{digest}\t{after[path]}")
    (PACK / "source_sha256.tsv").write_text("\n".join(hashes) + "\n")
    shutil.copy2(ROOT / "submission_srep/checklist_srep.md", PACK / "checklist_srep.md")
    side_by_side()
    shutil.rmtree(PACK / "_work")
    fig_dir = PACK / "figures"
    fig_dir.mkdir()
    for path in (ROOT / "submission_srep" / "figures").glob("*.png"):
        shutil.copy2(path, fig_dir / path.name)
    pack_bytes = sum(path.stat().st_size for path in PACK.rglob("*") if path.is_file())
    if pack_bytes > 25 * 1024 * 1024:
        raise SystemExit(f"review pack {pack_bytes} bytes")
    if (DIST / "release/models/F1s_lightgbm_rppa_procan.txt").exists():
        raise SystemExit("a second model must not be added")
    print(f"pack {pack_bytes} manifest {rel.sha256(manifest)}", flush=True)


if __name__ == "__main__":
    main()
