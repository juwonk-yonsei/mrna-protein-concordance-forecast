"""Build review pack v4. No new release: the training-row comparison is OPTION X and no license or DOI was supplied.

B_CPU_CAP is set before any import.
"""

import os

os.environ["B_CPU_CAP"] = "none"

import shutil
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / "manuscript" / "review_pack_srep_v4"


def before_after():
    rows = [
        ["D1 note 2", "The LightGBM version used for the earlier fit was not recorded and we could not identify the cause.", "OPTION X: LightGBM 4.7.0 was recorded before and after the fits, and the cause was not identified."],
        ["D2 methods", "the two-fit paragraph from the previous revision", "kept"],
        ["D3 limitations", "covered only the 7,432 genes whose CPTAC correlations were stored at development", "covered only the 7,449 genes, the 7,432 of the P set and 17 further antibody-drug conjugate targets"],
    ]
    lines = ["| Block | Before | After |", "|---|---|---|"]
    for block, before, after in rows:
        lines.append(f"| {block} | {before} | {after} |")
    return "\n".join(lines) + "\n"


def main():
    if not (ROOT / "dist" / "B_release_v1.2.0").exists():
        raise SystemExit("v1.2.0 is missing")
    if (ROOT / "release/models/F1s_lightgbm_rppa_procan.txt").exists():
        raise SystemExit("OPTION X must not add a second model")
    if PACK.exists():
        shutil.rmtree(PACK)
    PACK.mkdir(parents=True)
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep.md", PACK / "B_main_srep.md")
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep_claude.md", PACK / "B_main_srep_annotated.md")
    shutil.copy2(ROOT / "manuscript/srep/B_SI_srep.md", PACK / "B_SI_srep.md")
    shutil.copy2(ROOT / "manuscript/claude/B_revision_B14_claude.md", PACK / "B_revision_B14_claude.md")
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
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b14pack",
        "--convert-to", "pdf", "--outdir", str(work), str(work / "legends.docx"),
    ])
    subprocess.check_call([
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b14pack",
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
    pack_bytes = sum(path.stat().st_size for path in PACK.rglob("*") if path.is_file())
    if pack_bytes > 25 * 1024 * 1024:
        raise SystemExit(f"review pack {pack_bytes} bytes")
    print(f"pack {pack_bytes} remaining {len(remaining)}", flush=True)


if __name__ == "__main__":
    main()
