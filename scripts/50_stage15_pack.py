"""Build review pack v5. No new release: the CSV round trip is OPTION Z and no license or DOI was supplied.

B_CPU_CAP is set before any import.
"""

import os

os.environ["B_CPU_CAP"] = "none"

import shutil
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / "manuscript" / "review_pack_srep_v5"


def before_after():
    rows = [
        ["H1 note 2", "We could not identify the cause.", "kept; OPTION Z"],
        ["H2 methods", "the cause was not identified (Supplementary Note 2).", "kept; OPTION Z"],
        ["H3 note 2", "All results in the paper use the frozen predictions, which are released.", "same sentence, then: Code that fits a model on a table in memory and later reads the same table from disk should treat the two as different inputs."],
    ]
    lines = ["| Block | Before | After |", "|---|---|---|"]
    for block, before, after in rows:
        lines.append(f"| {block} | {before} | {after} |")
    return "\n".join(lines) + "\n"


def main():
    for name in ("B_release_v1.0.0", "B_release_v1.1.0", "B_release_v1.2.0"):
        if not (ROOT / "dist" / name).exists():
            raise SystemExit(f"{name} is missing")
    if (ROOT / "dist" / "B_release_v1.3.0").exists():
        raise SystemExit("v1.3.0 must not be created for OPTION Z")
    if (ROOT / "release/models/F1s_lightgbm_rppa_procan.txt").exists():
        raise SystemExit("OPTION Z must not add a second model")
    if PACK.exists():
        shutil.rmtree(PACK)
    PACK.mkdir(parents=True)
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep.md", PACK / "B_main_srep.md")
    shutil.copy2(ROOT / "manuscript/srep/B_main_srep_claude.md", PACK / "B_main_srep_annotated.md")
    shutil.copy2(ROOT / "manuscript/srep/B_SI_srep.md", PACK / "B_SI_srep.md")
    shutil.copy2(ROOT / "manuscript/claude/B_revision_B15_claude.md", PACK / "B_revision_B15_claude.md")
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
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b15pack",
        "--convert-to", "pdf", "--outdir", str(work), str(work / "legends.docx"),
    ])
    subprocess.check_call([
        soffice, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b15pack",
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
