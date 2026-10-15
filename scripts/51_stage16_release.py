"""Assemble release v1.3.0 and review pack v6. Earlier releases stay in place.

No second LightGBM file. Reloading the released model is compared with the frozen predictions and the measured gaps are written out.
"""

import os

os.environ["B_CPU_CAP"] = "none"

import hashlib
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    "stage10_release", Path(__file__).resolve().parents[1] / "scripts" / "39_stage10_release.py"
)
rel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rel)

ROOT = rel.ROOT
DIST = ROOT / "dist" / "B_release_v1.3.0"
PACK = ROOT / "manuscript" / "review_pack_srep_v6"
rel.DIST = DIST
rel.PACK = PACK

FEATURES = [
    "rho_mean", "rho_sd", "rho_min", "rho_max", "n_defined_cancers",
    "holdout_rna_mean_percentile", "holdout_rna_sd_percentile",
    "train_rna_mean_percentile_mean", "train_protein_rate_mean",
]


def write_readme(dist):
    history = subprocess.check_output(
        ["git", "log", "--format=%H %ci %s", "--",
         "config/prereg_B4.md", "config/prereg_B6.md",
         "config/analysis_plan_B3.yaml", "config/analysis_plan_B4.yaml",
         "config/analysis_plan_B6.yaml", "config/analysis_plan_B10.yaml",
         "config/analysis_plan_B11.yaml", "config/analysis_plan_B12.yaml",
         "config/analysis_plan_B13.yaml", "config/analysis_plan_B14.yaml",
         "config/analysis_plan_B15.yaml"],
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
        "# Project B release v1.3.0",
        "",
        "Local copy of the code, analysis plans, frozen outputs and Scientific Reports files for the mRNA–protein concordance forecast. v1.0.0, v1.1.0 and v1.2.0 are unchanged and are not in this folder.",
        "",
        "The code is under the MIT License. Derived result tables and model files are under Creative Commons Attribution 4.0 International. Raw data matrices are not redistributed.",
        "",
        "Raw matrices are not in this folder. The manuscript Data availability section lists the public locations. This folder does not redistribute them.",
        "",
        "Also absent: sample identifiers, the planted-signal and null checks, and any file over 100 MB. The three large files are listed in `zenodo_only/files.tsv`.",
        "",
        "Authors: Juwon Kang and Junjeong Choi. A GitHub URL and a Zenodo DOI are not filled in.",
        "",
        "## Forecast table v3",
        "",
        "`results/stage12/forecast_all_genes_v3.tsv.gz` adds `model_gene_features` to the v2 table. The other columns are the v2 columns.",
        "",
        "Of the 19,211 genes, 11,762 have model_gene_features 0: the gene-level features are the default values, so the forecast rests on the gene's RNA distribution in the query cohort.",
        "",
        "If model_gene_features is 0 and cptac_pooled_rho is available, the pooled CPTAC correlation was more accurate than the forecast in a post hoc comparison on APOLLO and should be preferred.",
        "",
        "## Forecast model",
        "",
        "One model file is released, `release/models/F1s_lightgbm.txt`. A post hoc grid over 1, 16 and 72 threads and over automatic, column-wise and row-wise histograms did not reproduce the frozen RPPA or ProCan predictions within 1e-12. The LightGBM version recorded before the RPPA and ProCan fit and after the APOLLO fit is 4.7.0, so no other version was installed. The cause was not identified, and no second model file was written.",
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


def section1_table():
    rows = [
        ["authors", "[AUTHOR TO FILL]", "Juwon Kang and Junjeong Choi, with ORCID, affiliation, corresponding author and e-mail on the main file; names and affiliation only on the SI"],
        ["LLM names", "[AUTHOR TO FILL: model and version]", "Claude Opus 5.5 (Anthropic); A coding agent in the Cursor editor running Grok 4.7 (xAI)"],
        ["funding", "[AUTHOR TO FILL: funding, with grant numbers.]", "This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors."],
        ["contributions", "[AUTHOR TO FILL]", "J.K. designed the study, wrote the analysis plans and the preregistrations, directed and reviewed the computational work, analysed the results and wrote the manuscript. J.C. supervised the study, reviewed the analysis plans and the manuscript, and provided critical revision. Both authors approved the submitted version."],
        ["competing interests", "The author(s) declare no competing interests. [AUTHOR TO CONFIRM]", "The authors declare no competing interests."],
        ["cover date", "a date placeholder", "[AUTHOR TO FILL: submission date]"],
        ["reviewers", "a suggested or excluded reviewer line", "deleted"],
        ["GitHub and Zenodo", "unfilled", "still unfilled"],
    ]
    lines = ["| Item | Before | After |", "|---|---|---|"]
    for item, before, after in rows:
        lines.append(f"| {item} | {before} | {after} |")
    return "\n".join(lines) + "\n"


def review_pack():
    if PACK.exists():
        shutil.rmtree(PACK)
    PACK.mkdir(parents=True)
    for name in ("B_main_srep.md", "B_SI_srep.md", "B_cover_letter_srep.md"):
        shutil.copy2(ROOT / "manuscript/srep" / name, PACK / name)
    shutil.copy2(ROOT / "manuscript/terms_of_use.md", PACK / "terms_of_use.md")
    for name in [
        "submission_srep/B_manuscript_srep.docx",
        "submission_srep/B_SI_srep.pdf",
        "submission_srep/B_cover_letter_srep.docx",
        "submission_srep/checklist_srep.md",
    ]:
        shutil.copy2(ROOT / name, PACK / Path(name).name)
    (PACK / "section1_before_after.md").write_text(section1_table())
    source = (ROOT / "scripts/08_stage3.py").read_text().splitlines()
    body = "\n".join(source[227:265])
    (PACK / "ridge_alpha_and_oof.md").write_text(
        "Verdict: ridge_alpha_and_oof does not modify its frame argument in place.\n\n"
        "```python\n" + body + "\n```\n"
    )
    ledger = pd.read_csv(ROOT / "manuscript/checks/ledger_srep.tsv", sep="\t")
    summary = ledger.groupby(["유형", "지위"]).size().reset_index(name="n")
    summary.to_csv(PACK / "ledger_summary.tsv", sep="\t", index=False)
    remaining = []
    for path in (ROOT / "manuscript/srep").glob("B_*srep.md"):
        if path.name.endswith("_claude.md"):
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if "[AUTHOR" in line or "[FILL" in line or "[CHECK" in line or "[OPTION" in line:
                remaining.append(f"{path.name}:{i}: {line.strip()[:240]}")
    (PACK / "remaining_placeholders.txt").write_text("\n".join(remaining) + "\n")
    fig_dir = PACK / "figures"
    fig_dir.mkdir()
    for path in (ROOT / "submission_srep" / "figures").glob("*.png"):
        shutil.copy2(path, fig_dir / path.name)


def predict_filled(booster, frame, med):
    filled = np.array(frame[FEATURES].to_numpy(dtype=np.float64), dtype=np.float64, copy=True)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(med, inds[1])
    return booster.predict(filled)


def reload_gap(booster):
    train = pd.read_csv(ROOT / "results/f1s_predictions.tsv", sep="\t")
    use = train[(train["in_P"] == 1) & np.isfinite(train["actual"].to_numpy(dtype=np.float64))]
    med = np.nanmedian(use[FEATURES].to_numpy(dtype=np.float64), axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    specs = [
        ("RPPA", ROOT / "results/stage9/cache_rppa_features.tsv", ROOT / "results/frozen_forecasts_rppa_F1s.tsv", ["xena_antibody", "cancer"]),
        ("ProCan", ROOT / "results/stage9/cache_procan_features.tsv", ROOT / "results/frozen_forecasts_procan_F1s.tsv", ["protein_id", "lineage"]),
    ]
    rows = []
    for dataset, features, frozen, keys in specs:
        frame = pd.read_csv(features, sep="\t")
        old = pd.read_csv(frozen, sep="\t")
        hat = frame[keys].copy()
        hat["f1s_lgbm"] = predict_filled(booster, frame, med)
        merged = hat.merge(old[keys + ["f1s_lgbm"]], on=keys, suffixes=("_new", "_old"))
        gap = float(np.max(np.abs(merged["f1s_lgbm_new"].to_numpy(dtype=np.float64) - merged["f1s_lgbm_old"].to_numpy(dtype=np.float64))))
        rows.append({"dataset": dataset, "n_aligned": int(len(merged)), "max_abs": gap})
        print("reload", dataset, gap, flush=True)
    sys.path.insert(0, str(ROOT))
    repro = importlib.import_module("scripts.43_stage12_repro")
    apollo = repro.apollo_features()
    old = pd.read_csv(ROOT / "results/stage6/frozen/apollo_gene_pred.tsv", sep="\t")
    hat = apollo[["gene_symbol"]].copy()
    hat["f1s_lgbm"] = predict_filled(booster, apollo, med)
    merged = hat.merge(old[["gene_symbol", "f1s_lgbm"]], on="gene_symbol", suffixes=("_new", "_old"))
    gap = float(np.max(np.abs(merged["f1s_lgbm_new"].to_numpy(dtype=np.float64) - merged["f1s_lgbm_old"].to_numpy(dtype=np.float64))))
    rows.append({"dataset": "APOLLO", "n_aligned": int(len(merged)), "max_abs": gap})
    print("reload", "APOLLO", gap, flush=True)
    out = ROOT / "results/stage16"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "reload_gap.tsv", sep="\t", index=False)
    return rows


def main():
    for name in ("B_release_v1.0.0", "B_release_v1.1.0", "B_release_v1.2.0"):
        if not (ROOT / "dist" / name).exists():
            raise SystemExit(f"{name} is missing")
    if DIST.exists():
        shutil.rmtree(DIST)
    DIST.mkdir(parents=True)
    copied, skipped = [], []
    for rel_path in [
        "scripts", "config",
        "results/stage4/confirm", "results/stage5/posthoc", "results/stage5/atlas",
        "results/stage6/confirm", "results/stage6/posthoc", "results/stage6/prep", "results/stage6/frozen",
        "results/stage7", "results/stage8", "results/stage9", "results/stage10",
        "results/stage11", "results/stage12", "results/stage15",
        "results/figures/source",
        "results/f1s_predictions.tsv", "results/reliability_predictions.tsv", "results/f1s_model_fit.tsv",
        "results/cptac_sample_counts.tsv", "results/cptac_set_sizes.tsv",
        "results/rppa_structure.tsv", "results/procan_sanger_by_lineage.tsv",
        "release/models/F1s_lightgbm.txt", "release/models/F1s_lightgbm_features.json",
        "release/model_cards.md", "release/NOTICE.md", "release/README.md",
        "manuscript/srep", "manuscript/claude", "manuscript/checks", "manuscript/terms_of_use.md",
        "LICENSE", "LICENSE-DATA",
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
    spec_fig = importlib.util.spec_from_file_location("figures", ROOT / "scripts" / "36_stage10_figures.py")
    figures = importlib.util.module_from_spec(spec_fig)
    spec_fig.loader.exec_module(figures)
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
    (DIST / ".gitignore").write_text("*\n!*/\n!.gitignore\n!README.md\n!MANIFEST.tsv\n!requirements.txt\n!excluded.tsv\n!smoke_redraw.tsv\n!zenodo_only/**\n!scripts_freeze/**\n")
    rows = ["path\tbytes\tsha256"]
    for path in sorted(p for p in DIST.rglob("*") if p.is_file() and p.name != "MANIFEST.tsv"):
        rows.append(f"{path.relative_to(DIST).as_posix()}\t{path.stat().st_size}\t{rel.sha256(path)}")
    manifest = DIST / "MANIFEST.tsv"
    manifest.write_text("\n".join(rows) + "\n")
    review_pack()
    pack_bytes = sum(path.stat().st_size for path in PACK.rglob("*") if path.is_file())
    if pack_bytes > 25 * 1024 * 1024:
        raise SystemExit(f"review pack {pack_bytes} bytes")
    readme = (DIST / "README.md").read_text()
    if "The code is under the MIT License" not in readme:
        raise SystemExit("README license paragraph missing")
    if (DIST / "release/models/F1s_lightgbm_rppa_procan.txt").exists():
        raise SystemExit("a second model must not be added")
    from lightgbm import Booster
    model_path = DIST / "release/models/F1s_lightgbm.txt"
    first = Booster(model_file=str(model_path))
    second = Booster(model_file=str(model_path))
    if first.model_to_string() != second.model_to_string():
        raise SystemExit("reloaded model text differs")
    if first.num_trees() != 300:
        raise SystemExit(f"tree count {first.num_trees()}")
    gaps = reload_gap(first)
    print(
        f"copied {len(copied)} skipped {len(skipped)} pack {pack_bytes} "
        f"reload_trees {first.num_trees()} manifest {rel.sha256(manifest)} gaps {gaps}",
        flush=True,
    )


if __name__ == "__main__":
    main()
