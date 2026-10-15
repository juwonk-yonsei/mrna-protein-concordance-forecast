"""Final figures 3, 5, and 6 at journal size.

Arial and Helvetica are not installed. Liberation Sans is the Arial-metric
substitute and is recorded as a judgment in the stage-7 report.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

from scripts.common import ROOT

font_manager.fontManager.addfont("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf")
font_manager.fontManager.addfont("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf")
plt.rcParams["font.family"] = "Liberation Sans"
plt.rcParams["font.size"] = 8
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["axes.titlesize"] = 8
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = False

OUT = ROOT / "results" / "figures"
SRC = OUT / "source"
MM = 1 / 25.4
DOUBLE = 175 * MM


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)


def write(name, frame):
    SRC.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SRC / name, sep="\t", index=False)


def fig3():
    diff = pd.read_csv(ROOT / "results" / "stage7" / "f0_diff.tsv", sep="\t")
    b4 = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "hypotheses.tsv", sep="\t")
    b6 = pd.read_csv(ROOT / "results" / "stage6" / "confirm" / "hypotheses.tsv", sep="\t")
    raw = pd.read_csv(ROOT / "results" / "stage6" / "posthoc" / "posthoc_05_raw.tsv", sep="\t")
    raw_map = dict(zip(raw["item"], raw["value"]))
    rows = []
    for name, cohort in (("H3", "RPPA"), ("H4", "ProCan")):
        rec = b4[b4["hypothesis"] == name].iloc[0]
        rows.append({"status": "preregistered", "cohort": cohort, "which": "F1s-LightGBM", "value": float(rec["S_F"]), "low": float(rec["S_F_ci_low"]), "high": float(rec["S_F_ci_high"]), "n": int(rec["n_pairs"])})
        rows.append({"status": "preregistered", "cohort": cohort, "which": "B_sd", "value": float(rec["S_B"]), "low": float(rec["S_B_ci_low"]), "high": float(rec["S_B_ci_high"]), "n": int(rec["n_pairs"])})
    rec = b6[b6["hypothesis"] == "P2"].iloc[0]
    rows.append({"status": "preregistered", "cohort": "APOLLO", "which": "F1s-LightGBM", "value": float(rec["S_F"]), "low": float(rec["S_F_ci_low"]), "high": float(rec["S_F_ci_high"]), "n": int(rec["n"])})
    rows.append({"status": "preregistered", "cohort": "APOLLO", "which": "B_sd", "value": float(rec["S_B"]), "low": float(rec["S_B_ci_low"]), "high": float(rec["S_B_ci_high"]), "n": int(rec["n"])})
    for cohort in ("RPPA", "ProCan", "APOLLO"):
        base = diff[(diff["dataset"] == cohort) & (diff["predictor"] == "f1s_lgbm") & (diff["stratum"] == "all")].iloc[0]
        rows.append({"status": "post hoc", "cohort": cohort, "which": "F0", "value": float(base["s_base"]), "low": np.nan, "high": np.nan, "n": int(base["n_pairs"])})
    rows.append({"status": "post hoc", "cohort": "APOLLO non-imputed", "which": "F1s-LightGBM", "value": float(raw_map["P2_S_F"]), "low": np.nan, "high": np.nan, "n": int(float(raw_map["P2_n"]))})
    rows.append({"status": "post hoc", "cohort": "APOLLO non-imputed", "which": "B_sd", "value": float(raw_map["P2_S_B"]), "low": np.nan, "high": np.nan, "n": int(float(raw_map["P2_n"]))})
    points = pd.DataFrame(rows)
    write("fig3_final_points.tsv", points)
    write("fig3_final_diff.tsv", diff)
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 115 * MM), gridspec_kw={"width_ratios": [1, 1.15]})
    cohorts = ["RPPA", "ProCan", "APOLLO", "APOLLO non-imputed"]
    colors = {"F1s-LightGBM": "#4C78A8", "B_sd": "#F58518", "F0": "#54A24B"}
    shifts = {"F1s-LightGBM": 0.18, "B_sd": 0.0, "F0": -0.18}
    ax = axes[0]
    for which, color in colors.items():
        sub = points[points["which"] == which]
        first = True
        for i, cohort in enumerate(cohorts):
            hit = sub[sub["cohort"] == cohort]
            if hit.empty:
                continue
            rec = hit.iloc[0]
            yy = len(cohorts) - 1 - i + shifts[which]
            marker = "o" if rec["status"] == "preregistered" else "s"
            ax.plot(rec["value"], yy, marker, color=color, ms=4, label=which if first else None)
            first = False
            if np.isfinite(rec["low"]):
                ax.plot([rec["low"], rec["high"]], [yy, yy], color=color, lw=0.8)
    ax.set_yticks(range(len(cohorts))[::-1], cohorts)
    ax.set_xlabel("Spearman with observed rho")
    ax.set_xlim(0, 0.9)
    ax = axes[1]
    order = []
    for dataset in ("RPPA", "ProCan", "APOLLO"):
        for stratum in ("all", "in_P", "outside_P"):
            order.append((dataset, stratum))
    ylabels = []
    for i, (dataset, stratum) in enumerate(order):
        yy = len(order) - 1 - i
        ylabels.append(f"{dataset} {stratum}")
        for pred, color, shift in (("f1s_lgbm", "#4C78A8", 0.12), ("f1s_ridge", "#B279A2", -0.12)):
            rec = diff[(diff["dataset"] == dataset) & (diff["predictor"] == pred) & (diff["stratum"] == stratum)].iloc[0]
            ax.plot(rec["diff"], yy + shift, "o", color=color, ms=3.5, label=pred if i == 0 else None)
            ax.plot([rec["ci_low"], rec["ci_high"]], [yy + shift, yy + shift], color=color, lw=0.8)
    ax.axvline(0, color="gray", lw=0.6)
    ax.set_yticks(range(len(order))[::-1], ylabels)
    ax.set_xlabel("Spearman(F1s) - Spearman(F0)")
    ax.set_title("Post hoc, unit bootstrap 95% interval")
    ax.legend(frameon=False, fontsize=7)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.28, 0.98))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "fig3")


def fig5():
    table = pd.read_csv(ROOT / "results" / "stage7" / "calibration_table.tsv", sep="\t")
    table = table[table["predictor"] == "f1s_lgbm"].copy()
    write("fig5_final_calibration.tsv", table)
    methods = [("none", "#4C78A8", "o"), ("delta", "#F58518", "s"), ("delta_quantile", "#54A24B", "D")]
    status = {"rppa": "post hoc", "procan": "post hoc", "APOLLO": "pre-specified descriptive"}
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 78 * MM), sharey=True)
    for ax, platform in zip(axes, ["rppa", "procan", "APOLLO"]):
        sub = table[table["platform"] == platform]
        ax.axhspan(0.75, 0.85, color="#dddddd", zorder=0)
        ax.axhline(0.80, color="gray", lw=0.5, ls="--")
        for method, color, marker in methods:
            hit = sub[sub["method"] == method].sort_values("K")
            ax.errorbar(
                hit["K"],
                hit["cov80_mean"],
                yerr=[hit["cov80_mean"] - hit["cov80_p2_5"], hit["cov80_p97_5"] - hit["cov80_mean"]],
                fmt=marker + "-",
                color=color,
                ms=3.5,
                lw=0.8,
                capsize=2,
                label=method,
            )
        ax.set_xscale("log")
        ax.set_xticks([10, 20, 50, 100], ["10", "20", "50", "100"])
        ax.set_title(f"{platform}\n{status[platform]}")
        ax.set_xlabel("K")
        ax.set_ylim(0.45, 1.0)
    axes[0].set_ylabel("80% coverage")
    axes[2].legend(frameon=False, fontsize=6, title="F1s-LightGBM")
    fig.tight_layout()
    save(fig, "fig5")


def fig6():
    atlas = pd.read_csv(ROOT / "results" / "stage5" / "atlas" / "adc_atlas_v2.tsv", sep="\t")
    conc_tab = pd.read_csv(ROOT / "results" / "stage7" / "adc_concordance.tsv", sep="\t")
    pairs = pd.read_csv(ROOT / "results" / "stage5" / "atlas" / "concordance_pairs.tsv", sep="\t")
    adc = pd.read_csv(ROOT / "results" / "stage6" / "confirm" / "descriptive_adc.tsv", sep="\t")
    extra = adc[np.isfinite(adc["actual"]) & np.isfinite(adc["rho_hat"])][["gene_symbol", "actual", "rho_hat"]].copy()
    extra.insert(0, "context", "LUAD")
    extra.insert(0, "source", "APOLLO")
    conc = pd.concat([pairs[["source", "gene_symbol", "context", "rho_hat", "actual"]], extra], ignore_index=True)
    conc.insert(0, "status", "post hoc")
    write("fig6_final_concordance.tsv", conc)
    mat = atlas.pivot(index="gene_symbol", columns="disease", values="rho_hat")
    genes = list(mat.index)
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 150 * MM), gridspec_kw={"width_ratios": [1.35, 1]})
    im = axes[0].imshow(mat.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    axes[0].set_yticks(range(len(genes)))
    labels = []
    for gene in genes:
        labels.append(gene + " *" if gene == "CLDN18" else gene)
    axes[0].set_yticklabels(labels, fontsize=6)
    for tick, gene in zip(axes[0].get_yticklabels(), genes):
        if gene == "CLDN18":
            tick.set_color("#E45756")
    axes[0].set_xticks(range(len(mat.columns)))
    axes[0].set_xticklabels(mat.columns, rotation=90, fontsize=5)
    axes[0].set_title("Post hoc rho_hat")
    fig.colorbar(im, ax=axes[0], fraction=0.03)
    colors = {"CPTAC": "#4C78A8", "RPPA": "#F58518", "ProCan": "#54A24B", "ProCan_frozen_forecast": "#72B7B2", "APOLLO": "#E45756"}
    ax = axes[1]
    for source, sub in conc.groupby("source"):
        rec = conc_tab[conc_tab["source"] == source]
        rho = float(rec["spearman_rho_hat_actual"].iloc[0]) if len(rec) else np.nan
        n = int(rec["n_pairs"].iloc[0]) if len(rec) else len(sub)
        ax.scatter(sub["rho_hat"], sub["actual"], s=12, c=colors.get(source, "gray"), label=f"{source} n={n}, rho={rho:.3f}", alpha=0.8)
    cldn = extra[extra["gene_symbol"] == "CLDN18"]
    if len(cldn):
        ax.scatter(cldn["rho_hat"], cldn["actual"], s=28, facecolors="none", edgecolors="#E45756", linewidths=0.8, label="CLDN18 APOLLO")
    ax.plot([-0.2, 1.2], [-0.2, 1.2], color="gray", lw=0.5)
    ax.set_xlabel("rho_hat")
    ax.set_ylabel("observed rho")
    ax.set_title("Post hoc concordance")
    ax.legend(frameon=False, fontsize=5.5, loc="upper left")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.text(
        0.01,
        0.006,
        "* CLDN18: gene-level assays cannot separate isoforms; the ADC target is CLDN18.2.",
        fontsize=6,
        color="#E45756",
    )
    save(fig, "fig6")


if __name__ == "__main__":
    fig3()
    fig5()
    fig6()
    print("figures 3 5 6 written", flush=True)
