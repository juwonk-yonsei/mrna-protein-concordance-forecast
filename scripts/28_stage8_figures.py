"""Redraw figures 1, 2, and 4 at journal size, and add F1s-L* to figure 3.

Figure content for 1, 2, and 4 is read from the existing source tables.
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
SINGLE = 85 * MM


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)


def fig1():
    delta = pd.read_csv(SRC / "fig1_loco_delta.tsv", sep="\t")
    rel = pd.read_csv(SRC / "fig1_loco_reliability.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 80 * MM))
    axes[0].bar(delta["cancer"], delta["mean_delta"], color="#4C78A8")
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_ylabel("Mean delta (rho_AR - rho_S0)")
    axes[0].set_title("Development")
    axes[0].tick_params(axis="x", rotation=60, labelsize=6)
    wide = rel.pivot(index="cancer", columns="predictor", values="spearman").reset_index()
    x = np.arange(len(wide))
    axes[1].bar(x - 0.18, wide["F1s-LightGBM"], width=0.36, label="F1s-LightGBM", color="#4C78A8")
    axes[1].bar(x + 0.18, wide["B_sd"], width=0.36, label="B_sd", color="#F58518")
    axes[1].set_xticks(x, wide["cancer"], rotation=60, fontsize=6)
    axes[1].set_ylabel("Spearman with observed rho")
    axes[1].set_title("Development")
    axes[1].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    save(fig, "fig1")


def fig2():
    frame = pd.read_csv(SRC / "fig2_effects.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(DOUBLE, 110 * MM))
    y = np.arange(len(frame))[::-1]
    for i, rec in enumerate(frame.itertuples(index=False)):
        color = "#4C78A8" if rec.status == "사전등록" else "#B8B8B8"
        marker = "o" if rec.status == "사전등록" else "s"
        ax.plot(rec.mean, y[i], marker, color=color, ms=4)
        if np.isfinite(rec.low):
            ax.plot([rec.low, rec.high], [y[i], y[i]], color=color, lw=0.8)
    ax.axvline(0, color="black", lw=0.6)
    labels = []
    for rec in frame.itertuples(index=False):
        status = "preregistered" if rec.status == "사전등록" else "post hoc"
        labels.append(f"{rec.item} {rec.contrast} ({status})")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Mean unit difference")
    fig.tight_layout()
    save(fig, "fig2")


def fig3():
    points = pd.read_csv(SRC / "fig3_final_points.tsv", sep="\t")
    diff = pd.read_csv(SRC / "fig3_final_diff.tsv", sep="\t")
    look = pd.read_csv(ROOT / "results" / "stage8" / "spearman_diff.tsv", sep="\t")
    look = look[(look["left"] == "f1s_lgbm") & (look["right"] == "Lstar") & (look["stratum"].isin(["all", "L1"]))]
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 125 * MM), gridspec_kw={"width_ratios": [1, 1.25]})
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
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.28, 0.98))
    ax = axes[1]
    order = []
    for dataset in ("RPPA", "ProCan", "APOLLO"):
        for stratum in ("all", "in_P", "outside_P"):
            order.append(("f0", dataset, stratum))
        for stratum in ("all", "L1"):
            order.append(("Lstar", dataset, stratum))
    ylabels = []
    for i, (kind, dataset, stratum) in enumerate(order):
        yy = len(order) - 1 - i
        ylabels.append(f"{dataset} {stratum} vs {kind}")
        if kind == "f0":
            for pred, color, shift in (("f1s_lgbm", "#4C78A8", 0.12), ("f1s_ridge", "#B279A2", -0.12)):
                rec = diff[(diff["dataset"] == dataset) & (diff["predictor"] == pred) & (diff["stratum"] == stratum)].iloc[0]
                ax.plot(rec["diff"], yy + shift, "o", color=color, ms=3, label=pred if i == 0 else None)
                ax.plot([rec["ci_low"], rec["ci_high"]], [yy + shift, yy + shift], color=color, lw=0.7)
        else:
            rec = look[(look["dataset"] == dataset) & (look["stratum"] == stratum)].iloc[0]
            ax.plot(rec["diff"], yy, "D", color="#E45756", ms=3.5, label="f1s_lgbm - L*" if i == 3 else None)
            if np.isfinite(rec["ci_low"]):
                ax.plot([rec["ci_low"], rec["ci_high"]], [yy, yy], color="#E45756", lw=0.7)
    ax.axvline(0, color="gray", lw=0.6)
    ax.set_yticks(range(len(order))[::-1], ylabels, fontsize=6)
    ax.set_xlabel("Spearman difference")
    ax.set_title("Post hoc")
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "fig3")


def fig4():
    frame = pd.read_csv(SRC / "fig4_scope.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(SINGLE, 80 * MM))
    post = frame[frame["status"] == "사후"]
    pre = frame[frame["status"] == "사전등록"]
    ax.scatter(post["s"], post["median_delta"], c="#4C78A8", s=16, label=f"post hoc n={len(post)}")
    ax.scatter(pre["s"], pre["median_delta"], c="#E45756", marker="D", s=28, label="preregistered APOLLO")
    ax.axvline(0.8830509478307224, color="gray", lw=0.6, ls="--")
    ax.axvline(0.9788744329350166, color="gray", lw=0.6, ls="--")
    ax.set_xlabel("s")
    ax.set_ylabel("Median delta")
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    save(fig, "fig4")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    print("figures 1-4 restyled", flush=True)
