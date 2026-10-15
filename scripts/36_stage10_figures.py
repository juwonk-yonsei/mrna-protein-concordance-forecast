"""Draw the Scientific Reports figures from stored result tables."""
import os

os.environ["B_CPU_CAP"] = "none"

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "submission_srep" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
SRC = ROOT / "results" / "figures" / "source"
FORMAT = []

font_manager.fontManager.addfont("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf")
font_manager.fontManager.addfont("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf")
plt.rcParams["font.family"] = "Liberation Sans"
plt.rcParams["font.size"] = 8
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = False

MM = 1 / 25.4
DOUBLE = 175 * MM


def letter(ax, name):
    ax.set_title(name, loc="left", fontweight="bold", fontsize=9, pad=2)


def save(fig, name):
    fig.canvas.draw()
    sizes = []
    for ax in fig.axes:
        sizes.append(float(ax.xaxis.label.get_size()))
        sizes.append(float(ax.yaxis.label.get_size()))
        for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            sizes.append(float(lab.get_size()))
    w_in, h_in = fig.get_size_inches()
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.tif", dpi=300)
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)
    FORMAT.append({
        "figure": name,
        "width_mm": w_in * 25.4,
        "height_mm": h_in * 25.4,
        "tiff_width_px": int(round(w_in * 300)),
        "tiff_height_px": int(round(h_in * 300)),
        "png_width_px": int(round(w_in * 200)),
        "png_height_px": int(round(h_in * 200)),
        "tiff_dpi": 300,
        "png_dpi": 200,
        "min_axis_font_pt": min(sizes) if sizes else "",
    })


def fig1():
    delta = pd.read_csv(SRC / "fig1_loco_delta.tsv", sep="\t")
    rel = pd.read_csv(SRC / "fig1_loco_reliability.tsv", sep="\t")
    wide = rel.pivot(index="cancer", columns="predictor", values="spearman").reset_index()
    order = wide.sort_values("F1s-LightGBM", ascending=False)["cancer"].tolist()
    wide = wide.set_index("cancer").loc[order].reset_index()
    delta = delta.set_index("cancer").loc[order].reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 88 * MM))
    x = np.arange(len(wide))
    axes[0].bar(x - 0.18, wide["F1s-LightGBM"], width=0.36, label="Forecast", color="#4C78A8")
    axes[0].bar(x + 0.18, wide["B_sd"], width=0.36, label="RNA variance", color="#F58518")
    axes[0].set_xticks(x, wide["cancer"], rotation=45, ha="right")
    axes[0].set_ylabel("Spearman correlation with observed ρ")
    axes[0].legend(frameon=False, fontsize=7, loc="upper right")
    axes[1].bar(delta["cancer"], delta["mean_delta"], color="#B8B8B8")
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].set_ylabel("Mean ρ_AR − ρ_S0")
    axes[1].tick_params(axis="x", rotation=45)
    for lab in axes[1].get_xticklabels():
        lab.set_ha("right")
    letter(axes[0], "a")
    letter(axes[1], "b")
    fig.tight_layout()
    save(fig, "Fig1")


def span(rows):
    vals = []
    for value, low, high in rows:
        vals.append(float(value))
        if np.isfinite(low):
            vals.append(float(low))
        if np.isfinite(high):
            vals.append(float(high))
    lo, hi = min(vals), max(vals)
    pad = 0.08 * (hi - lo if hi > lo else 1.0)
    return lo - pad, hi + pad


def draw_rows(ax, rows):
    y = np.arange(len(rows))[::-1]
    for i, (label, value, low, high, marker) in enumerate(rows):
        color = "#4C78A8" if marker == "o" else "#E45756"
        ax.plot(value, y[i], marker, color=color, ms=4)
        if np.isfinite(low) and np.isfinite(high):
            ax.plot([low, high], [y[i], y[i]], color=color, lw=0.8)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y, [row[0] for row in rows], fontsize=7)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(*span([(row[1], row[2], row[3]) for row in rows]))
    ax.set_xlabel("Difference in Spearman correlation")


def fig2():
    points = pd.read_csv(SRC / "fig3_final_points.tsv", sep="\t")
    points = points[points["cohort"] != "APOLLO non-imputed"]
    spear = pd.read_csv(ROOT / "results" / "stage8" / "spearman.tsv", sep="\t")
    diff = pd.read_csv(ROOT / "results" / "stage8" / "spearman_diff.tsv", sep="\t")
    hyp4 = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "hypotheses.tsv", sep="\t")
    hyp6 = pd.read_csv(ROOT / "results" / "stage6" / "confirm" / "hypotheses.tsv", sep="\t")
    fig = plt.figure(figsize=(DOUBLE, 132 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.6], width_ratios=[1, 1.25], hspace=0.55, wspace=0.85)
    ax = fig.add_subplot(gs[:, 0])
    cohorts = ["RPPA", "ProCan", "APOLLO"]
    series = [
        ("F1s-LightGBM", "o", "#4C78A8", 0.18, True, "Forecast (preregistered)"),
        ("B_sd", "o", "#F58518", 0.06, True, "RNA variance (preregistered)"),
        ("F0", "s", "#54A24B", -0.06, False, "CPTAC mean lookup (post hoc)"),
    ]
    for which, marker, color, shift, interval, label in series:
        for i, cohort in enumerate(cohorts):
            rec = points[(points["which"] == which) & (points["cohort"] == cohort)].iloc[0]
            yy = len(cohorts) - 1 - i + shift
            ax.plot(rec["value"], yy, marker, color=color, ms=4)
            if interval and np.isfinite(rec["low"]):
                ax.plot([rec["low"], rec["high"]], [yy, yy], color=color, lw=0.8)
    for i, cohort in enumerate(cohorts):
        rec = spear[(spear.dataset == cohort) & (spear.stratum == "all") & (spear.score == "Lstar")].iloc[0]
        ax.plot(rec["spearman"], len(cohorts) - 1 - i - 0.18, "s", color="#B279A2", ms=4)
    ax.plot([], [], "o", color="#4C78A8", label="Forecast (preregistered)")
    ax.plot([], [], "o", color="#F58518", label="RNA variance (preregistered)")
    ax.plot([], [], "s", color="#54A24B", label="CPTAC mean lookup (post hoc)")
    ax.plot([], [], "s", color="#B279A2", label="Cascade lookup (post hoc)")
    ax.set_yticks(range(len(cohorts))[::-1], cohorts)
    ax.set_xlabel("Spearman correlation with observed ρ")
    ax.set_xlim(0, 0.9)
    handles, labels = ax.get_legend_handles_labels()
    letter(ax, "a")
    pre_rows = []
    for cohort, key in (("RPPA", "H3"), ("ProCan", "H4"), ("APOLLO", "P2")):
        src = hyp4 if key.startswith("H") else hyp6
        rec = src[src.hypothesis == key].iloc[0]
        pre_rows.append((f"{cohort} forecast − RNA variance", rec["diff"], rec["diff_ci_low"], rec["diff_ci_high"], "o"))
    post_rows = []
    for stratum, label in (("all", "forecast − cascade, all"), ("L1", "forecast − cascade, matched")):
        for cohort in cohorts:
            rec = diff[(diff.dataset == cohort) & (diff.stratum == stratum) & (diff.left == "f1s_lgbm") & (diff.right == "Lstar")].iloc[0]
            post_rows.append((f"{cohort} {label}", rec["diff"], rec["ci_low"], rec["ci_high"], "D"))
    for cohort in cohorts:
        rec = diff[(diff.dataset == cohort) & (diff.stratum == "all") & (diff.left == "HYB") & (diff.right == "f1s_lgbm")].iloc[0]
        post_rows.append((f"{cohort} hybrid − forecast", rec["diff"], rec["ci_low"], rec["ci_high"], "D"))
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    draw_rows(ax_b, pre_rows)
    draw_rows(ax_c, post_rows)
    ax_b.set_title("Preregistered", loc="center", fontsize=8)
    ax_c.set_title("Post hoc", loc="center", fontsize=8)
    ax_b.yaxis.tick_right()
    ax_c.yaxis.tick_right()
    letter(ax_b, "b")
    letter(ax_c, "c")
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.06, 0.01), frameon=False, fontsize=7, ncol=2)
    fig.subplots_adjust(left=0.10, right=0.68, bottom=0.18, top=0.94, hspace=0.55, wspace=0.12)
    save(fig, "Fig2")


def fig3():
    cal = pd.read_csv(ROOT / "results" / "stage7" / "calibration_table.tsv", sep="\t")
    cal = cal[cal.predictor == "f1s_lgbm"]
    methods = [("none", "No adjustment", "#4C78A8"), ("delta", "Shift (δ)", "#F58518"), ("delta_quantile", "Shift + quantile", "#54A24B")]
    offsets = {"none": -1.6, "delta": 0.0, "delta_quantile": 1.6}
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 78 * MM), sharey=True)
    lows, highs = [], []
    for ax, platform, title, name in zip(axes, ("rppa", "procan", "APOLLO"), ("RPPA", "ProCan", "APOLLO"), ("a", "b", "c")):
        for method, label, color in methods:
            sub = cal[(cal.platform == platform) & (cal.method == method)].sort_values("K")
            ax.errorbar(
                sub["K"] + offsets[method],
                sub["cov80_mean"],
                yerr=[sub["cov80_mean"] - sub["cov80_p2_5"], sub["cov80_p97_5"] - sub["cov80_mean"]],
                fmt="o-",
                color=color,
                ms=3,
                lw=0.8,
                label=label,
            )
            lows.append(float(sub["cov80_p2_5"].min()))
            highs.append(float(sub["cov80_p97_5"].max()))
        ax.axhspan(0.75, 0.85, color="#D0D0D0", zorder=0)
        ax.axhline(0.80, color="black", ls="--", lw=0.6)
        ax.set_title(title, loc="center", fontsize=8)
        letter(ax, name)
        ax.set_xlabel("K")
        ax.set_xticks([10, 20, 50, 100])
        if platform == "APOLLO":
            ax.text(0.03, 0.06, "already within the band", transform=ax.transAxes, fontsize=7)
    pad = 0.03
    for ax in axes:
        ax.set_ylim(min(lows) - pad, max(highs) + pad)
    axes[0].set_ylabel("Coverage of nominal 80% intervals")
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "Fig3")


def draw_effects(ax, effects):
    y = np.arange(len(effects))[::-1]
    for i, rec in enumerate(effects.itertuples(index=False)):
        color = "#4C78A8" if rec.status == "사전등록" else "#B8B8B8"
        marker = "o" if rec.status == "사전등록" else "s"
        ax.plot(rec.mean, y[i], marker, color=color, ms=4)
        if np.isfinite(rec.low):
            ax.plot([rec.low, rec.high], [y[i], y[i]], color=color, lw=0.8)
    ax.axvline(0, color="black", lw=0.6)
    labels = []
    for rec in effects.itertuples(index=False):
        if rec.status == "사전등록":
            labels.append(f"{rec.item} {rec.cohort} ({rec.contrast})")
        else:
            labels.append(f"{rec.item} (post hoc)")
    ax.set_yticks(y, labels, fontsize=7)
    vals = [float(v) for v in effects["mean"]]
    for rec in effects.itertuples(index=False):
        if np.isfinite(rec.low):
            vals.extend([float(rec.low), float(rec.high)])
    lo, hi = min(vals), max(vals)
    pad = 0.08 * (hi - lo if hi > lo else 1.0)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xlabel("Mean difference")


def fig4():
    effects = pd.read_csv(SRC / "fig2_effects.tsv", sep="\t")
    scope = pd.read_csv(SRC / "fig4_scope.tsv", sep="\t")
    pre = effects[effects["status"] == "사전등록"]
    post_fx = effects[effects["status"] == "사후"]
    fig = plt.figure(figsize=(DOUBLE, 118 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.05], width_ratios=[1.25, 1], hspace=0.5, wspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    draw_effects(ax_a, pre)
    draw_effects(ax_b, post_fx)
    ax_a.set_title("Preregistered", loc="center", fontsize=8)
    ax_b.set_title("Post hoc; no intervals computed.", loc="center", fontsize=8)
    letter(ax_a, "a")
    letter(ax_b, "b")
    post = scope[scope["status"] == "사후"]
    pre_s = scope[scope["status"] == "사전등록"]
    ax = fig.add_subplot(gs[:, 1])
    ax.scatter(post["s"], post["median_delta"], c="#4C78A8", s=16, label="RPPA and ProCan")
    ax.scatter(pre_s["s"], pre_s["median_delta"], c="#E45756", marker="D", s=28, label="APOLLO")
    low = post.loc[post["s"].idxmin()]
    ax.annotate(
        "Haematopoietic and lymphoid",
        (low["s"], low["median_delta"]),
        textcoords="offset points",
        xytext=(18, 10),
        fontsize=7,
        arrowprops={"arrowstyle": "-", "lw": 0.4},
    )
    ax.axvline(0.8830509478307224, color="gray", lw=0.6, ls="--")
    ax.axvline(0.9788744329350166, color="gray", lw=0.6, ls="--")
    ax.axhline(0, color="black", lw=0.6)
    ax.text(0.03, 0.78, "Spearman ρ = 0.521", transform=ax.transAxes, fontsize=7)
    ax.set_xlabel("Similarity to CPTAC")
    ax.set_ylabel("Median ρ_AR − ρ_S0")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    letter(ax, "c")
    fig.subplots_adjust(left=0.32, right=0.98, bottom=0.12, top=0.94, hspace=0.48, wspace=0.32)
    save(fig, "Fig4")


def fig5():
    atlas = pd.read_csv(ROOT / "results" / "stage7" / "adc_atlas_final.tsv", sep="\t")
    pairs = pd.read_csv(SRC / "fig6_final_concordance.tsv", sep="\t")
    counts = pairs.groupby("source").size()
    print("fig5b source counts", counts.to_dict(), flush=True)
    fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 168 * MM), gridspec_kw={"width_ratios": [1.45, 1]})
    mat = atlas.pivot(index="gene_symbol", columns="disease", values="rho_hat")
    mat = mat.reindex(sorted(mat.index, key=lambda g: (g != "CLDN18", g)))
    plotted = mat.to_numpy(dtype=float)
    finite = plotted[np.isfinite(plotted)]
    vmin, vmax = float(np.min(finite)), float(np.max(finite))
    print("fig5a colour", vmin, vmax, "below_0.2", int(np.sum(finite < 0.2)), "above_0.9", int(np.sum(finite > 0.9)), flush=True)
    im = axes[0].imshow(plotted, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
    names = ["CLDN18*" if gene == "CLDN18" else gene for gene in mat.index]
    axes[0].set_yticks(range(len(names)), names, fontsize=7)
    for lab in axes[0].get_yticklabels():
        if lab.get_text() == "CLDN18*":
            lab.set_fontweight("bold")
    axes[0].set_xticks(range(len(mat.columns)), mat.columns, rotation=90, fontsize=7)
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.02, label="Forecast")
    colors = {"CPTAC": "#4C78A8", "RPPA": "#F58518", "ProCan_frozen_forecast": "#54A24B", "APOLLO": "#E45756"}
    labels = {"CPTAC": "CPTAC leave-one-cancer-out", "RPPA": "RPPA", "ProCan_frozen_forecast": "ProCan frozen forecast", "APOLLO": "APOLLO"}
    zorder = {"CPTAC": 2, "RPPA": 3, "ProCan_frozen_forecast": 4, "APOLLO": 5}
    ax = axes[1]
    for source, color in colors.items():
        sub = pairs[pairs["source"] == source]
        ax.scatter(sub["rho_hat"], sub["actual"], s=12, c=color, label=labels[source], alpha=0.75, zorder=zorder[source])
    cldn = pairs[(pairs.source == "APOLLO") & (pairs.gene_symbol == "CLDN18")]
    print("fig5b APOLLO", int((pairs.source == "APOLLO").sum()), "CLDN18", int(len(cldn)), flush=True)
    if len(cldn):
        x0 = float(cldn["rho_hat"].iloc[0])
        y0 = float(cldn["actual"].iloc[0])
        ax.scatter([x0], [y0], s=36, facecolors="none", edgecolors="black", linewidths=0.8, zorder=6)
        ax.annotate(
            "CLDN18",
            (x0, y0),
            textcoords="offset points",
            xytext=(14, -18),
            fontsize=7,
            arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "black"},
            zorder=7,
        )
    ax.plot([0, 1], [0, 1], color="gray", lw=0.6, zorder=1)
    ax.text(0.58, 0.66, "y = x", color="gray", fontsize=7, rotation=42)
    ax.set_xlabel("Forecast")
    ax.set_ylabel("Observed concordance")
    ax.legend(frameon=False, fontsize=7)
    letter(axes[0], "a")
    letter(axes[1], "b")
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.text(0.01, 0.008, "*Gene-level assays do not separate CLDN18.2.", fontsize=7)
    save(fig, "Fig5")


def fig_s1():
    diff = pd.read_csv(ROOT / "results" / "stage7" / "f0_diff.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(DOUBLE, 88 * MM))
    order = []
    for dataset in ("RPPA", "ProCan", "APOLLO"):
        for stratum in ("all", "in_P", "outside_P"):
            order.append((dataset, stratum))
    y = np.arange(len(order))[::-1]
    ax.axhspan(5.5, 8.5, color="#F2F2F2", zorder=0)
    ax.axhspan(-0.5, 2.5, color="#F2F2F2", zorder=0)
    ax.axhline(5.5, color="#B8B8B8", lw=0.6, zorder=0)
    ax.axhline(2.5, color="#B8B8B8", lw=0.6, zorder=0)
    for pred, color, shift in (("f1s_lgbm", "#4C78A8", 0.12), ("f1s_ridge", "#B279A2", -0.12)):
        for i, (dataset, stratum) in enumerate(order):
            rec = diff[(diff.dataset == dataset) & (diff.predictor == pred) & (diff.stratum == stratum)].iloc[0]
            ax.plot(rec["diff"], y[i] + shift, "o", color=color, ms=3.5, zorder=2)
            ax.plot([rec["ci_low"], rec["ci_high"]], [y[i] + shift, y[i] + shift], color=color, lw=0.7, zorder=2)
    ax.axvline(0, color="black", lw=0.6, zorder=1)
    ax.set_yticks(y, [f"{d} {s}" for d, s in order], fontsize=7)
    ax.set_xlabel("Forecast − CPTAC mean lookup")
    ax.plot([], [], "o", color="#4C78A8", label="F1s-LightGBM")
    ax.plot([], [], "o", color="#B279A2", label="F1s-ridge")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "FigS1")


def fig_s2():
    frame = pd.read_csv(ROOT / "results" / "stage8" / "within_context.tsv", sep="\t")
    sub = frame[frame.dataset == "RPPA"].copy()
    below = sub[sub["f1s_lgbm"] < sub["Lstar"]]
    above = sub[sub["f1s_lgbm"] > sub["Lstar"]]
    print("figs2 n", len(sub), "below", len(below), "above", len(above), below["context"].tolist(), flush=True)
    fig, ax = plt.subplots(figsize=(89 * MM, 89 * MM))
    ax.scatter(sub["Lstar"], sub["f1s_lgbm"], s=16, c="#4C78A8")
    lo = min(float(sub["Lstar"].min()), float(sub["f1s_lgbm"].min())) - 0.08
    hi = max(float(sub["Lstar"].max()), float(sub["f1s_lgbm"].max())) + 0.08
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.6)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    places = {
        "cholangiocarcinoma": (0.28, 0.30),
        "lung squamous cell carcinoma": (0.18, 0.58),
        "uterine corpus endometrioid carcinoma": (0.46, 0.24),
    }
    wrapped = {
        "cholangiocarcinoma": "cholangiocarcinoma",
        "lung squamous cell carcinoma": "lung squamous\ncell carcinoma",
        "uterine corpus endometrioid carcinoma": "uterine corpus\nendometrioid carcinoma",
    }
    for rec in below.itertuples(index=False):
        ax.annotate(
            wrapped[rec.context],
            (rec.Lstar, rec.f1s_lgbm),
            textcoords="data",
            xytext=places[rec.context],
            fontsize=7,
            arrowprops={"arrowstyle": "-", "lw": 0.4},
            annotation_clip=False,
        )
    ax.text(0.03, 0.97, "28 of 31 above the line", transform=ax.transAxes, va="top", fontsize=7)
    ax.set_xlabel("Cascade lookup")
    ax.set_ylabel("Forecast")
    fig.tight_layout()
    save(fig, "FigS2")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig_s1()
    fig_s2()
    out = ROOT / "results" / "stage17"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(FORMAT).to_csv(out / "figure_format.tsv", sep="\t", index=False)
    print("stage10 figures done", flush=True)
