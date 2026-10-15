"""Draft figures and tables for the APOLLO confirmation.

Source tables are written before the plots. Preregistered results and
post-hoc results are marked in the source files and in the plots.
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

font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42

OUT = ROOT / "results" / "figures"
SRC = OUT / "source"
CONFIRM = ROOT / "results" / "stage6" / "confirm"


def average_ranks(x):
    x = np.asarray(x, dtype=np.float64)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    n = int(x.size)
    distinct = np.empty(n, dtype=bool)
    distinct[0] = True
    if n > 1:
        distinct[1:] = xs[1:] != xs[:-1]
    groups = np.cumsum(distinct) - 1
    pos = np.arange(1, n + 1, dtype=np.float64)
    counts = np.bincount(groups)
    ranks = np.empty(n, dtype=np.float64)
    ranks[order] = (np.bincount(groups, weights=pos) / counts)[groups]
    return ranks


def finite_spearman(a, b, min_n=3):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    n = int(m.sum())
    if n < min_n:
        return np.nan, n
    x = a[m]
    y = b[m]
    if x.max() == x.min() or y.max() == y.min():
        return np.nan, n
    rx = average_ranks(x) - average_ranks(x).mean()
    ry = average_ranks(y) - average_ranks(y).mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if den == 0:
        return np.nan, n
    return float((rx * ry).sum() / den), n


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)


def write(name, frame):
    SRC.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SRC / name, sep="\t", index=False)


def fig1():
    pairs = pd.read_csv(ROOT / "results" / "ar_pairs.tsv", sep="\t", usecols=["heldout_cancer", "delta_ar"])
    rows = []
    for cancer, sub in pairs.groupby("heldout_cancer", sort=True):
        d = sub["delta_ar"].to_numpy(dtype=np.float64)
        d = d[np.isfinite(d)]
        rows.append({"status": "개발", "cancer": cancer, "n": int(d.size), "mean_delta": float(d.mean()), "median_delta": float(np.median(d))})
    delta = pd.DataFrame(rows)
    write("fig1_loco_delta.tsv", delta)
    pred = pd.read_csv(
        ROOT / "results" / "f1s_predictions.tsv",
        sep="\t",
        usecols=["heldout_cancer", "actual", "f1s_lgbm", "holdout_rna_sd_percentile"],
    )
    rel_rows = []
    for cancer, sub in pred.groupby("heldout_cancer", sort=True):
        for name, col in (("F1s-LightGBM", "f1s_lgbm"), ("B_sd", "holdout_rna_sd_percentile")):
            rho, n = finite_spearman(sub[col], sub["actual"], 3)
            rel_rows.append({"status": "개발", "cancer": cancer, "predictor": name, "n": n, "spearman": rho})
    rel = pd.DataFrame(rel_rows)
    write("fig1_loco_reliability.tsv", rel)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].bar(delta["cancer"], delta["mean_delta"], color="#4C78A8")
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_ylabel("평균 Δ (ρ_AR − ρ_S0)")
    axes[0].set_title("개발: CPTAC LOCO")
    axes[0].tick_params(axis="x", rotation=45)
    wide = rel.pivot(index="cancer", columns="predictor", values="spearman").reset_index()
    x = np.arange(len(wide))
    axes[1].bar(x - 0.18, wide["F1s-LightGBM"], width=0.36, label="F1s-LightGBM", color="#4C78A8")
    axes[1].bar(x + 0.18, wide["B_sd"], width=0.36, label="B_sd", color="#F58518")
    axes[1].set_xticks(x, wide["cancer"], rotation=45)
    axes[1].set_ylabel("Spearman(예측, 실제 ρ)")
    axes[1].set_title("개발: 보류 암종별")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    save(fig, "fig1")


def fig2():
    b4 = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "hypotheses.tsv", sep="\t")
    b6 = pd.read_csv(CONFIRM / "hypotheses.tsv", sep="\t")
    rows = []
    for name, cohort in (("H1", "TCGA RPPA"), ("H2", "ProCan")):
        rec = b4[b4["hypothesis"] == name].iloc[0]
        rows.append(
            {
                "status": "사전등록",
                "item": name,
                "cohort": cohort,
                "contrast": "AR−S0",
                "n": int(rec["n_units"]),
                "mean": float(rec["mean"]),
                "low": float(rec["mean_ci_low"]),
                "high": float(rec["mean_ci_high"]),
            }
        )
    for name, contrast in (("P1", "AR−S0"), ("Q1", "AR−R1"), ("Q2", "AR−R2")):
        rec = b6[b6["hypothesis"] == name].iloc[0]
        rows.append(
            {
                "status": "사전등록",
                "item": name,
                "cohort": "APOLLO",
                "contrast": contrast,
                "n": int(rec["n"]),
                "mean": float(rec["mean"]),
                "low": float(rec["mean_ci_low"]),
                "high": float(rec["mean_ci_high"]),
            }
        )
    post = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "posthoc_23_summary.tsv", sep="\t")
    for rec in post.itertuples(index=False):
        if rec.model == "AR":
            continue
        rows.append(
            {
                "status": "사후",
                "item": f"{rec.platform}-{rec.model}",
                "cohort": rec.platform,
                "contrast": rec.model,
                "n": int(rec.n_units),
                "mean": float(rec.mean_of_unit_means),
                "low": np.nan,
                "high": np.nan,
            }
        )
    frame = pd.DataFrame(rows)
    write("fig2_effects.tsv", frame)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    y = np.arange(len(frame))[::-1]
    for i, rec in enumerate(frame.itertuples(index=False)):
        color = "#4C78A8" if rec.status == "사전등록" else "#B8B8B8"
        marker = "o" if rec.status == "사전등록" else "s"
        yy = y[i]
        ax.plot(rec.mean, yy, marker, color=color)
        if np.isfinite(rec.low):
            ax.plot([rec.low, rec.high], [yy, yy], color=color)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y, [f"{r.item} {r.contrast} ({r.status})" for r in frame.itertuples(index=False)])
    ax.set_xlabel("단위 평균 Δ")
    fig.tight_layout()
    save(fig, "fig2")


def fig3():
    b4 = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "hypotheses.tsv", sep="\t")
    b6 = pd.read_csv(CONFIRM / "hypotheses.tsv", sep="\t")
    raw = pd.read_csv(ROOT / "results" / "stage6" / "posthoc" / "posthoc_05_raw.tsv", sep="\t")
    raw_map = dict(zip(raw["item"], raw["value"]))
    rows = []
    for name, cohort in (("H3", "TCGA RPPA"), ("H4", "ProCan")):
        rec = b4[b4["hypothesis"] == name].iloc[0]
        rows.append({"status": "사전등록", "item": name, "cohort": cohort, "which": "S_F", "value": float(rec["S_F"]), "low": float(rec["S_F_ci_low"]), "high": float(rec["S_F_ci_high"]), "n": int(rec["n_pairs"])})
        rows.append({"status": "사전등록", "item": name, "cohort": cohort, "which": "S_B", "value": float(rec["S_B"]), "low": float(rec["S_B_ci_low"]), "high": float(rec["S_B_ci_high"]), "n": int(rec["n_pairs"])})
    rec = b6[b6["hypothesis"] == "P2"].iloc[0]
    rows.append({"status": "사전등록", "item": "P2", "cohort": "APOLLO", "which": "S_F", "value": float(rec["S_F"]), "low": float(rec["S_F_ci_low"]), "high": float(rec["S_F_ci_high"]), "n": int(rec["n"])})
    rows.append({"status": "사전등록", "item": "P2", "cohort": "APOLLO", "which": "S_B", "value": float(rec["S_B"]), "low": float(rec["S_B_ci_low"]), "high": float(rec["S_B_ci_high"]), "n": int(rec["n"])})
    rows.append({"status": "사후", "item": "비대치", "cohort": "APOLLO", "which": "S_F", "value": float(raw_map["P2_S_F"]), "low": np.nan, "high": np.nan, "n": int(float(raw_map["P2_n"]))})
    rows.append({"status": "사후", "item": "비대치", "cohort": "APOLLO", "which": "S_B", "value": float(raw_map["P2_S_B"]), "low": np.nan, "high": np.nan, "n": int(float(raw_map["P2_n"]))})
    frame = pd.DataFrame(rows)
    write("fig3_reliability.tsv", frame)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    items = ["H3", "H4", "P2", "비대치"]
    for which, color, shift in (("S_F", "#4C78A8", -0.12), ("S_B", "#F58518", 0.12)):
        sub = frame[frame["which"] == which]
        for i, item in enumerate(items):
            rec = sub[sub["item"] == item].iloc[0]
            yy = len(items) - 1 - i + shift
            ax.plot(rec["value"], yy, "o" if rec["status"] == "사전등록" else "s", color=color, label=which if i == 0 else None)
            if np.isfinite(rec["low"]):
                ax.plot([rec["low"], rec["high"]], [yy, yy], color=color)
    ax.set_yticks(range(len(items))[::-1], [f"{item}" for item in items])
    ax.set_xlabel("Spearman")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "fig3")


def fig4():
    sim = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "posthoc_22_similarity.tsv", sep="\t")
    use = sim[sim["in_confirmation_34"] == 1][["cohort", "context", "s", "median_delta", "n_pairs"]].copy()
    use.insert(0, "status", "사후")
    ar = pd.read_csv(CONFIRM / "descriptive_ar.tsv", sep="\t")
    s_tab = pd.read_csv(CONFIRM / "descriptive_s.tsv", sep="\t")
    apollo = pd.DataFrame(
        [
            {
                "status": "사전등록",
                "cohort": "APOLLO",
                "context": "LUAD",
                "s": float(s_tab["value"].iloc[0]),
                "median_delta": float(ar.loc[ar["item"] == "median", "value"].iloc[0]),
                "n_pairs": int(ar.loc[ar["item"] == "n", "value"].iloc[0]),
            }
        ]
    )
    frame = pd.concat([use, apollo], ignore_index=True)
    write("fig4_scope.tsv", frame)
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    post = frame[frame["status"] == "사후"]
    pre = frame[frame["status"] == "사전등록"]
    ax.scatter(post["s"], post["median_delta"], c="#4C78A8", label=f"사후 맥락 n={len(post)}", s=28)
    ax.scatter(pre["s"], pre["median_delta"], c="#E45756", marker="D", s=50, label="사전등록 APOLLO")
    ax.axvline(0.8830509478307224, color="gray", lw=0.6, ls="--")
    ax.axvline(0.9788744329350166, color="gray", lw=0.6, ls="--")
    ax.set_xlabel("s")
    ax.set_ylabel("Δ 중앙값")
    ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "fig4")


def fig5():
    old = pd.read_csv(ROOT / "results" / "stage5" / "posthoc" / "posthoc_24_summary.tsv", sep="\t")
    old = old[old["predictor"] == "f1s_lgbm"][["platform", "K", "cov80_none", "cov80_bias", "cov80_both", "n_draws"]].copy()
    old.insert(0, "status", "사후")
    new = pd.read_csv(CONFIRM / "descriptive_calibration.tsv", sep="\t")
    new = new[new["predictor"] == "f1s_lgbm"][["K", "cov80_none", "cov80_bias", "cov80_both", "n_draws"]].copy()
    new.insert(0, "platform", "APOLLO")
    new.insert(0, "status", "사전등록")
    frame = pd.concat([old, new], ignore_index=True)
    write("fig5_calibration.tsv", frame)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    for ax, platform in zip(axes, ["rppa", "procan", "APOLLO"]):
        sub = frame[frame["platform"] == platform].sort_values("K")
        status = sub["status"].iloc[0]
        ax.plot(sub["K"], sub["cov80_none"], "o-", label="보정 없음")
        ax.plot(sub["K"], sub["cov80_bias"], "s-", label="δ만")
        ax.plot(sub["K"], sub["cov80_both"], "D-", label="δ와 분위수")
        ax.axhline(0.8, color="gray", lw=0.6, ls="--")
        ax.set_title(f"{platform} ({status})")
        ax.set_xlabel("K")
        ax.set_xticks([10, 20, 50, 100])
    axes[0].set_ylabel("80% 포함률")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, "fig5")


def fig6():
    atlas = pd.read_csv(ROOT / "results" / "stage5" / "atlas" / "adc_atlas_v2.tsv", sep="\t")
    heat = atlas[["gene_symbol", "disease", "rho_hat"]].copy()
    heat.insert(0, "status", "사후")
    write("fig6_heatmap.tsv", heat)
    pairs = pd.read_csv(ROOT / "results" / "stage5" / "atlas" / "concordance_pairs.tsv", sep="\t")
    pairs = pairs[pairs["source"].isin(["CPTAC", "RPPA", "ProCan"])][["source", "gene_symbol", "context", "rho_hat", "actual"]].copy()
    pairs.insert(0, "status", "사후")
    adc = pd.read_csv(CONFIRM / "descriptive_adc.tsv", sep="\t")
    extra = adc.rename(columns={"rho_hat": "rho_hat", "actual": "actual"})
    extra = extra[np.isfinite(extra["actual"]) & np.isfinite(extra["rho_hat"])][["gene_symbol", "actual", "rho_hat"]].copy()
    extra.insert(0, "context", "LUAD")
    extra.insert(0, "source", "APOLLO")
    extra.insert(0, "status", "사후")
    conc = pd.concat([pairs, extra], ignore_index=True)
    write("fig6_concordance.tsv", conc)
    mat = heat.pivot(index="gene_symbol", columns="disease", values="rho_hat")
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.2), gridspec_kw={"width_ratios": [1.4, 1]})
    im = axes[0].imshow(mat.to_numpy(), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    axes[0].set_yticks(range(len(mat.index)), mat.index, fontsize=7)
    axes[0].set_xticks(range(len(mat.columns)), mat.columns, rotation=90, fontsize=6)
    axes[0].set_title("사후: rho_hat")
    fig.colorbar(im, ax=axes[0], fraction=0.03)
    colors = {"CPTAC": "#4C78A8", "RPPA": "#F58518", "ProCan": "#54A24B", "APOLLO": "#E45756"}
    for source, sub in conc.groupby("source"):
        axes[1].scatter(sub["rho_hat"], sub["actual"], s=18 if source != "APOLLO" else 36, c=colors[source], label=f"{source} n={len(sub)}", alpha=0.75)
    axes[1].plot([-0.2, 1.2], [-0.2, 1.2], color="gray", lw=0.6)
    axes[1].set_xlabel("rho_hat")
    axes[1].set_ylabel("실제 ρ")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].set_title("사후: 일치")
    fig.tight_layout()
    save(fig, "fig6")


def tables():
    cptac = pd.read_csv(ROOT / "results" / "cptac_sample_counts.tsv", sep="\t")
    tumors = cptac[cptac["tissue"] == "tumor"]
    ccle = pd.read_csv(ROOT / "results" / "ccle_lineage_metrics.tsv", sep="\t")
    rppa = pd.read_csv(ROOT / "results" / "rppa_structure.tsv", sep="\t").iloc[0]
    pro = pd.read_csv(ROOT / "results" / "procan_structure.tsv", sep="\t").iloc[0]
    counts = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "counts.tsv", sep="\t")
    cmap = dict(zip(counts["item"], counts["value"]))
    rows = [
        {"dataset": "CPTAC", "role": "개발", "assay": "proteome and RNA, 10 tumor types", "n_samples": int(tumors["n_both"].sum()), "n_units": 7432, "unit": "P genes", "detail": "tumor RNA∩protein samples summed over 10 cancers"},
        {"dataset": "CCLE Lung", "role": "개발", "assay": "cell-line proteome and RNA", "n_samples": int(ccle.loc[ccle["lineage"] == "Lung", "n_cell_lines"].iloc[0]), "n_units": int(ccle.loc[ccle["lineage"] == "Lung", "n_s0"].iloc[0]), "unit": "genes with S0", "detail": ""},
        {"dataset": "CCLE Skin", "role": "개발", "assay": "cell-line proteome and RNA", "n_samples": int(ccle.loc[ccle["lineage"] == "Skin", "n_cell_lines"].iloc[0]), "n_units": int(ccle.loc[ccle["lineage"] == "Skin", "n_s0"].iloc[0]), "unit": "genes with S0", "detail": ""},
        {"dataset": "TCGA RPPA", "role": "외부", "assay": "RPPA antibodies", "n_samples": int(rppa["n_samples"]), "n_units": int(cmap["n_h1_units"]), "unit": "H1 antibodies in AR", "detail": f"confirmation pairs {cmap['n_rppa_pairs']}"},
        {"dataset": "ProCan", "role": "외부", "assay": "mass spectrometry protein matrix", "n_samples": int(pro["n_samples"]), "n_units": int(cmap["n_h2_units"]), "unit": "H2 proteins", "detail": f"protein columns {int(pro['n_protein_columns'])}; confirmation pairs {cmap['n_procan_pairs']}"},
        {"dataset": "APOLLO LUAD", "role": "외부", "assay": "TMT11plex protein and gene-level TPM RNA", "n_samples": 87, "n_units": 7029, "unit": "P1 proteins with defined d", "detail": "U_P1 upper bound 7148; U_P2 upper bound 7318; P2 n 7196"},
    ]
    write("table1_data.tsv", pd.DataFrame(rows))
    b4 = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "hypotheses.tsv", sep="\t")
    b4h = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "holm.tsv", sep="\t")
    b6 = pd.read_csv(CONFIRM / "hypotheses.tsv", sep="\t")
    b6h = pd.read_csv(CONFIRM / "holm.tsv", sep="\t")
    hyp = []
    for name in ("H1", "H2", "H3", "H4"):
        rec = b4[b4["hypothesis"] == name].iloc[0]
        call = b4h[b4h["hypothesis"] == name]["call"].iloc[0]
        hyp.append(
            {
                "hypothesis": name,
                "family": "B4",
                "n": int(rec["n_pairs"]) if name in ("H3", "H4") else int(rec["n_units"]),
                "estimate": float(rec["mean"]) if name in ("H1", "H2") else float(rec["diff"]),
                "pvalue": float(rec["pvalue"]),
                "call": call,
                "prereg_commit": "f7befd8324bcef63d89cebb57de36c0232b5b18f",
                "prereg_time": "2026-10-02 09:57:46 +0900",
                "output_commit": "a4764533d08bd09f404f827f7e01726c6a897cff",
                "output_time": "2026-10-02 10:29:53 +0900",
            }
        )
    for name in ("P1", "P2", "Q1", "Q2"):
        rec = b6[b6["hypothesis"] == name].iloc[0]
        call = b6h[b6h["hypothesis"] == name]["call"].iloc[0]
        estimate = float(rec["mean"]) if name != "P2" else float(rec["diff"])
        hyp.append(
            {
                "hypothesis": name,
                "family": "B6 primary" if name in ("P1", "P2") else "B6 secondary",
                "n": int(rec["n"]),
                "estimate": estimate,
                "pvalue": float(rec["pvalue"]),
                "call": call,
                "prereg_commit": "87967f098ea1434f39c3cd3109d2edeea1c1a0ef",
                "prereg_time": "2026-10-02 15:50:05 +0900",
                "output_commit": "57a3a3816bec154ae0df4f38dfdae5155eedc464",
                "output_time": "2026-10-02 16:59:47 +0900",
            }
        )
    table = pd.DataFrame(hyp)
    write("table2_hypotheses.tsv", table)
    fig, ax = plt.subplots(figsize=(11, 2.6))
    ax.axis("off")
    show = pd.DataFrame(rows)[["dataset", "role", "n_samples", "n_units", "unit"]]
    ax.table(cellText=show.values, colLabels=list(show.columns), loc="center")
    fig.tight_layout()
    save(fig, "table1")
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.axis("off")
    show = table[["hypothesis", "n", "estimate", "pvalue", "call", "output_time"]]
    ax.table(cellText=np.round(show[["n", "estimate", "pvalue"]].to_numpy(), 6) if False else show.astype(str).values, colLabels=list(show.columns), loc="center")
    fig.tight_layout()
    save(fig, "table2")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    tables()
    print("figures done", flush=True)


if __name__ == "__main__":
    main()
