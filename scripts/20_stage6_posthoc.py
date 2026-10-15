"""Post-hoc tables after the APOLLO confirmation output was committed.

Every row is labeled 사후. Does not rewrite the preregistered confirmation
directory or adc_atlas.tsv.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.common import ROOT

LABEL = "사후"
OUT = ROOT / "results" / "stage6" / "posthoc"
CONFIRM = ROOT / "results" / "stage6" / "confirm"
ATLAS = ROOT / "results" / "stage5" / "atlas" / "adc_atlas.tsv"
ATLAS_V2 = ROOT / "results" / "stage5" / "atlas" / "adc_atlas_v2.tsv"


def average_ranks(x):
    x = np.asarray(x, dtype=np.float64)
    n = int(x.size)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
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
    if int(m.sum()) < min_n:
        return np.nan, int(m.sum())
    x = a[m]
    y = b[m]
    if x.max() == x.min() or y.max() == y.min():
        return np.nan, int(m.sum())
    rx = average_ranks(x)
    ry = average_ranks(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    if den == 0:
        return np.nan, int(m.sum())
    return float((rx * ry).sum() / den), int(m.sum())


def write(path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def extremes():
    genes = pd.read_csv(CONFIRM / "genes.tsv", sep="\t")
    pairs = pd.read_csv(ROOT / "results" / "ar_pairs.tsv", sep="\t", usecols=["heldout_cancer", "gene_symbol", "delta_ar"])
    loco = pairs[pairs["heldout_cancer"] == "LUAD"][["gene_symbol", "delta_ar"]]
    use = genes[np.isfinite(genes["d_p1"])][["gene_symbol", "d_p1", "actual", "rho_ar", "rho_s0"]].copy()
    use = use.merge(loco, on="gene_symbol", how="left")
    use = use.sort_values(["d_p1", "gene_symbol"], ascending=[False, True])
    top = use.head(15).copy()
    top["end"] = "top"
    bottom = use.tail(15).copy()
    bottom["end"] = "bottom"
    both = pd.concat([top, bottom], ignore_index=True)
    both.insert(0, "label", LABEL)
    write(OUT / "posthoc_01_extremes.tsv", both)
    scatter = use.rename(columns={"d_p1": "apollo_d", "delta_ar": "cptac_luad_d"})
    scatter.insert(0, "label", LABEL)
    rho, n = finite_spearman(scatter["apollo_d"], scatter["cptac_luad_d"], 3)
    write(OUT / "posthoc_01_scatter.tsv", scatter)
    write(
        OUT / "posthoc_01_summary.tsv",
        pd.DataFrame(
            [
                {
                    "label": LABEL,
                    "n_apollo_d": int(use.shape[0]),
                    "n_both": int(np.isfinite(use["delta_ar"]).sum()),
                    "spearman_d": rho,
                    "n_spearman": n,
                }
            ]
        ),
    )


def bias_row(name, actual, hat):
    actual = np.asarray(actual, dtype=np.float64)
    hat = np.asarray(hat, dtype=np.float64)
    m = np.isfinite(actual) & np.isfinite(hat)
    err = actual[m] - hat[m]
    return {
        "label": LABEL,
        "cohort": name,
        "n": int(m.sum()),
        "mean_actual_minus_f1s_lgbm": float(err.mean()) if m.any() else np.nan,
    }


def bias():
    f1s = pd.read_csv(ROOT / "results" / "f1s_predictions.tsv", sep="\t", usecols=["heldout_cancer", "gene_symbol", "actual", "f1s_lgbm"])
    cptac = f1s[f1s["heldout_cancer"] == "LUAD"]
    genes = pd.read_csv(CONFIRM / "genes.tsv", sep="\t", usecols=["gene_symbol", "in_U_P2", "actual"])
    level = pd.read_csv(CONFIRM / "gene_level.tsv", sep="\t")
    apollo = genes[(genes["in_U_P2"] == 1)].merge(level, on="gene_symbol", how="left")
    rppa = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "rppa_pairs.tsv", sep="\t", usecols=["gene_symbol", "cancer", "actual", "f1s_lgbm"])
    rppa = rppa[rppa["cancer"] == "lung adenocarcinoma"]
    pro = pd.read_csv(ROOT / "results" / "stage4" / "confirm" / "procan_pairs.tsv", sep="\t", usecols=["gene_symbol", "lineage", "actual", "f1s_lgbm"])
    pro = pro[pro["lineage"] == "Lung"]
    rows = [
        bias_row("CPTAC_LUAD_holdout", cptac["actual"], cptac["f1s_lgbm"]),
        bias_row("APOLLO", apollo["actual"], apollo["f1s_lgbm"]),
        bias_row("TCGA_RPPA_LUAD", rppa["actual"], rppa["f1s_lgbm"]),
        bias_row("ProCan_Lung", pro["actual"], pro["f1s_lgbm"]),
    ]
    write(OUT / "posthoc_02_bias.tsv", pd.DataFrame(rows))


def common():
    genes = pd.read_csv(CONFIRM / "genes.tsv", sep="\t", usecols=["gene_symbol", "actual", "d_p1"])
    rppa = pd.read_csv(
        ROOT / "results" / "stage4" / "confirm" / "rppa_pairs.tsv",
        sep="\t",
        usecols=["xena_antibody", "gene_symbol", "cancer", "actual", "delta"],
    )
    rppa = rppa[rppa["cancer"] == "lung adenocarcinoma"]
    pro = pd.read_csv(
        ROOT / "results" / "stage4" / "confirm" / "procan_pairs.tsv",
        sep="\t",
        usecols=["protein_id", "gene_symbol", "lineage", "actual", "delta"],
    )
    pro = pro[pro["lineage"] == "Lung"]
    r_counts = rppa.groupby("gene_symbol").size()
    p_counts = pro.groupby("gene_symbol").size()
    shared = sorted(set(r_counts.index) & set(p_counts.index) & set(genes["gene_symbol"]))
    gmap = genes.set_index("gene_symbol")
    rows = []
    n_single = 0
    for gene in shared:
        n_r = int(r_counts.loc[gene])
        n_p = int(p_counts.loc[gene])
        row = {
            "label": LABEL,
            "gene_symbol": gene,
            "n_rppa": n_r,
            "n_procan": n_p,
            "apollo_actual": float(gmap.loc[gene, "actual"]) if np.isfinite(gmap.loc[gene, "actual"]) else np.nan,
            "apollo_d": float(gmap.loc[gene, "d_p1"]) if np.isfinite(gmap.loc[gene, "d_p1"]) else np.nan,
            "rppa_actual": np.nan,
            "rppa_d": np.nan,
            "procan_actual": np.nan,
            "procan_d": np.nan,
        }
        if n_r == 1 and n_p == 1:
            n_single += 1
            rr = rppa[rppa["gene_symbol"] == gene].iloc[0]
            pp = pro[pro["gene_symbol"] == gene].iloc[0]
            row["rppa_actual"] = float(rr["actual"]) if np.isfinite(rr["actual"]) else np.nan
            row["rppa_d"] = float(rr["delta"]) if np.isfinite(rr["delta"]) else np.nan
            row["procan_actual"] = float(pp["actual"]) if np.isfinite(pp["actual"]) else np.nan
            row["procan_d"] = float(pp["delta"]) if np.isfinite(pp["delta"]) else np.nan
        rows.append(row)
    write(OUT / "posthoc_03_common.tsv", pd.DataFrame(rows))
    write(
        OUT / "posthoc_03_note.tsv",
        pd.DataFrame(
            [
                {
                    "label": LABEL,
                    "n_genes_in_all_three": len(shared),
                    "n_genes_one_row_each": n_single,
                    "rule": "Values are filled only when RPPA and ProCan each have one row for the gene. Multiple rows are not averaged.",
                }
            ]
        ),
    )


def atlas_v2():
    before = ATLAS.read_bytes()
    adc = pd.read_csv(CONFIRM / "descriptive_adc.tsv", sep="\t")
    values = {}
    for rec in adc.itertuples(index=False):
        values[rec.gene_symbol] = (rec.actual, rec.d, rec.rho_hat)
    lines = ATLAS.read_text().splitlines()
    header = lines[0].split("\t")
    gene_i = header.index("gene_symbol")
    disease_i = header.index("disease")
    out_lines = ["\t".join(header + ["apollo_actual", "apollo_d", "apollo_rho_hat"])]
    for line in lines[1:]:
        fields = line.split("\t")
        extra = ["", "", ""]
        if fields[disease_i] == "lung adenocarcinoma":
            actual, d_value, hat = values.get(fields[gene_i], (np.nan, np.nan, np.nan))
            extra = [("" if not np.isfinite(v) else repr(float(v))) for v in (actual, d_value, hat)]
        out_lines.append(line + "\t" + "\t".join(extra))
    ATLAS_V2.write_text("\n".join(out_lines) + "\n")
    if ATLAS.read_bytes() != before:
        raise RuntimeError("adc_atlas.tsv changed")
    pairs = pd.read_csv(ROOT / "results" / "stage5" / "atlas" / "concordance_pairs.tsv", sep="\t")
    rows = []
    for source, sub in pairs.groupby("source", sort=True):
        rho, n = finite_spearman(sub["rho_hat"], sub["actual"], 3)
        rows.append({"label": LABEL, "source": source, "n_pairs": n, "spearman_rho_hat_actual": rho})
    rho, n = finite_spearman(adc["rho_hat"], adc["actual"], 3)
    rows.append({"label": LABEL, "source": "APOLLO", "n_pairs": n, "spearman_rho_hat_actual": rho})
    write(OUT / "posthoc_04_concordance.tsv", pd.DataFrame(rows))


def main():
    extremes()
    bias()
    common()
    atlas_v2()
    print("posthoc 1-4 done", flush=True)


if __name__ == "__main__":
    main()
