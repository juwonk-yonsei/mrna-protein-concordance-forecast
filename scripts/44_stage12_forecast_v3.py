"""Add model_gene_features to the all-gene forecast. No model is fitted."""

import os

os.environ["B_CPU_CAP"] = "none"

import gzip
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

stage3 = importlib.import_module("scripts.08_stage3")

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
SRC = RES / "stage11" / "forecast_all_genes_v2.tsv.gz"
OUT = RES / "stage12" / "forecast_all_genes_v3.tsv.gz"
TOL = 1e-8


def gene_rows():
    genes = [g for g in (RES / "G_RNA_genes.txt").read_text().split() if g]
    p_genes = [g for g in (RES / "P_genes.txt").read_text().split() if g]
    rho = stage3.defined_rho(genes)
    adc = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    adc_loco = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    pset = set(p_genes)
    for rec in adc_loco.itertuples(index=False):
        if rec.gene_symbol in pset or rec.gene_symbol not in set(adc):
            continue
        if int(rec.n) >= 30 and np.isfinite(rec.s0):
            rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
    gmean = stage3.global_mean(stage3.CANCERS, rho, p_genes)
    stored = float(pd.read_csv(RES / "stage8" / "lookup_build.tsv", sep="\t")["global_mean"].iloc[0])
    empty = {}
    rows = {}
    for gene in genes:
        rec = stage3.f1s_feature_row(gene, stage3.CANCERS, rho, empty, empty, empty, gmean, np.nan, np.nan)
        rows[gene] = rec
    return rows, gmean, stored


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows, gmean, stored = gene_rows()
    flags = {}
    for gene, rec in rows.items():
        close = all(abs(float(rec[name]) - stored) <= TOL for name in ("rho_mean", "rho_min", "rho_max"))
        default = int(rec["n_defined_cancers"]) == 0 and float(rec["rho_sd"]) == 0.0 and close
        flags[gene] = 0 if default else 1
    n_mismatch = 0
    n_flag0 = 0
    n_flag0_pooled = 0
    seen = set()
    pred_max = 0.0
    with gzip.open(SRC, "rt") as src, gzip.open(OUT, "wt", newline="") as dst:
        header = src.readline().rstrip("\n")
        names = header.split("\t")
        dst.write(header + "\tmodel_gene_features\n")
        i_gene = names.index("gene_symbol")
        i_n = names.index("n_defined_cancers")
        i_rho = names.index("cptac_pooled_rho")
        i_forecast = names.index("forecast")
        for line in src:
            raw = line.rstrip("\n")
            parts = raw.split("\t")
            gene = parts[i_gene]
            rec = rows[gene]
            if int(parts[i_n]) != int(rec["n_defined_cancers"]):
                n_mismatch += 1
            flag = flags[gene]
            if gene not in seen:
                seen.add(gene)
                if flag == 0:
                    n_flag0 += 1
                    if parts[i_rho] != "":
                        n_flag0_pooled += 1
            dst.write(raw + "\t" + str(flag) + "\n")
            pred_max = max(pred_max, 0.0)
            _ = parts[i_forecast]
    if n_mismatch:
        raise RuntimeError(f"n_defined mismatch {n_mismatch}")
    # Byte check of the shared columns, and prediction max abs.
    max_abs = 0.0
    with gzip.open(SRC, "rt") as src, gzip.open(OUT, "rt") as dst:
        src.readline()
        dst.readline()
        i_forecast = names.index("forecast")
        for left, right in zip(src, dst):
            left = left.rstrip("\n")
            right = right.rstrip("\n")
            prefix, _flag = right.rsplit("\t", 1)
            if prefix != left:
                raise RuntimeError("shared columns are not byte-identical")
            a = float(left.split("\t")[i_forecast])
            b = float(prefix.split("\t")[i_forecast])
            max_abs = max(max_abs, abs(a - b))
    pd.DataFrame([{
        "n_rows": int(sum(1 for _ in gzip.open(OUT, "rt")) - 1),
        "n_genes": len(seen),
        "n_genes_model_gene_features_0": n_flag0,
        "n_flag0_with_cptac_pooled_rho": n_flag0_pooled,
        "feature_gmean": gmean,
        "stored_global_mean": stored,
        "gmean_abs_diff": abs(gmean - stored),
        "n_defined_mismatch": n_mismatch,
        "prediction_max_abs": max_abs,
        "bytes": OUT.stat().st_size,
    }]).to_csv(RES / "stage12" / "forecast_v3_check.tsv", sep="\t", index=False)
    print("genes", len(seen), "flag0", n_flag0, "flag0_pooled", n_flag0_pooled, "gmean_diff", abs(gmean - stored), flush=True)
    print("stage12 v3 done", flush=True)


if __name__ == "__main__":
    main()
