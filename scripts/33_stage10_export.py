"""Export the public F1s-LightGBM text model from a frozen-parameter refit.

Attempt 1 parameters only. The APOLLO comparison is the identity check.
Confirmation files and frozen predictions are not rewritten.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import importlib

import numpy as np
import pandas as pd
from lightgbm import Booster

stage9 = importlib.import_module("scripts.31_stage9_export")

ROOT = stage9.ROOT
OUT = ROOT / "results" / "stage10"
REL = ROOT / "release" / "models"
LIMIT = 1e-12


def filled_matrix(fit, frame):
    x = frame[stage9.F1S_FEATURES].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    med = np.nanmedian(fit["x_tr"], axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(med, inds[1])
    return filled, med


def max_abs(a, b):
    both = np.isfinite(a) & np.isfinite(b)
    if not both.any():
        return np.nan
    return float(np.max(np.abs(a[both] - b[both])))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    REL.mkdir(parents=True, exist_ok=True)
    params = stage9.frozen_params()
    print("frozen params", params, flush=True)
    _train, x_f, y = stage9.training_matrix()
    genes, p_genes, adc = stage9.stage3.load_lists()
    rho, mean_pct, sd_pct, rate, gmean = stage9.rho_tables(genes, p_genes, adc)
    cache_r = ROOT / "results" / "stage9" / "cache_rppa_features.tsv"
    cache_p = ROOT / "results" / "stage9" / "cache_procan_features.tsv"
    print("loading cached features", flush=True)
    rppa = pd.read_csv(cache_r, sep="\t")
    pro = pd.read_csv(cache_p, sep="\t")
    print("fit", flush=True)
    fit = stage9.fit_lgbm(params, x_f, y)
    print("APOLLO features", flush=True)
    u_p2 = pd.read_csv(ROOT / "results" / "stage6" / "prep" / "U_P2_genes.tsv", sep="\t")["gene_symbol"].tolist()
    mat = pd.read_parquet(ROOT / "data" / "processed" / "apollo_luad" / "rna_log2tpm.parquet")
    samples = mat.columns.astype(str).tolist()
    table = stage9.stage3_reliability(genes, p_genes, u_p2, mat, samples, fit)
    apollo_old = pd.read_csv(ROOT / "results" / "stage6" / "frozen" / "apollo_gene_pred.tsv", sep="\t")
    merged = table.merge(apollo_old, on="gene_symbol", suffixes=("_new", "_old"))
    apollo_gap = max_abs(merged["f1s_lgbm_new"].to_numpy(dtype=np.float64), merged["f1s_lgbm_old"].to_numpy(dtype=np.float64))
    print("apollo max abs", apollo_gap, "n", len(merged), flush=True)
    if not np.isfinite(apollo_gap) or apollo_gap > LIMIT:
        pd.DataFrame([{"apollo_max_abs": apollo_gap, "n": len(merged), "exported": 0}]).to_csv(OUT / "export_check.tsv", sep="\t", index=False)
        raise SystemExit("APOLLO identity check failed")
    text_path = REL / "F1s_lightgbm.txt"
    fit["lgbm"].booster_.save_model(str(text_path))
    _filled_train, med = filled_matrix(fit, rppa.iloc[:1])
    import lightgbm
    meta = {
        "features": list(stage9.F1S_FEATURES),
        "median_fill": [float(v) for v in med],
        "lightgbm_version": lightgbm.__version__,
    }
    json_path = REL / "F1s_lightgbm_features.json"
    json_path.write_text(json.dumps(meta, indent=2) + "\n")
    booster = Booster(model_file=str(text_path))
    rows = []
    for name, frame in (("APOLLO", table), ("RPPA", rppa), ("ProCan", pro)):
        filled, _med = filled_matrix(fit, frame)
        a = fit["lgbm"].predict(filled)
        b = booster.predict(filled)
        gap = max_abs(a, b)
        rows.append({"dataset": name, "reload_max_abs": gap, "n": int(len(frame)), "apollo_identity_max_abs": apollo_gap})
        print("reload", name, gap, flush=True)
        if not np.isfinite(gap) or gap > LIMIT:
            pd.DataFrame(rows).to_csv(OUT / "export_check.tsv", sep="\t", index=False)
            raise SystemExit(f"reload check failed for {name}")
    pd.DataFrame(rows).to_csv(OUT / "export_check.tsv", sep="\t", index=False)
    print("stage10 export done", flush=True)


if __name__ == "__main__":
    main()
