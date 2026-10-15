"""Refit F1s-LightGBM on the training rows stored at the B4 freeze commit.

The rows are written only under /tmp/b11_b4. The repository tree is not changed
by the extraction. num_threads is not set. A model file is written only if the
maximum absolute difference from the frozen predictions is at most 1e-12.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import hashlib
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import importlib

import numpy as np
import pandas as pd
import yaml
from lightgbm import Booster, LGBMRegressor
from scipy.stats import rankdata

stage3 = importlib.import_module("scripts.08_stage3")
lookup = importlib.import_module("scripts.25_stage8_lookup")

ROOT = stage3.ROOT
OUT = ROOT / "results" / "stage11"
TMP = Path("/tmp/b11_b4")
COMMIT = "224779714bffd58b02059f51de48d92d1324a24f"
LIMIT = 1e-12
FEATURES = stage3.F1S_FEATURES


def frozen_params():
    b2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
    stage3.B2 = b2
    stage3.B2_LGBM = dict(b2["reliability"]["lightgbm"])
    return stage3.lgbm_params()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def training_matrix(path):
    preds = pd.read_csv(path, sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    train = preds[preds["in_P"] & np.isfinite(preds["actual"])].copy()
    x = train[FEATURES].to_numpy(dtype=np.float64)
    y = train["actual"].to_numpy(dtype=np.float64)
    filled, _, _ = stage3.fill_median(x, x[:1])
    return train, filled, y


def predict(model, x_train, frame):
    x = frame[FEATURES].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    med = np.nanmedian(x_train, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(med, inds[1])
    return model.predict(filled)


def compare(name, new, old, keys):
    merged = new.merge(old[keys + ["f1s_lgbm"]], on=keys, suffixes=("_new", "_old"))
    if len(merged) != len(old) or len(merged) != len(new):
        raise SystemExit(f"{name} alignment {len(new)} {len(old)} {len(merged)}")
    a = merged["f1s_lgbm_new"].to_numpy(dtype=np.float64)
    b = merged["f1s_lgbm_old"].to_numpy(dtype=np.float64)
    both = np.isfinite(a) & np.isfinite(b)
    diff = np.abs(a[both] - b[both])
    ra = rankdata(a[both], method="average")
    rb = rankdata(b[both], method="average")
    return {
        "dataset": name,
        "n_aligned": int(len(merged)),
        "n_finite_pred": int(both.sum()),
        "max_abs": float(diff.max()) if diff.size else np.nan,
        "spearman_new_vs_frozen": lookup.finite_spearman(a, b, 3),
        "n_rank_changed": int(np.sum(ra != rb)),
        "limit": LIMIT,
        "passed": int(diff.size and float(diff.max()) <= LIMIT),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    dest = TMP / "f1s_predictions.tsv"
    blob = subprocess.check_output(["git", "show", f"{COMMIT}:results/f1s_predictions.tsv"], cwd=ROOT)
    dest.write_bytes(blob)
    current = sha256(ROOT / "results" / "f1s_predictions.tsv")
    extracted = sha256(dest)
    params = frozen_params()
    if "num_threads" in params:
        raise SystemExit("frozen params contain num_threads")
    print("params", params, flush=True)
    print("sha256 extracted", extracted, "current", current, "equal", int(extracted == current), flush=True)
    _train, x_f, y = training_matrix(dest)
    print("training rows", len(y), flush=True)
    model = LGBMRegressor(**params)
    model.fit(x_f, y)
    rppa_x = pd.read_csv(ROOT / "results" / "stage9" / "cache_rppa_features.tsv", sep="\t")
    pro_x = pd.read_csv(ROOT / "results" / "stage9" / "cache_procan_features.tsv", sep="\t")
    rppa_x = rppa_x.copy()
    pro_x = pro_x.copy()
    rppa_x["f1s_lgbm"] = predict(model, x_f, rppa_x)
    pro_x["f1s_lgbm"] = predict(model, x_f, pro_x)
    rppa_old = pd.read_csv(ROOT / "results" / "frozen_forecasts_rppa_F1s.tsv", sep="\t")
    pro_old = pd.read_csv(ROOT / "results" / "frozen_forecasts_procan_F1s.tsv", sep="\t")
    rows = [
        compare("RPPA", rppa_x, rppa_old, ["xena_antibody", "gene_symbol", "cancer"]),
        compare("ProCan", pro_x, pro_old, ["protein_id", "gene_symbol", "lineage"]),
    ]
    frame = pd.DataFrame(rows)
    frame.insert(0, "training_commit", COMMIT)
    frame.insert(1, "training_sha256", extracted)
    frame.insert(2, "current_training_sha256", current)
    frame.insert(3, "training_bytes_equal", int(extracted == current))
    frame.insert(4, "n_training_rows", int(len(y)))
    frame.to_csv(OUT / "b4_refit_compare.tsv", sep="\t", index=False)
    print(frame.to_string(index=False), flush=True)
    if int(frame["passed"].min()) == 1:
        rel = ROOT / "release" / "models"
        rel.mkdir(parents=True, exist_ok=True)
        model.booster_.save_model(str(rel / "F1s_lightgbm_B4.txt"))
        reloaded = Booster(model_file=str(rel / "F1s_lightgbm_B4.txt"))
        reload_rows = []
        for name, features in (("RPPA", rppa_x), ("ProCan", pro_x)):
            x = features[FEATURES].to_numpy(dtype=np.float64)
            filled = np.array(x, dtype=np.float64, copy=True)
            med = np.nanmedian(x_f, axis=0)
            med = np.where(np.isfinite(med), med, 0.0)
            inds = np.where(~np.isfinite(filled))
            if inds[0].size:
                filled[inds] = np.take(med, inds[1])
            hat = reloaded.predict(filled)
            gap = np.max(np.abs(hat - features["f1s_lgbm"].to_numpy(dtype=np.float64)))
            reload_rows.append({"dataset": name, "reload_max_abs": float(gap)})
        pd.DataFrame(reload_rows).to_csv(OUT / "b4_reload_check.tsv", sep="\t", index=False)
        print(reload_rows, flush=True)
    else:
        print("limit not met; model file not written", flush=True)
    print("stage11 b4 refit done", flush=True)


if __name__ == "__main__":
    main()
