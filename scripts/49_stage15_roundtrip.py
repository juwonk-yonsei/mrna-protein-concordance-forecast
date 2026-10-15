"""CSV round-trip check for the frozen F1s training table.

The criterion is config/analysis_plan_B15.yaml, committed before this run.
A fit happens only when that file's overall rule says float32_origin.
num_threads is not set. B_CPU_CAP is set before any import.
"""

import os

os.environ["B_CPU_CAP"] = "none"

import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage15"
TMP = Path("/tmp/b15")
COMMIT = "224779714bffd58b02059f51de48d92d1324a24f"
BLOB = "80d0d9935a037c4949760aa3e01ef5c9020009b4"
FEATURES = [
    "rho_mean",
    "rho_sd",
    "rho_min",
    "rho_max",
    "n_defined_cancers",
    "holdout_rna_mean_percentile",
    "holdout_rna_sd_percentile",
    "train_rna_mean_percentile_mean",
    "train_protein_rate_mean",
]
INT_RE = re.compile(r"^[+-]?\d+(?:\.0+)?(?:[eE][+-]?\d+)?$")
MISSING = {"", "na", "nan", "none", "null"}


def sig_digits(token):
    text = token.strip()
    if text.lower() in MISSING:
        return None
    body = text[1:] if text[0] in "+-" else text
    mantissa = re.split(r"[eE]", body, maxsplit=1)[0]
    if re.fullmatch(r"\d+\.\d*|\.\d+|\d+", mantissa or "") is None:
        return None
    digits = mantissa.replace(".", "").lstrip("0")
    return 1 if digits == "" else len(digits)


def exact_integer_text(token):
    text = token.strip()
    if INT_RE.fullmatch(text) is None:
        return False
    try:
        value = float(text)
    except ValueError:
        return False
    return bool(np.isfinite(value) and float(value).is_integer())


def check_counter():
    assert sig_digits("0.20749223") == 8
    assert sig_digits("1.2300") == 5
    assert sig_digits("-0.00123") == 3
    assert sig_digits("1.23e-4") == 3
    assert sig_digits("0") == 1
    assert sig_digits("100") == 3
    assert sig_digits("100.0") == 4
    assert sig_digits("nan") is None
    assert exact_integer_text("9")
    assert exact_integer_text("9.0")
    assert exact_integer_text("-3")
    assert exact_integer_text("1e1")
    assert not exact_integer_text("12e-1")
    assert not exact_integer_text("9.5")
    assert not exact_integer_text("0.20749223")


def extract():
    TMP.mkdir(parents=True, exist_ok=True)
    found = subprocess.check_output(
        ["git", "rev-parse", f"{COMMIT}:results/f1s_predictions.tsv"], cwd=ROOT, text=True
    ).strip()
    if found != BLOB:
        raise SystemExit(f"blob {found}")
    dest = TMP / "f1s_predictions.tsv"
    dest.write_bytes(subprocess.check_output(["git", "show", f"{COMMIT}:results/f1s_predictions.tsv"], cwd=ROOT))
    return dest


def text_stats(path):
    counts = {name: [] for name in FEATURES}
    non_integer = {name: [] for name in FEATURES}
    n_integer = {name: 0 for name in FEATURES}
    n_rows = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").rstrip("\r").split("\t")
        index = {name: header.index(name) for name in FEATURES}
        width = len(header)
        for line in handle:
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if len(parts) == 1 and parts[0] == "":
                continue
            if len(parts) != width:
                raise SystemExit(f"width {len(parts)} != {width}")
            n_rows += 1
            for name, pos in index.items():
                token = parts[pos]
                digits = sig_digits(token)
                if digits is None:
                    continue
                counts[name].append(digits)
                if exact_integer_text(token):
                    n_integer[name] += 1
                else:
                    non_integer[name].append(digits)
    return n_rows, counts, non_integer, n_integer


def summarize(values):
    if not values:
        return np.nan, np.nan, np.nan
    arr = np.asarray(values, dtype=np.float64)
    return int(arr.min()), float(np.median(arr)), int(arr.max())


def column_class(sig_max_noninteger, n_noninteger, float32_fraction, repr_fraction):
    if n_noninteger == 0:
        return "exact_integer"
    if sig_max_noninteger < 17 and float32_fraction < 1:
        return "digit_truncation"
    if sig_max_noninteger < 17 and float32_fraction == 1:
        return "float32_origin"
    if sig_max_noninteger >= 17 and repr_fraction == 1:
        return "no_loss"
    return "unrecovered"


def overall_form(classes):
    if any(name in ("digit_truncation", "unrecovered") for name in classes):
        return "digit_truncation_or_unrecovered"
    if any(name == "float32_origin" for name in classes):
        return "float32_origin"
    return "no_loss"


def file_check(path):
    n_rows, counts, non_integer, n_integer = text_stats(path)
    frame = pd.read_csv(path, sep="\t")
    if len(frame) != n_rows:
        raise SystemExit(f"row count text {n_rows} pandas {len(frame)}")
    rows = []
    for name in FEATURES:
        parsed = frame[name].to_numpy(dtype=np.float64)
        finite = np.isfinite(parsed)
        n_finite = int(finite.sum())
        n_nonfinite = int((~finite).sum())
        if n_finite != len(counts[name]):
            raise SystemExit(f"{name} finite text {len(counts[name])} parsed {n_finite}")
        values = parsed[finite]
        float32_equal = np.float32(values).astype(np.float64) == values
        repr_equal = np.fromiter((float(repr(float(v))) == float(v) for v in values), dtype=bool, count=n_finite)
        sig_min, sig_median, sig_max = summarize(counts[name])
        ni_min, ni_median, ni_max = summarize(non_integer[name])
        float32_fraction = float(float32_equal.mean()) if n_finite else np.nan
        repr_fraction = float(repr_equal.mean()) if n_finite else np.nan
        kind = column_class(
            -1 if not non_integer[name] else int(ni_max),
            len(non_integer[name]),
            float32_fraction,
            repr_fraction,
        )
        rows.append({
            "feature": name,
            "n_rows": n_rows,
            "n_finite": n_finite,
            "n_nonfinite": n_nonfinite,
            "n_exact_integer_text": n_integer[name],
            "sig_min": sig_min,
            "sig_median": sig_median,
            "sig_max": sig_max,
            "sig_min_noninteger": ni_min,
            "sig_median_noninteger": ni_median,
            "sig_max_noninteger": ni_max,
            "float32_equal_fraction": float32_fraction,
            "repr_equal_fraction": repr_fraction,
            "dtype": str(frame[name].dtype),
            "column_class": kind,
        })
    table = pd.DataFrame(rows)
    form = overall_form(table["column_class"].tolist())
    table["overall_form"] = form
    table["action"] = "fit_float32_roundtrip" if form == "float32_origin" else "no_fit_option_Z"
    return table


def fit_recovered(path, table):
    import importlib

    import yaml
    from lightgbm import Booster, LGBMRegressor
    from scipy.stats import rankdata, spearmanr

    sys.path.insert(0, str(ROOT))
    stage3 = importlib.import_module("scripts.08_stage3")
    b2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
    stage3.B2 = b2
    stage3.B2_LGBM = dict(b2["reliability"]["lightgbm"])
    params = stage3.lgbm_params()
    if "num_threads" in params or "force_col_wise" in params or "force_row_wise" in params:
        raise SystemExit(f"forbidden arg in params {sorted(params)}")
    preds = pd.read_csv(path, sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    train = preds[preds["in_P"].astype(bool) & np.isfinite(preds["actual"])].copy()
    raw = train[FEATURES].to_numpy(dtype=np.float64)
    recovered = raw.astype(np.float32).astype(np.float64)
    changed = np.isfinite(raw) & (recovered != raw)
    y = train["actual"].to_numpy(dtype=np.float64)
    x_f, _, _ = stage3.fill_median(recovered, recovered[:1])
    rppa = pd.read_csv(ROOT / "results" / "stage9" / "cache_rppa_features.tsv", sep="\t")
    pro = pd.read_csv(ROOT / "results" / "stage9" / "cache_procan_features.tsv", sep="\t")
    apollo = apollo_features(stage3)
    frozen = {
        "RPPA": (rppa, pd.read_csv(ROOT / "results" / "frozen_forecasts_rppa_F1s.tsv", sep="\t"), ["xena_antibody", "cancer"]),
        "ProCan": (pro, pd.read_csv(ROOT / "results" / "frozen_forecasts_procan_F1s.tsv", sep="\t"), ["protein_id", "lineage"]),
        "APOLLO": (apollo, pd.read_csv(ROOT / "results" / "stage6" / "frozen" / "apollo_gene_pred.tsv", sep="\t"), ["gene_symbol"]),
    }
    saved = {}
    models = []
    records = []
    for repeat in (1, 2):
        model = LGBMRegressor(**params)
        model.fit(x_f, y)
        models.append(model)
        hats = {}
        for dataset, (frame, old, keys) in frozen.items():
            hat = predict(model, x_f, frame)
            hats[dataset] = hat
            rec = compare(frame.assign(f1s_lgbm=hat), old, keys, rankdata, spearmanr)
            rec.update({"dataset": dataset, "repeat": repeat})
            records.append(rec)
        saved[repeat] = hats
    out = pd.DataFrame(records)
    pair = []
    for dataset in ("RPPA", "ProCan", "APOLLO"):
        gap = float(np.max(np.abs(saved[1][dataset] - saved[2][dataset])))
        pair.append({"dataset": dataset, "repeat_pair_max_abs": gap})
    out = out.merge(pd.DataFrame(pair), on="dataset")
    out["n_train"] = int(len(train))
    out["n_feature_values_changed_by_cast"] = int(changed.sum())
    both = np.isfinite(raw) & np.isfinite(recovered)
    out["max_abs_cast"] = float(np.max(np.abs(raw[both] - recovered[both]))) if both.any() else np.nan
    out["limit"] = 1e-12
    limit = 1e-12
    ok = {}
    for dataset in ("RPPA", "ProCan", "APOLLO"):
        gaps = out.loc[out["dataset"] == dataset, "max_abs"].tolist()
        ok[dataset] = gaps
    y_pass = (
        all(v <= limit for v in ok["RPPA"])
        and all(v <= limit for v in ok["ProCan"])
        and all(v > limit for v in ok["APOLLO"])
    )
    out["option_Y"] = int(y_pass)
    if y_pass:
        export(models[0], x_f, frozen, params, int(len(train)), saved[1])
    return out


def predict(model, x_train, frame):
    x = frame[FEATURES].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    med = np.nanmedian(x_train, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(med, inds[1])
    return model.predict(filled)


def compare(new, old, keys, rankdata, spearmanr):
    merged = new.merge(old[keys + ["f1s_lgbm"]], on=keys, suffixes=("_new", "_old"))
    if len(merged) != len(new) or len(merged) != len(old):
        raise SystemExit(f"alignment {len(merged)} {len(new)} {len(old)}")
    a = merged["f1s_lgbm_new"].to_numpy(dtype=np.float64)
    b = merged["f1s_lgbm_old"].to_numpy(dtype=np.float64)
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise SystemExit("non-finite prediction")
    return {
        "n_aligned": int(len(merged)),
        "max_abs": float(np.max(np.abs(a - b))),
        "spearman": float(spearmanr(a, b).statistic),
        "n_rank_changed": int(np.sum(rankdata(a, method="average") != rankdata(b, method="average"))),
    }


def apollo_features(stage3):
    genes, p_genes, adc = stage3.load_lists()
    targets = pd.read_csv(ROOT / "results" / "stage6" / "prep" / "U_P2_genes.tsv", sep="\t")["gene_symbol"].tolist()
    mat = pd.read_parquet(ROOT / "data" / "processed" / "apollo_luad" / "rna_log2tpm.parquet")
    samples = mat.columns.astype(str).tolist()
    rho = stage3.defined_rho(p_genes)
    adc_loco = pd.read_csv(ROOT / "results" / "adc_loco.tsv", sep="\t")
    pset = set(p_genes)
    for rec in adc_loco.itertuples(index=False):
        if rec.gene_symbol in pset or rec.gene_symbol not in set(adc):
            continue
        if int(rec.n) >= 30 and np.isfinite(rec.s0):
            rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
    rho = stage3.forecast_rho(rho, targets)
    _stats, mean_pct, sd_pct, _ = stage3.cptac_percentile_tables(genes)
    rates = pd.read_csv(ROOT / "results" / "cptac_tumor_protein_rate.tsv", sep="\t")
    rate = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
    gmean = stage3.gmean_all(rho, p_genes)
    mu, sd = stage3.gene_mean_sd(mat.reindex(index=genes), samples)
    mu_pct = stage3.percentile_map(mu)
    sd_pct_ext = stage3.percentile_map(sd)
    rows = []
    for gene in targets:
        rows.append(
            stage3.f1s_feature_row(
                gene, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean,
                mu_pct.get(gene, (np.nan, np.nan))[1],
                sd_pct_ext.get(gene, (np.nan, np.nan))[1],
            ) | {"gene_symbol": gene}
        )
    return pd.DataFrame(rows)


def export(model, x_train, frozen, params, n_train, hats):
    import json

    from lightgbm import Booster

    dest = ROOT / "release" / "models"
    dest.mkdir(parents=True, exist_ok=True)
    model_path = dest / "F1s_lightgbm_rppa_procan.txt"
    model.booster_.save_model(str(model_path))
    reloaded = Booster(model_file=str(model_path))
    gaps = []
    for dataset, (frame, _old, _keys) in frozen.items():
        hat = predict(reloaded, x_train, frame)
        gap = float(np.max(np.abs(hat - hats[dataset])))
        gaps.append({"dataset": dataset, "reload_max_abs": gap})
        if gap != 0:
            raise SystemExit(f"reload {dataset} {gap}")
    settings = {
        "recovery": "parsed feature columns cast to float32 and back to float64 before fill_median",
        "training_rows": "in_P and finite actual",
        "n_train": n_train,
        "features": FEATURES,
        "params": params,
        "limit": 1e-12,
        "seed": 20261001,
        "used_for": {
            "F1s_lightgbm_rppa_procan.txt": ["RPPA", "ProCan"],
            "F1s_lightgbm.txt": ["APOLLO"],
        },
        "reload_max_abs": gaps,
    }
    (dest / "F1s_lightgbm_rppa_procan_settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    pd.DataFrame(gaps).to_csv(OUT / "roundtrip_reload.tsv", sep="\t", index=False)


def main():
    check_counter()
    path = extract()
    table = file_check(path)
    OUT.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT / "roundtrip_check.tsv", sep="\t", index=False)
    form = table["overall_form"].iloc[0]
    print(table.to_string(index=False), flush=True)
    print("overall", form, table["action"].iloc[0], flush=True)
    if form != "float32_origin":
        print("no fit", flush=True)
        return
    refit = fit_recovered(path, table)
    refit.to_csv(OUT / "roundtrip_refit.tsv", sep="\t", index=False)
    print(refit.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
