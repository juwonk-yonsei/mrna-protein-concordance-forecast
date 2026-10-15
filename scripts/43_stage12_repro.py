"""Post-hoc LightGBM thread and histogram grid against the frozen forecasts.

num_threads is set only here. The script does not cap OpenMP, and it does not
import scripts.04_models or scripts.02_cptac_audit.
"""

import os

os.environ["B_CPU_CAP"] = "none"

import importlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from lightgbm import Booster, LGBMRegressor
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

stage3 = importlib.import_module("scripts.08_stage3")

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "stage12"
LIMIT = 1e-12
COMMIT = "224779714bffd58b02059f51de48d92d1324a24f"
F1S_FEATURES = list(stage3.F1S_FEATURES)
THREADS = [1, 16, 72]
HISTOGRAMS = ["auto", "force_col_wise", "force_row_wise"]


def frozen_params():
    b2 = yaml.safe_load(open(ROOT / "config" / "analysis_plan_B2.yaml"))
    stage3.B2 = b2
    stage3.B2_LGBM = dict(b2["reliability"]["lightgbm"])
    return stage3.lgbm_params()


def training_matrix():
    dest = Path("/tmp/b12_train")
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / "f1s_predictions.tsv"
    blob = subprocess.check_output(["git", "show", f"{COMMIT}:results/f1s_predictions.tsv"])
    path.write_bytes(blob)
    preds = pd.read_csv(path, sep="\t")
    if preds["in_P"].dtype == object:
        preds["in_P"] = preds["in_P"].astype(str).str.lower().isin(["true", "1"])
    train = preds[preds["in_P"].astype(bool) & np.isfinite(preds["actual"])].copy()
    x = train[F1S_FEATURES].to_numpy(dtype=np.float64)
    y = train["actual"].to_numpy(dtype=np.float64)
    filled, _, _ = stage3.fill_median(x, x[:1])
    return filled, y, int(len(train))


def predict(model, frame):
    x = frame[F1S_FEATURES].to_numpy(dtype=np.float64)
    filled = np.array(x, dtype=np.float64, copy=True)
    med = np.nanmedian(model._x_tr, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    inds = np.where(~np.isfinite(filled))
    if inds[0].size:
        filled[inds] = np.take(med, inds[1])
    return model.predict(filled)


def fit_one(x, y, params):
    model = LGBMRegressor(**params)
    model.fit(x, y)
    model._x_tr = x
    return model


def compare(new, old, keys):
    merged = new.merge(old[keys + ["f1s_lgbm"]], on=keys, suffixes=("_new", "_old"))
    if len(merged) != len(new) or len(merged) != len(old):
        raise RuntimeError(f"alignment {len(merged)} {len(new)} {len(old)}")
    a = merged["f1s_lgbm_new"].to_numpy(dtype=np.float64)
    b = merged["f1s_lgbm_old"].to_numpy(dtype=np.float64)
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise RuntimeError("non-finite prediction")
    ra = rankdata(a, method="average")
    rb = rankdata(b, method="average")
    return {
        "n_aligned": int(len(merged)),
        "max_abs": float(np.max(np.abs(a - b))),
        "spearman": float(spearmanr(a, b).statistic),
        "n_rank_changed": int(np.sum(ra != rb)),
    }


def apollo_features():
    genes, p_genes, adc = stage3.load_lists()
    targets = pd.read_csv(RES / "stage6" / "prep" / "U_P2_genes.tsv", sep="\t")["gene_symbol"].tolist()
    mat = pd.read_parquet(ROOT / "data" / "processed" / "apollo_luad" / "rna_log2tpm.parquet")
    samples = mat.columns.astype(str).tolist()
    rho = stage3.defined_rho(p_genes)
    adc_loco = pd.read_csv(RES / "adc_loco.tsv", sep="\t")
    pset = set(p_genes)
    for rec in adc_loco.itertuples(index=False):
        if rec.gene_symbol in pset or rec.gene_symbol not in set(adc):
            continue
        if int(rec.n) >= 30 and np.isfinite(rec.s0):
            rho[(rec.heldout_cancer, rec.gene_symbol)] = float(rec.s0)
    rho = stage3.forecast_rho(rho, targets)
    _stats, mean_pct, sd_pct, _ = stage3.cptac_percentile_tables(genes)
    rates = pd.read_csv(RES / "cptac_tumor_protein_rate.tsv", sep="\t")
    rate = {(r.cancer, r.gene_symbol): float(r.rate) for r in rates.itertuples(index=False)}
    gmean = stage3.gmean_all(rho, p_genes)
    mu, sd = stage3.gene_mean_sd(mat.reindex(index=genes), samples)
    mu_pct = stage3.percentile_map(mu)
    sd_pct_ext = stage3.percentile_map(sd)
    rows = []
    for gene in targets:
        ext_mu_pct = mu_pct.get(gene, (np.nan, np.nan))[1]
        hold_pct = sd_pct_ext.get(gene, (np.nan, np.nan))[1]
        rows.append(stage3.f1s_feature_row(gene, stage3.CANCERS, rho, mean_pct, sd_pct, rate, gmean, ext_mu_pct, hold_pct) | {"gene_symbol": gene})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    x, y, n_train = training_matrix()
    print("train", n_train, flush=True)
    rppa = pd.read_csv(RES / "stage9" / "cache_rppa_features.tsv", sep="\t")
    pro = pd.read_csv(RES / "stage9" / "cache_procan_features.tsv", sep="\t")
    print("APOLLO features", flush=True)
    apollo = apollo_features()
    frozen = {
        "RPPA": (rppa, pd.read_csv(RES / "frozen_forecasts_rppa_F1s.tsv", sep="\t"), ["xena_antibody", "cancer"]),
        "ProCan": (pro, pd.read_csv(RES / "frozen_forecasts_procan_F1s.tsv", sep="\t"), ["protein_id", "lineage"]),
        "APOLLO": (apollo, pd.read_csv(RES / "stage6" / "frozen" / "apollo_gene_pred.tsv", sep="\t"), ["gene_symbol"]),
    }
    base = frozen_params()
    rows = []
    saved = {}
    for threads in THREADS:
        for hist in HISTOGRAMS:
            pair = {}
            for repeat in (1, 2):
                params = dict(base)
                params.pop("n_jobs", None)
                params["num_threads"] = int(threads)
                if hist == "force_col_wise":
                    params["force_col_wise"] = True
                elif hist == "force_row_wise":
                    params["force_row_wise"] = True
                print("fit", threads, hist, repeat, flush=True)
                model = fit_one(x, y, params)
                hats = {}
                for dataset, (frame, old, keys) in frozen.items():
                    hat = frame.copy()
                    hat["f1s_lgbm"] = predict(model, frame)
                    hats[dataset] = hat["f1s_lgbm"].to_numpy(dtype=np.float64)
                    rec = compare(hat, old, keys)
                    rec.update({
                        "num_threads": int(threads),
                        "histogram": hist,
                        "repeat": int(repeat),
                        "dataset": dataset,
                        "n_train": n_train,
                        "reproduced": int(rec["max_abs"] <= LIMIT),
                    })
                    rows.append(rec)
                    print(dataset, rec["max_abs"], rec["reproduced"], flush=True)
                pair[repeat] = hats
                saved[(threads, hist, repeat)] = model
            for rec in rows:
                if rec["num_threads"] == threads and rec["histogram"] == hist:
                    a = pair[1][rec["dataset"]]
                    b = pair[2][rec["dataset"]]
                    rec["repeat_pair_max_abs"] = float(np.max(np.abs(a - b)))
    table = pd.DataFrame(rows)
    cols = ["num_threads", "histogram", "repeat", "dataset", "n_train", "n_aligned", "max_abs", "spearman", "n_rank_changed", "reproduced", "repeat_pair_max_abs"]
    table[cols].to_csv(OUT / "repro_grid.tsv", sep="\t", index=False)
    both = table[table["dataset"].isin(["RPPA", "ProCan"])]
    grouped = both.groupby(["num_threads", "histogram"], sort=False)
    winners = []
    for key, sub in grouped:
        if int(sub["reproduced"].min()) == 1 and float(sub["repeat_pair_max_abs"].max()) <= LIMIT:
            winners.append({"num_threads": int(key[0]), "histogram": key[1]})
    pd.DataFrame(winners, columns=["num_threads", "histogram"]).to_csv(OUT / "repro_winners.tsv", sep="\t", index=False)
    print("winners", winners, flush=True)
    if len(winners) == 1:
        chosen = winners[0]
        model = saved[(chosen["num_threads"], chosen["histogram"], 1)]
        text = ROOT / "release" / "models" / "F1s_lightgbm_rppa_procan.txt"
        model.booster_.save_model(str(text))
        settings = {
            "num_threads": chosen["num_threads"],
            "histogram": chosen["histogram"],
            "repeat_written": 1,
            "n_train": n_train,
            "training_commit": COMMIT,
            "base_params": {k: (int(v) if isinstance(v, (np.integer,)) else v) for k, v in base.items()},
        }
        import json
        (ROOT / "release" / "models" / "F1s_lightgbm_rppa_procan_settings.json").write_text(json.dumps(settings, indent=2) + "\n")
        booster = Booster(model_file=str(text))
        reloads = []
        for dataset in ("RPPA", "ProCan"):
            frame, old, keys = frozen[dataset]
            filled = frame[F1S_FEATURES].to_numpy(dtype=np.float64)
            med = np.nanmedian(model._x_tr, axis=0)
            med = np.where(np.isfinite(med), med, 0.0)
            copy = np.array(filled, dtype=np.float64, copy=True)
            inds = np.where(~np.isfinite(copy))
            if inds[0].size:
                copy[inds] = np.take(med, inds[1])
            hat = pd.DataFrame({**{k: frame[k] for k in keys}, "f1s_lgbm": booster.predict(copy)})
            rec = compare(hat, old, keys)
            rec["dataset"] = dataset
            reloads.append(rec)
            print("reload", dataset, rec["max_abs"], flush=True)
        pd.DataFrame(reloads).to_csv(OUT / "repro_reload.tsv", sep="\t", index=False)
    print("stage12 repro done", flush=True)


if __name__ == "__main__":
    main()
