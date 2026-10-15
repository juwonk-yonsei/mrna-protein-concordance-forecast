"""Post-hoc sensitivity on the public non-imputed APOLLO matrix.

Uses the 87 samples in config/apollo_sample_map.tsv. The other sample columns
are not in that map and are not used. Labeled 사후.
"""
import os

os.environ["B_CPU_CAP"] = "none"

import hashlib
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from scripts.common import ROOT

RAW = ROOT / "data" / "raw" / "apollo_luad" / "APOLLO1_Level3_GProt_RAW.txt"
OUT = ROOT / "results" / "stage6" / "posthoc"
LABEL = "사후"


def load_confirm():
    path = ROOT / "scripts" / "17_confirm_B6.py"
    spec = importlib.util.spec_from_file_location("confirm_b6", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    data = RAW.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    mod = load_confirm()
    sample_map = pd.read_csv(ROOT / "config" / "apollo_sample_map.tsv", sep="\t")
    with np.load(ROOT / "results" / "stage6" / "confirm" / "used_predictions.npz", allow_pickle=False) as raw:
        used = {key: np.array(raw[key]) for key in raw.files}
    mids = [str(s) for s in used["samples"]]
    col_of = dict(zip(sample_map["middle_id"].astype(str), sample_map["protein_column"].astype(str)))
    columns = [col_of[mid] for mid in mids]
    frame = pd.read_csv(RAW, sep="\t", usecols=["Gene", *columns])
    frame["Gene"] = frame["Gene"].astype(str)
    counts = frame["Gene"].value_counts()
    dropped = sorted(counts[counts > 1].index)
    frame = frame[~frame["Gene"].isin(set(dropped))]
    frame = frame.drop_duplicates("Gene").set_index("Gene")
    frame = frame[columns]
    level = pd.read_csv(ROOT / "results" / "stage6" / "confirm" / "gene_level.tsv", sep="\t").set_index("gene_symbol")
    u_p1 = set(pd.read_csv(ROOT / "results" / "stage6" / "prep" / "U_P1_genes.tsv", sep="\t")["gene_symbol"])
    u_p2 = set(pd.read_csv(ROOT / "results" / "stage6" / "prep" / "U_P2_genes.tsv", sep="\t")["gene_symbol"])
    index = {str(g): i for i, g in enumerate(used["genes"])}
    rows = []
    for symbol, protein in frame.iterrows():
        if symbol not in index:
            continue
        i = index[symbol]
        yy = protein.to_numpy(dtype=np.float64)
        actual = mod.finite_spearman(used["log2"][i], yy, 30)
        rho_s0 = mod.finite_spearman(used["rna_z"][i], yy, 30) if symbol in u_p1 else np.nan
        rho_ar = mod.finite_spearman(used["ar"][i], yy, 30) if symbol in u_p1 else np.nan
        rec = level.loc[symbol] if symbol in level.index else None
        rows.append(
            {
                "gene_symbol": symbol,
                "in_U_P1": int(symbol in u_p1),
                "in_U_P2": int(symbol in u_p2),
                "actual": actual,
                "d_p1": rho_ar - rho_s0 if np.isfinite(rho_ar) and np.isfinite(rho_s0) else np.nan,
                "f1s_lgbm": float(rec["f1s_lgbm"]) if rec is not None else np.nan,
                "b_sd": float(rec["b_sd"]) if rec is not None else np.nan,
            }
        )
    scored = pd.DataFrame(rows)
    p1 = scored.loc[scored["in_U_P1"] == 1, "d_p1"].to_numpy(dtype=np.float64)
    p1_row = mod.location_row("P1", p1)
    p2_rows = scored[(scored["in_U_P2"] == 1) & np.isfinite(scored["actual"]) & np.isfinite(scored["f1s_lgbm"]) & np.isfinite(scored["b_sd"])]
    p2 = mod.paired_spearman(p2_rows["actual"].to_numpy(), p2_rows["f1s_lgbm"].to_numpy(), p2_rows["b_sd"].to_numpy())
    header = pd.read_csv(RAW, sep="\t", nrows=0).columns.tolist()
    extra = sorted(set(header) - set(columns) - {"Accession", "Protein", "Gene"})
    note = pd.DataFrame(
        [
            {"label": LABEL, "item": "file", "value": str(RAW.relative_to(ROOT))},
            {"label": LABEL, "item": "sha256", "value": digest},
            {"label": LABEL, "item": "size_bytes", "value": str(len(data))},
            {"label": LABEL, "item": "n_sample_columns", "value": str(len(header) - 3)},
            {"label": LABEL, "item": "n_samples_used", "value": "87"},
            {"label": LABEL, "item": "unused_columns", "value": ",".join(extra)},
            {"label": LABEL, "item": "duplicate_symbols_dropped", "value": ",".join(dropped)},
            {"label": LABEL, "item": "P1_n", "value": str(p1_row["n"])},
            {"label": LABEL, "item": "P1_mean", "value": str(p1_row.get("mean", ""))},
            {"label": LABEL, "item": "P1_pvalue", "value": str(p1_row["pvalue"])},
            {"label": LABEL, "item": "P2_n", "value": str(p2["n"])},
            {"label": LABEL, "item": "P2_S_F", "value": str(p2["S_F"])},
            {"label": LABEL, "item": "P2_S_B", "value": str(p2["S_B"])},
            {"label": LABEL, "item": "P2_diff", "value": str(p2["diff"])},
            {"label": LABEL, "item": "P2_pvalue", "value": str(p2["pvalue"])},
        ]
    )
    OUT.mkdir(parents=True, exist_ok=True)
    note.to_csv(OUT / "posthoc_05_raw.tsv", sep="\t", index=False)
    print(note.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
