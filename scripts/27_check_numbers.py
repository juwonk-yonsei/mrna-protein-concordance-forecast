"""Check the stage-8 headline numbers in the manuscript against the result tables."""
import os
import re
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sig3(value):
    return f"{float(value):.3g}"


def main():
    diff = pd.read_csv(ROOT / "results" / "stage8" / "spearman_diff.tsv", sep="\t")
    text = (ROOT / "manuscript" / "B_manuscript_draft.md").read_text()
    rows = ["check\tfile_value\tmanuscript_has_3sig\tstatus"]
    wanted = diff[
        (diff["left"] == "f1s_lgbm")
        & (diff["right"] == "Lstar")
        & (diff["stratum"].isin(["all", "L1"]))
    ]
    for rec in wanted.itertuples(index=False):
        shown = sig3(rec.diff)
        status = "ok" if shown in text or shown.replace("-", "−") in text else "missing"
        rows.append(f"{rec.dataset} {rec.stratum} f1s-L*\t{rec.diff}\t{shown}\t{status}")
    out = ROOT / "manuscript" / "checks"
    out.mkdir(parents=True, exist_ok=True)
    (out / "numbers.tsv").write_text("\n".join(rows) + "\n")
    print("numbers", sum(1 for row in rows if row.endswith("\tok")), "ok", flush=True)


if __name__ == "__main__":
    main()
