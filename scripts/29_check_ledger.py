"""Check manuscript numbers against manuscript/checks/ledger_B.tsv.

The manuscript path is an argument so a later draft can be checked with the same ledger.
"""
import os
import re
import sys
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "manuscript" / "checks" / "ledger_B.tsv"
CLAIMS = ROOT / "manuscript" / "review_pack" / "claims_ledger.md"
NUM = re.compile(r"(?<![\w.])(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
PATH = re.compile(r"(results/[\w./-]+\.tsv|config/[\w./-]+\.(?:tsv|yaml))")
YEAR = re.compile(r"^(19|20)\d{2}$")


def sig_digits(text, p_value):
    value = float(text.replace("−", "-"))
    digits = 2 if p_value else 3
    return f"{value:.{digits}g}"


def is_p(sentence, start):
    window = sentence[max(0, start - 12):start + 4].lower()
    return "p =" in window or "p=" in window or " p " in window


def build_ledger():
    rows = []
    lines = CLAIMS.read_text().splitlines()[1:]
    key = 0
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        section, sentence, source, status = parts[0], parts[1], parts[2], parts[3]
        paths = PATH.findall(source)
        source_file = paths[0] if paths else ""
        for match in NUM.finditer(sentence):
            token = match.group(1)
            if YEAR.match(token):
                continue
            key += 1
            rows.append(
                {
                    "key": f"n{key:04d}",
                    "location": section,
                    "manuscript_notation": token,
                    "source_file": source_file,
                    "column": "",
                    "row_condition": "",
                    "rounding": "p_2sig" if is_p(sentence, match.start()) else ("exact_integer" if "." not in token and "e" not in token.lower() else "3sig"),
                    "status": status,
                    "sentence": sentence[:240],
                }
            )
    frame = pd.DataFrame(rows)
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(LEDGER, sep="\t", index=False)
    print("ledger rows", len(frame), flush=True)
    return frame


def load_numeric(path):
    frame = pd.read_csv(ROOT / path, sep="\t")
    cols = []
    for col in frame.columns:
        values = pd.to_numeric(frame[col], errors="coerce")
        if values.notna().any():
            cols.append((col, values.to_numpy(dtype=np.float64)))
    return cols


def check(manuscript):
    ledger = pd.read_csv(LEDGER, sep="\t", dtype=str).fillna("")
    cache = {}
    out = []
    for rec in ledger.itertuples(index=False):
        target = rec.manuscript_notation if rec.rounding == "exact_integer" else sig_digits(rec.manuscript_notation, rec.rounding == "p_2sig")
        status = "FAIL"
        detail = "no source file"
        if rec.source_file and (ROOT / rec.source_file).exists():
            if rec.source_file not in cache:
                cache[rec.source_file] = load_numeric(rec.source_file)
            detail = "no matching cell"
            for col, values in cache[rec.source_file]:
                finite = values[np.isfinite(values)]
                if finite.size == 0:
                    continue
                if rec.rounding == "exact_integer":
                    shown = np.array([str(int(v)) if float(v).is_integer() else "" for v in finite])
                else:
                    shown = np.array([f"{v:.2g}" if rec.rounding == "p_2sig" else f"{v:.3g}" for v in finite])
                if np.any(shown == target):
                    status = "PASS"
                    detail = col
                    break
        out.append(
            {
                "key": rec.key,
                "location": rec.location,
                "manuscript_notation": rec.manuscript_notation,
                "source_file": rec.source_file,
                "rounding": rec.rounding,
                "rounded": target,
                "result": status,
                "detail": detail,
            }
        )
    check_path = ROOT / "manuscript" / "checks" / "ledger_check_B9.tsv"
    pd.DataFrame(out).to_csv(check_path, sep="\t", index=False)
    text = Path(manuscript).read_text()
    body = text.split("# Introduction", 1)[-1].split("# Methods", 1)[0]
    abstract = text.split("## Motivation", 1)[-1].split("## Availability", 1)[0]
    known = set(ledger["manuscript_notation"])
    missing = []
    seen = set()
    for chunk, where in ((abstract, "abstract"), (body, "body")):
        for match in NUM.finditer(chunk):
            token = match.group(1)
            if YEAR.match(token) or token in known or token in seen:
                continue
            if token in {"1", "2", "3", "4", "5", "6"} and "Figure" in chunk[max(0, match.start() - 12):match.start() + 8]:
                continue
            seen.add(token)
            missing.append({"location": where, "number": token, "context": chunk[max(0, match.start() - 40):match.end() + 40].replace("\n", " ")})
    miss_path = ROOT / "manuscript" / "checks" / "unledgered_numbers_B9.tsv"
    pd.DataFrame(missing).to_csv(miss_path, sep="\t", index=False)
    n_fail = sum(row["result"] == "FAIL" for row in out)
    print("check", len(out), "fail", n_fail, "unledgered", len(missing), flush=True)


if __name__ == "__main__":
    manuscript = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "manuscript" / "B_manuscript_draft.md")
    if not LEDGER.exists():
        build_ledger()
    check(manuscript)
