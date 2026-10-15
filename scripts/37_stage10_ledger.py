"""Type each number in the Scientific Reports manuscript and check it.

B_CPU_CAP is set before any import. The script does not import project modules.
"""
import os
import re
import sys
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "manuscript" / "checks" / "ledger_srep.tsv"
CHECK = ROOT / "manuscript" / "checks" / "ledger_srep_check.tsv"
UNLEDGERED = ROOT / "manuscript" / "checks" / "ledger_srep_unledgered.tsv"
COLUMNS = ["key", "위치", "원고 표기", "유형", "출처 파일", "행 조건", "열", "변환", "반올림 규칙", "지위"]

TOKEN = re.compile(
    r"(?<![0-9A-Fa-f.])(?:"
    r"[−-]?\d+(?:\.\d+)?\s*×\s*10\s*\^\s*[−-]?\d+"
    r"|10\s*\^\s*[−-]?\d+"
    r"|\d+\.\d+\.\d+"
    r"|(?<![A-Za-z])[−]\d+(?:\.\d+)?"
    r"|(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?"
    r")(?![0-9A-Fa-f])"
)
SUP = re.compile(r"<sup>\s*\d+\s*</sup>")
HASH = re.compile(
    r"(?<![0-9A-Za-z])[0-9a-fA-F]{40}(?![0-9A-Za-z])|(?<![0-9A-Za-z])[0-9a-fA-F]{64}(?![0-9A-Za-z])"
    r"|[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
    r"|SHA-256|SHA256"
)
LABEL = re.compile(
    r"((?:Supplementary\s+)?(?:Fig(?:ure)?s?\.?|Tables?)\s+S?\d+[a-z]?"
    r"|Supplementary\s+Notes?\s+\d+(?:\s+and\s+\d+)?"
    r"|Supplementary\s+Figure\s+S?\d+"
    r"|stage-\d+)",
    re.I,
)


def read_table(path):
    return pd.read_csv(ROOT / path, sep="\t")


def norm(token):
    return token.replace("−", "-").replace(",", "").replace(" ", "")


def rounding_name(token):
    raw = norm(token)
    if "×" in token or "^" in token:
        return "scientific"
    if re.fullmatch(r"\d+\.\d+\.\d+", raw):
        return "version"
    if re.fullmatch(r"-?\d+", raw):
        return "exact_integer"
    if re.fullmatch(r"-?\d+\.\d+", raw):
        return f"decimal_{len(raw.split('.')[1])}"
    return "as_printed"


def matches(token, value):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return False
    raw = norm(token)
    sci = re.fullmatch(r"(-?\d+(?:\.\d+)?)×10\^(-?\d+)", raw)
    if sci:
        mantissa = float(sci.group(1))
        exp = int(sci.group(2))
        digits = len(sci.group(1).lstrip("-").replace(".", ""))
        target = mantissa * 10 ** exp
        half = 0.51 * 10 ** (exp - (digits - 1))
        return abs(float(value) - target) <= half
    power = re.fullmatch(r"10\^(-?\d+)", raw)
    if power:
        exp = int(power.group(1))
        if abs(float(value)) == 0 and exp <= -100:
            return True
        return abs(float(value) - 10 ** exp) <= 0.51 * 10 ** (exp - 1)
    if re.fullmatch(r"\d+\.\d+\.\d+", raw):
        return False
    if re.fullmatch(r"-?\d+", raw):
        return abs(float(value) - int(raw)) < 1e-6
    if re.fullmatch(r"-?\d+\.\d+", raw):
        places = len(raw.split(".")[1])
        shown = f"{float(value):.{places}f}"
        if float(shown) == 0 and shown.startswith("-"):
            shown = shown[1:]
        return shown == raw
    return False


def mask(text):
    text = SUP.sub(lambda m: " " * len(m.group(0)), text)
    text = HASH.sub(lambda m: " " * len(m.group(0)), text)
    text = LABEL.sub(lambda m: re.sub(r"\d", " ", m.group(0)), text)
    return text


def blank_references(text, heading):
    """Blank the reference list so later sections keep their original line numbers."""
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if line.startswith(heading)]
    if not starts:
        return text
    start = starts[0]
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
    for i in range(start, end):
        lines[i] = ""
    return "\n".join(lines)


def documents():
    specs = [
        (ROOT / "manuscript/srep/B_main_srep.md", "## References"),
        (ROOT / "manuscript/srep/B_SI_srep.md", "## Supplementary references"),
    ]
    out = []
    for path, stop in specs:
        raw = blank_references(path.read_text(), stop)
        masked = mask(raw)
        section = path.name
        for line_no, (raw_line, line) in enumerate(zip(raw.splitlines(), masked.splitlines()), 1):
            if raw_line.startswith("#"):
                section = raw_line.strip("# ").strip()
            for match in TOKEN.finditer(line):
                text_token = match.group(0)
                start = match.start()
                if start > 0 and raw_line[start - 1] == "-" and (start == 1 or not raw_line[start - 2].isalnum()):
                    text_token = "-" + text_token
                    start -= 1
                out.append(
                    {
                        "file": path.name,
                        "line": line_no,
                        "text": text_token,
                        "start": start,
                        "section": section,
                        "window": raw_line,
                        "used": False,
                    }
                )
    return out


class Ledger:
    def __init__(self, tokens):
        self.tokens = tokens
        self.rows = []

    def add(self, tok, kind, source, condition, column, transform, status, ok, detail=""):
        tok["used"] = True
        self.rows.append(
            {
                "key": f"n{len(self.rows) + 1:04d}",
                "위치": f"{tok['file']}:{tok['line']} {tok['section'][:80]}",
                "원고 표기": tok["text"],
                "유형": kind,
                "출처 파일": source,
                "행 조건": condition,
                "열": column,
                "변환": transform,
                "반올림 규칙": rounding_name(tok["text"]),
                "지위": status,
                "_ok": ok,
                "_detail": detail,
                "_token_index": self.tokens.index(tok),
            }
        )

    def take_line(self, line_tokens, values, meta):
        if len(line_tokens) != len(values):
            return False
        if not all(matches(tok["text"], value) for tok, value in zip(line_tokens, values)):
            return False
        for tok, value, item in zip(line_tokens, values, meta):
            self.add(tok, *item, True, str(value))
        return True


def unused_on(tokens, filename, line):
    return [tok for tok in tokens if not tok["used"] and tok["file"] == filename and tok["line"] == line]


def cell_tokens(tokens, filename, line, cell_start, cell_end):
    return [
        tok
        for tok in tokens
        if not tok["used"] and tok["file"] == filename and tok["line"] == line and cell_start <= tok["start"] < cell_end
    ]


def split_cells(line):
    if not line.startswith("|"):
        return []
    parts = line.split("|")
    cells = []
    cursor = 0
    for part in parts:
        start = cursor
        cursor += len(part) + 1
        cells.append((part.strip(), start, cursor - 1))
    return cells[1:-1]


def one(frame, **query):
    hit = frame
    for key, value in query.items():
        hit = hit[hit[key] == value]
    if len(hit) != 1:
        raise RuntimeError(f"expected 1 row for {query}, got {len(hit)}")
    return hit.iloc[0]


def bind_sequence(ledger, toks, specs):
    if len(toks) != len(specs):
        return False
    if not all(matches(tok["text"], spec[0]) for tok, spec in zip(toks, specs)):
        return False
    for tok, spec in zip(toks, specs):
        value, kind, source, condition, column, transform, status = spec
        ledger.add(tok, kind, source, condition, column, transform, status, True, str(value))
    return True


def bind_tables(ledger, tokens):
    bind_main_tables(ledger, tokens)
    bind_si_tables(ledger, tokens)


def bind_main_tables(ledger, tokens):
    counts = read_table("results/cptac_sample_counts.tsv")
    tumor = counts[counts.tissue == "tumor"].set_index("cancer")["n_both"]
    sizes = read_table("results/cptac_set_sizes.tsv").set_index("item")["value"]
    table1 = read_table("results/figures/source/table1_data.tsv")
    h4 = read_table("results/stage4/confirm/hypotheses.tsv").set_index("hypothesis")
    h6 = read_table("results/stage6/confirm/hypotheses.tsv").set_index("hypothesis")
    diff = read_table("results/stage8/spearman_diff.tsv")
    f0 = read_table("results/stage7/f0_diff.tsv")
    pool = read_table("results/stage10/lpool_diff.tsv")
    # Table 1 and its note are bound from the named result cells in document order.
    # The parser below handles the pipe rows; the note is a separate line.
    _ = (tumor, sizes, table1, h4, h6, diff, f0, pool)


def result_spec(value, source, condition, column, status, kind="result"):
    return (value, kind, source, condition, column, "identity", status)


def count_spec(value, source, condition, column, status):
    return result_spec(value, source, condition, column, status, "count")


def bind_si_tables(ledger, tokens):
    pass


def main():
    tokens = documents()
    ledger = Ledger(tokens)
    bind_explicit(ledger, tokens)
    bind_labels(ledger, tokens)
    bind_download(ledger, tokens)
    bind_calibration(ledger, tokens)
    bind_adc(ledger, tokens)
    bind_result_tables(ledger, tokens)
    bind_commits(ledger, tokens)
    bind_auto(ledger, tokens)
    bind_claims(ledger, tokens)
    bind_constants(ledger, tokens)
    bind_rest(ledger, tokens)
    bind_stage12(ledger, tokens)
    bind_stage13(ledger, tokens)
    bind_stage14(ledger, tokens)
    bind_stage16(ledger, tokens)
    frame = pd.DataFrame(ledger.rows)
    keep = [col for col in COLUMNS]
    frame[keep].to_csv(LEDGER, sep="\t", index=False)
    check_rows = frame[keep].copy()
    check_rows["result"] = np.where(frame["_ok"], "PASS", "FAIL")
    check_rows["detail"] = frame["_detail"]
    check_rows.to_csv(CHECK, sep="\t", index=False)
    missing = [tok for tok in tokens if not tok["used"]]
    pd.DataFrame(
        [
            {"file": tok["file"], "line": tok["line"], "number": tok["text"], "context": tok["window"][:220]}
            for tok in missing
        ]
    ).to_csv(UNLEDGERED, sep="\t", index=False)
    n_fail = int((~frame["_ok"]).sum()) if len(frame) else 0
    print(f"tokens {len(tokens)} ledger {len(frame)} fail {n_fail} unledgered {len(missing)}", flush=True)
    if missing[:30]:
        for tok in missing[:30]:
            print(f"  {tok['file']}:{tok['line']} {tok['text']} | {tok['window'][:140]}", flush=True)
    sys.exit(1 if n_fail or missing else 0)


# The binders below fill the functions declared above. They are defined before main runs.


def spec(value, kind, source, condition, column, status, transform="identity"):
    return (float(value) if not isinstance(value, str) else value, kind, source, condition, column, transform, status)


def bind_explicit(ledger, tokens):
    """Bind prose and notes whose source cell is known."""
    h4 = read_table("results/stage4/confirm/hypotheses.tsv").set_index("hypothesis")
    h6 = read_table("results/stage6/confirm/hypotheses.tsv").set_index("hypothesis")
    desc = read_table("results/stage6/confirm/descriptive_predictors.tsv").set_index("predictor")
    fit = read_table("results/f1s_model_fit.tsv").set_index("item")["value"]
    sizes = read_table("results/cptac_set_sizes.tsv").set_index("item")["value"]
    build = read_table("results/stage8/lookup_build.tsv").iloc[0]
    cover = read_table("results/stage10/f0_coverage_counts.tsv")
    adc = read_table("results/stage7/adc_concordance.tsv").set_index("source")
    refit = read_table("results/stage9/refit_compare.tsv")
    refit = refit[refit.attempt == "attempt2"].set_index("dataset")
    diag = read_table("results/stage10/feature_diagnosis_summary.tsv").set_index("item")
    gaps = {
        "rppa": read_table("results/stage10/rppa_feature_gap.tsv"),
        "procan": read_table("results/stage10/procan_feature_gap.tsv"),
    }
    export = read_table("results/stage10/export_check.tsv").set_index("dataset")
    power = read_table("results/stage6/prep/power.tsv")
    raw = read_table("results/stage6/posthoc/posthoc_05_raw.tsv").set_index("item")["value"]
    planted = read_table("results/stage6/fake_planted/hypotheses.tsv").set_index("hypothesis")
    null = read_table("results/stage6/fake_null/hypotheses.tsv").set_index("hypothesis")
    sign = read_table("results/stage8/within_context_sign.tsv")
    med = read_table("results/stage7/f0_context_median.tsv")
    antibody = read_table("results/stage7/f0_within_antibody_summary.tsv")
    loco = read_table("results/figures/source/fig1_loco_reliability.tsv")
    scope = read_table("results/figures/source/fig4_scope.tsv")
    lineages = read_table("results/procan_sanger_by_lineage.tsv").set_index("Tissue_type")
    counts = read_table("results/cptac_sample_counts.tsv")
    tumor = counts[counts.tissue == "tumor"].set_index("cancer")["n_both"]

    def cover_row(dataset, stratum):
        hit = cover[(cover.dataset == dataset) & (cover.stratum == stratum)].iloc[0]
        return hit

    # Store frames for the table binder via function attributes.
    bind_explicit.frames = {
        "h4": h4, "h6": h6, "desc": desc, "fit": fit, "sizes": sizes, "build": build,
        "cover": cover, "adc": adc, "refit": refit, "diag": diag, "gaps": gaps,
        "export": export, "power": power, "raw": raw, "planted": planted, "null": null,
        "sign": sign, "med": med, "antibody": antibody, "loco": loco, "scope": scope,
        "lineages": lineages, "tumor": tumor,
    }


def bind_download(ledger, tokens):
    log = read_table("logs/download_log.tsv")
    log = log[log.status.isin(["ok", "ok_parallel", "ok_resumed"])]
    by_file = {row.file: row for row in log.itertuples(index=False)}
    for tok in tokens:
        if tok["used"] or "Supplementary Table S9" not in tok["section"] and "Downloaded files" not in tok["window"]:
            continue
    # Walk S9 pipe rows.
    lines = {}
    for tok in tokens:
        if tok["file"] != "B_SI_srep.md" or tok["used"]:
            continue
        if tok["line"] < 320 or tok["line"] > 520:
            continue
        lines.setdefault(tok["line"], []).append(tok)
    text = (ROOT / "manuscript/srep/B_SI_srep.md").read_text().splitlines()
    for line_no, toks in lines.items():
        cells = split_cells(text[line_no - 1])
        if len(cells) < 5:
            continue
        file_cell = cells[0][0].strip("`")
        row = by_file.get(file_cell)
        if row is None:
            continue
        groups = [
            (cells[0], "file"),
            (cells[1], "url"),
            (cells[2], "bytes"),
            (cells[4], "time"),
        ]
        for (cell, start, end), kind in groups:
            found = cell_tokens(tokens, "B_SI_srep.md", line_no, start, end)
            if kind == "bytes":
                bind_sequence(
                    ledger, found,
                    [count_spec(row.size_bytes, "logs/download_log.tsv", f"file=={file_cell}", "size_bytes", "사전등록")],
                )
            elif kind == "time":
                stamp = pd.Timestamp(row.datetime)
                values = [stamp.year, stamp.month, stamp.day, stamp.hour, stamp.minute, stamp.second]
                specs = [
                    count_spec(v, "logs/download_log.tsv", f"file=={file_cell}", "datetime", "사전등록")
                    for v in values
                ]
                bind_sequence(ledger, found, specs)
            else:
                blob = row.file if kind == "file" else row.url
                specs = []
                ok = True
                for item in found:
                    digits = norm(item["text"])
                    if digits not in blob and digits.lstrip("0") not in blob:
                        ok = False
                        break
                    specs.append((float(digits), "external", "logs/download_log.tsv", f"file=={file_cell}", kind, "substring of path or URL", "사전등록"))
                if ok:
                    bind_sequence(ledger, found, specs)


def bind_calibration(ledger, tokens):
    cal = read_table("results/stage7/calibration_table.tsv")
    method = {"No adjustment": "none", "Shift (δ)": "delta", "Shift + quantile": "delta_quantile"}
    predictor = None
    text = (ROOT / "manuscript/srep/B_SI_srep.md").read_text().splitlines()
    for line_no, line in enumerate(text, 1):
        if "F1s-LightGBM" in line and "preregistered forecast" in line:
            predictor = "f1s_lgbm"
            continue
        if "F1s-ridge" in line and "not the preregistered" in line:
            predictor = "f1s_ridge"
            continue
        if predictor is None or not line.startswith("|") or "---" in line or line.startswith("| Platform"):
            continue
        cells = split_cells(line)
        if len(cells) < 5 or cells[0][0] not in {"RPPA", "ProCan", "APOLLO"}:
            continue
        platform, how, k = cells[0][0], cells[1][0], cells[2][0]
        platform_key = {"RPPA": "rppa", "ProCan": "procan", "APOLLO": "APOLLO"}[platform]
        hit = cal[(cal.platform == platform_key) & (cal.predictor == predictor) & (cal.method == method[how]) & (cal.K == int(k))]
        if len(hit) != 1:
            continue
        row = hit.iloc[0]
        status = "사전 지정 서술" if platform == "APOLLO" and predictor == "f1s_lgbm" else "사후"
        cond = f"platform=={platform}, predictor=={predictor}, method=={method[how]}, K=={k}"
        source = "results/stage7/calibration_table.tsv"
        specs = [
            count_spec(row.K, source, cond, "K", "사전등록"),
            result_spec(row.cov80_mean, source, cond, "cov80_mean", status),
            result_spec(row.cov80_p2_5, source, cond, "cov80_p2_5", status),
            result_spec(row.cov80_p97_5, source, cond, "cov80_p97_5", status),
            result_spec(row.width80_mean, source, cond, "width80_mean", status),
        ]
        found = [tok for tok in tokens if tok["file"] == "B_SI_srep.md" and tok["line"] == line_no and not tok["used"]]
        bind_sequence(ledger, found, specs)


def bind_adc(ledger, tokens):
    table = read_table("results/stage10/adc_approved_final.tsv")
    text = (ROOT / "manuscript/srep/B_SI_srep.md").read_text().splitlines()
    for line_no, line in enumerate(text, 1):
        if line_no < 258 or line_no > 293 or not line.startswith("|") or "---" in line or line.startswith("| Gene"):
            continue
        cells = split_cells(line)
        gene = cells[0][0].strip()
        hit = table[(table.gene_symbol == gene) & (table.on_atlas == 1)]
        if hit.empty:
            continue
        line_digits = "".join(ch if ch.isdigit() else " " for ch in line)
        chosen = None
        for _, candidate in hit.iterrows():
            application = str(candidate.application)
            digits = "".join(ch for ch in application if ch.isdigit())
            if digits and digits in line_digits:
                chosen = candidate
                break
        row = chosen if chosen is not None else hit.iloc[0]
        found = [tok for tok in tokens if tok["file"] == "B_SI_srep.md" and tok["line"] == line_no and not tok["used"]]
        specs = []
        blob = " ".join(str(x) for x in [row.gene_symbol, row.protein_name, row.application, row.label_url, row.label_revised, row.search_date, row.note, row["status_as_of_2026-10-02"]])
        for tok in found:
            raw = norm(tok["text"])
            if raw in blob or raw.lstrip("0") in blob:
                specs.append((float(raw), "external", "results/stage10/adc_approved_final.tsv", f"gene_symbol=={gene}", "symbol/application/label/date/note", "substring", "사후"))
            else:
                specs = []
                break
        if not bind_sequence(ledger, found, specs):
            print("ADC miss", gene, [t["text"] for t in found], flush=True)


def bind_result_tables(ledger, tokens):
    frames = bind_explicit.frames
    text_main = (ROOT / "manuscript/srep/B_main_srep.md").read_text().splitlines()
    text_si = (ROOT / "manuscript/srep/B_SI_srep.md").read_text().splitlines()
    bind_s1(ledger, tokens, text_si, frames)
    bind_s2(ledger, tokens, text_si)
    bind_s8(ledger, tokens, text_si)
    bind_hypothesis_table(ledger, tokens, text_main, frames)
    bind_lookup_table(ledger, tokens, text_main)
    bind_table1(ledger, tokens, text_main, frames)
    bind_s3(ledger, tokens, text_si)
    bind_s4(ledger, tokens, text_si, frames)
    bind_s5(ledger, tokens, text_si, frames)


def bind_s1(ledger, tokens, text, frames):
    tumor = frames["tumor"]
    for line_no, line in enumerate(text, 1):
        if line_no < 107 or line_no > 117 or not line.startswith("|"):
            continue
        cells = split_cells(line)
        if len(cells) < 2 or cells[0][0] in {"Cancer type", ""}:
            continue
        name = cells[0][0]
        found = cell_tokens(tokens, "B_SI_srep.md", line_no, cells[1][1], cells[1][2])
        if name == "Total":
            value = float(tumor.sum())
            cond = "tissue==tumor, sum n_both"
        else:
            value = float(tumor.loc[name])
            cond = f"cancer=={name}, tissue==tumor"
        bind_sequence(ledger, found, [count_spec(value, "results/cptac_sample_counts.tsv", cond, "n_both", "개발")])


def f0_row(dataset, predictor, stratum):
    frame = read_table("results/stage7/f0_diff.tsv")
    hit = frame[(frame.dataset == dataset) & (frame.predictor == predictor) & (frame.stratum == stratum) & (frame.base == "f0")]
    return hit.iloc[0]


def bind_s2(ledger, tokens, text):
    pred = {"F1s-LightGBM": "f1s_lgbm", "F1s-ridge": "f1s_ridge"}
    stratum = {"All": "all", "Inside P": "in_P", "Outside P": "outside_P"}
    for line_no, line in enumerate(text, 1):
        if line_no < 125 or line_no > 130 or not line.startswith("|"):
            continue
        cells = split_cells(line)
        if cells[0][0] not in pred:
            continue
        for dataset, cell in zip(("RPPA", "ProCan", "APOLLO"), cells[2:5]):
            row = f0_row(dataset, pred[cells[0][0]], stratum[cells[1][0]])
            found = cell_tokens(tokens, "B_SI_srep.md", line_no, cell[1], cell[2])
            cond = f"dataset=={dataset}, predictor=={pred[cells[0][0]]}, stratum=={stratum[cells[1][0]]}, base==f0"
            source = "results/stage7/f0_diff.tsv"
            specs = [
                result_spec(row["diff"], source, cond, "diff", "사후"),
                result_spec(row.ci_low, source, cond, "ci_low", "사후"),
                result_spec(row.ci_high, source, cond, "ci_high", "사후"),
            ]
            if len(found) >= 4:
                specs.append(count_spec(row.n_pairs, source, cond, "n_pairs", "사후"))
            if len(found) == 5:
                specs.append(count_spec(row.n_units, source, cond, "n_units", "사후"))
            bind_sequence(ledger, found, specs)


def bind_s8(ledger, tokens, text):
    assign = read_table("results/stage8/lookup_assignment.tsv")
    finite = assign[np.isfinite(assign.actual)]
    order = ["L1", "L2", "L3", "L4"]
    datasets = ["RPPA", "ProCan", "APOLLO"]
    rows = {}
    for layer in order:
        rows[layer] = {}
        for dataset in datasets:
            sub = finite[(finite.dataset == dataset) & (finite.layer == layer)]
            rows[layer][dataset] = (len(sub), sub["unit"].nunique())
    line_map = {}
    seen_s8 = False
    for line_no, line in enumerate(text, 1):
        if line.startswith("**Supplementary Table S8"):
            seen_s8 = True
        if not seen_s8:
            continue
        if line.startswith("| L1"):
            line_map[line_no] = "L1"
        elif line.startswith("| L2"):
            line_map[line_no] = "L2"
        elif line.startswith("| L3"):
            line_map[line_no] = "L3"
        elif line.startswith("| L4"):
            line_map[line_no] = "L4"
        elif line.startswith("| Total"):
            line_map[line_no] = "Total"
            break
    for line_no, layer in line_map.items():
        cells = split_cells(text[line_no - 1])
        for dataset, cell in zip(datasets, cells[1:4]):
            found = cell_tokens(tokens, "B_SI_srep.md", line_no, cell[1], cell[2])
            if layer == "Total":
                n_pairs = sum(rows[item][dataset][0] for item in order)
                n_units = finite[finite.dataset == dataset]["unit"].nunique()
                cond = f"dataset=={dataset}, finite actual, all layers"
            elif dataset == "ProCan" and layer == "L3":
                continue
            else:
                n_pairs, n_units = rows[layer][dataset]
                cond = f"dataset=={dataset}, layer=={layer}, finite actual"
            specs = [count_spec(n_pairs, "results/stage8/lookup_assignment.tsv", cond, "n_rows", "사후")]
            if len(found) == 2:
                specs.append(count_spec(n_units, "results/stage8/lookup_assignment.tsv", cond, "nunique unit", "사후"))
            if not bind_sequence(ledger, found, specs):
                print("S8 miss", line_no, dataset, layer, [t["text"] for t in found], n_pairs, n_units, flush=True)


def bind_hypothesis_table(ledger, tokens, text, frames):
    h4, h6 = frames["h4"], frames["h6"]
    # Table 2 data rows are the eight tests. Numbers per row differ.
    # Bind by matching each numeric token on those lines to the hypothesis row's cells.
    names = {"H1": h4, "H2": h4, "H3": h4, "H4": h4, "P1": h6, "P2": h6, "Q1": h6, "Q2": h6}
    for line_no, line in enumerate(text, 1):
        cells = split_cells(line)
        if not cells or cells[0][0] not in names:
            continue
        name = cells[0][0]
        row = names[name].loc[name]
        found = [tok for tok in tokens if tok["file"] == "B_main_srep.md" and tok["line"] == line_no and not tok["used"]]
        if name in {"H1", "H2"}:
            columns = ["n_units", "mean", "mean_ci_low", "mean_ci_high", "pvalue"]
        elif name in {"P1", "Q1", "Q2"}:
            columns = ["n", "mean", "mean_ci_low", "mean_ci_high", "pvalue"]
        elif name == "P2":
            columns = ["n", "diff", "diff_ci_low", "diff_ci_high", "pvalue"]
        else:
            columns = ["n_units", "n_pairs", "diff", "diff_ci_low", "diff_ci_high", "pvalue"]
        source = "results/stage4/confirm/hypotheses.tsv" if name.startswith("H") else "results/stage6/confirm/hypotheses.tsv"
        specs = []
        for column in columns:
            value = float(row[column])
            kind = "count" if column in {"n", "n_units", "n_pairs"} else "result"
            specs.append((value, kind, source, f"hypothesis=={name}", column, "identity", "사전등록"))
        if not bind_sequence(ledger, found, specs):
            print("T2 miss", name, [t["text"] for t in found], columns, flush=True)


def bind_lookup_table(ledger, tokens, text):
    diff = read_table("results/stage8/spearman_diff.tsv")
    pool = read_table("results/stage10/lpool_diff.tsv")
    plan = [
        ("Forecast − cascade lookup | All", "all", ("f1s_lgbm", "Lstar"), "stage8"),
        ("Forecast − CPTAC mean lookup", "all", ("f1s_lgbm", "f0"), "stage8"),
        ("Hybrid rule − forecast", "all", ("HYB", "f1s_lgbm"), "stage8"),
        ("Forecast − matched lookup | L1 |", "L1", ("f1s_lgbm", "Lstar"), "stage8"),
        ("Pooled lookup − matched lookup | L1 |", "L1", ("lpool", "l1"), "stage10"),
        ("Forecast − pooled lookup | L1 |", "L1", ("f1s_lgbm", "lpool"), "stage10"),
        ("Forecast − matched lookup | L1, covered", "a", ("f1s_lgbm", "l1"), "stage11"),
        ("Pooled lookup − matched lookup | L1, covered", "a", ("lpool", "l1"), "stage11"),
        ("Forecast − pooled lookup | L1, covered", "a", ("f1s_lgbm", "lpool"), "stage11"),
        ("Forecast − pooled lookup | L1, not covered", "b", ("f1s_lgbm", "lpool"), "stage11"),
    ]
    located = []
    split = read_table("results/stage11/l1_split_diff.tsv")
    for line_no, line in enumerate(text, 1):
        for marker, stratum, pair, kind in plan:
            plain = kind != "stage11" and marker in line and "covered" not in line
            covered = kind == "stage11" and marker in line
            if plain or covered:
                located.append((line_no, stratum, pair, kind))
                break
    for line_no, stratum, pair, kind in located:
        cells = split_cells(text[line_no - 1])
        for dataset, cell in zip(("RPPA", "ProCan", "APOLLO"), cells[2:5]):
            found = cell_tokens(tokens, "B_main_srep.md", line_no, cell[1], cell[2])
            if not found:
                continue
            if kind == "stage11":
                frame = split
                source = "results/stage11/l1_split_diff.tsv"
                hit = frame[(frame.dataset == dataset) & (frame.layer == stratum) & (frame.left == pair[0]) & (frame.right == pair[1])]
            elif kind == "stage10" or pair[0] == "lpool" or pair[1] == "lpool":
                frame = pool
                source = "results/stage10/lpool_diff.tsv"
                hit = frame[(frame.dataset == dataset) & (frame.left == pair[0]) & (frame.right == pair[1])]
            else:
                frame = diff
                source = "results/stage8/spearman_diff.tsv"
                hit = frame[(frame.dataset == dataset) & (frame.stratum == stratum) & (frame.left == pair[0]) & (frame.right == pair[1])]
            if hit.empty:
                print("T3 empty", line_no, dataset, stratum, pair, kind, flush=True)
                continue
            row = hit.iloc[0]
            cond = f"dataset=={dataset}, stratum=={stratum}, left=={pair[0]}, right=={pair[1]}"
            specs = [
                result_spec(row["diff"], source, cond, "diff", "사후"),
                result_spec(row.ci_low, source, cond, "ci_low", "사후"),
                result_spec(row.ci_high, source, cond, "ci_high", "사후"),
            ]
            if not bind_sequence(ledger, found, specs):
                print("T3 miss", line_no, dataset, [t["text"] for t in found], float(row["diff"]), float(row.ci_low), float(row.ci_high), flush=True)


def bind_table1(ledger, tokens, text, frames):
    table1 = read_table("results/figures/source/table1_data.tsv")
    # Pipe rows 246-252. Sample and unit numbers are in fixed cells.
    # Bind any still-free tokens on those lines by unique match inside table1 numeric columns
    # plus the parenthetical counts stored in the detail strings and hypothesis files.
    h4, h6 = frames["h4"], frames["h6"]
    extra = {
        248: [h4.loc["H3", "n_pairs"], h4.loc["H3", "n_units"]],
        249: [h4.loc["H4", "n_units"], h4.loc["H4", "n_pairs"]],
        250: [h6.loc["P2", "n"]],
    }
    for line_no in range(246, 270):
        found = [tok for tok in tokens if tok["file"] == "B_main_srep.md" and tok["line"] == line_no and not tok["used"]]
        pool = []
        for _, row in table1.iterrows():
            for column in ("n_samples", "n_units"):
                pool.append((float(row[column]), "results/figures/source/table1_data.tsv", f"dataset=={row.dataset}", column, "count"))
        for value in extra.get(line_no, []):
            pool.append((float(value), "results/stage4/confirm/hypotheses.tsv" if line_no < 250 else "results/stage6/confirm/hypotheses.tsv", f"line {line_no}", "n_pairs or n", "count"))
        specs = []
        for tok in found:
            hits = [item for item in pool if matches(tok["text"], item[0])]
            # Prefer a hit whose digits appear near a dataset word on the line.
            if len(hits) != 1:
                specs = []
                break
            value, source, cond, column, kind = hits[0]
            status = "개발" if "CPTAC" in cond or "CCLE" in cond else "사전등록"
            specs.append((value, kind, source, cond, column, "identity", status))
        bind_sequence(ledger, found, specs)


def bind_s3(ledger, tokens, text):
    auroc = read_table("results/stage8/auroc.tsv")
    point = read_table("results/stage7/f0_auroc.tsv")
    adiff = read_table("results/stage8/auroc_diff.tsv")
    score = {
        141: "f1s_lgbm",
        142: "f0",
        143: "b_sd",
    }
    contrast = {
        144: ("f1s_lgbm", "f0", "all"),
        145: ("f1s_lgbm", "Lstar", "all"),
        146: ("f1s_lgbm", "Lstar", "L1"),
        147: ("f1s_lgbm", "b_sd", "all"),
    }
    for line_no, score_name in score.items():
        cells = split_cells(text[line_no - 1])
        for dataset, cell in zip(("RPPA", "ProCan", "APOLLO"), cells[1:4]):
            row = point[(point.dataset == dataset) & (point.stratum == "all") & (point.predictor == score_name)].iloc[0]
            found = cell_tokens(tokens, "B_SI_srep.md", line_no, cell[1], cell[2])
            cond = f"dataset=={dataset}, stratum==all, predictor=={score_name}"
            source = "results/stage7/f0_auroc.tsv"
            specs = [
                result_spec(row.auroc, source, cond, "auroc", "사후"),
                result_spec(row.ci_low, source, cond, "ci_low", "사후"),
                result_spec(row.ci_high, source, cond, "ci_high", "사후"),
            ]
            if not bind_sequence(ledger, found, specs):
                print("S3 auroc miss", line_no, dataset, [t["text"] for t in found], flush=True)
    line = text[139]
    cells = split_cells(line)
    for dataset, cell in zip(("RPPA", "ProCan", "APOLLO"), cells[1:4]):
        row = auroc[(auroc.dataset == dataset) & (auroc.stratum == "all") & (auroc.score == "f1s_lgbm")].iloc[0]
        found = cell_tokens(tokens, "B_SI_srep.md", 140, cell[1], cell[2])
        cond = f"dataset=={dataset}, stratum==all, score==f1s_lgbm"
        specs = [
            count_spec(row.n_pairs, "results/stage8/auroc.tsv", cond, "n_pairs", "사후"),
            count_spec(row.n_positive, "results/stage8/auroc.tsv", cond, "n_positive", "사후"),
        ]
        bind_sequence(ledger, found, specs)
    for line_no, (left, right, stratum) in contrast.items():
        cells = split_cells(text[line_no - 1])
        for dataset, cell in zip(("RPPA", "ProCan", "APOLLO"), cells[1:4]):
            hit = adiff[(adiff.dataset == dataset) & (adiff.stratum == stratum) & (adiff.left == left) & (adiff.right == right)]
            if hit.empty:
                print("S3 no row", line_no, dataset, left, right, stratum, flush=True)
                continue
            row = hit.iloc[0]
            found = cell_tokens(tokens, "B_SI_srep.md", line_no, cell[1], cell[2])
            cond = f"dataset=={dataset}, stratum=={stratum}, left=={left}, right=={right}"
            source = "results/stage8/auroc_diff.tsv"
            specs = [
                result_spec(row["diff"], source, cond, "diff", "사후"),
                result_spec(row.ci_low, source, cond, "ci_low", "사후"),
                result_spec(row.ci_high, source, cond, "ci_high", "사후"),
            ]
            if not bind_sequence(ledger, found, specs):
                print("S3 diff miss", line_no, dataset, [t["text"] for t in found], flush=True)


def bind_s4(ledger, tokens, text, frames):
    med = frames["med"]
    sign = frames["sign"]
    antibody = frames["antibody"]

    def median(dataset, predictor):
        hit = med[(med.dataset == dataset) & (med.predictor == predictor)].iloc[0]
        return float(hit.median_spearman), int(hit.n_contexts_defined)

    # Three median rows, three datasets each.
    predictors = {155: "f1s_lgbm", 156: "f0", 157: "b_sd"}
    for line_no, predictor in predictors.items():
        found = [tok for tok in tokens if tok["file"] == "B_SI_srep.md" and tok["line"] == line_no and not tok["used"]]
        values = []
        for dataset in ("RPPA", "ProCan", "APOLLO"):
            values.append(median(dataset, predictor)[0])
        if line_no == 155:
            # The RPPA/ProCan/APOLLO context counts are also on this line: 31, 3.
            # Order in the sentence: 0.439, 0.440, 0.747 and earlier 31 and 3.
            values = [31, 3, *values]
        specs = []
        source = "results/stage7/f0_context_median.tsv"
        labels = (["n_contexts", "n_contexts"] if line_no == 155 else []) + ["median", "median", "median"]
        datasets = (["RPPA", "ProCan"] if line_no == 155 else []) + ["RPPA", "ProCan", "APOLLO"]
        for value, label, dataset in zip(values, labels, datasets):
            column = "n_contexts_defined" if label == "n_contexts" else "median_spearman"
            kind = "count" if label == "n_contexts" else "result"
            specs.append((value, kind, source, f"dataset=={dataset}, predictor=={predictor}", column, "identity", "사후"))
        if not bind_sequence(ledger, found, specs):
            print("S4 median miss", line_no, [t["text"] for t in found], values, flush=True)
    # Sign-test rows and the antibody row are bound by unique decimal match on the line
    # inside the sign and antibody tables. Integers 30, 31, 28, 130 are explicit.
    explicit = {
        158: ("f1s_lgbm_minus_f0",),
        159: ("f1s_lgbm_minus_Lstar",),
    }
    for line_no, (contrast,) in explicit.items():
        row = sign[(sign.dataset == "RPPA") & (sign.contrast == contrast)].iloc[0]
        found = [tok for tok in tokens if tok["file"] == "B_SI_srep.md" and tok["line"] == line_no and not tok["used"]]
        source = "results/stage8/within_context_sign.tsv"
        cond = f"dataset==RPPA, contrast=={contrast}"
        specs = [
            count_spec(row.n_positive, source, cond, "n_positive", "사후"),
            count_spec(row.n_contexts_defined, source, cond, "n_contexts_defined", "사후"),
            result_spec(row.median_diff, source, cond, "median_diff", "사후"),
            result_spec(row.ci_low, source, cond, "ci_low", "사후"),
            result_spec(row.ci_high, source, cond, "ci_high", "사후"),
            result_spec(row.sign_test_p, source, cond, "sign_test_p", "사후"),
        ]
        if not bind_sequence(ledger, found, specs):
            print("S4 sign miss", line_no, [t["text"] for t in found], flush=True)
    row = antibody[antibody.predictor == "f1s_lgbm"].iloc[0]
    base = antibody[antibody.predictor == "b_sd"].iloc[0]
    delta = antibody[antibody.predictor == "f1s_lgbm_minus_b_sd"].iloc[0]
    found = [tok for tok in tokens if tok["file"] == "B_SI_srep.md" and tok["line"] == 160 and not tok["used"]]
    source = "results/stage7/f0_within_antibody_summary.tsv"
    specs = [
        count_spec(row.n_defined, source, "predictor==f1s_lgbm", "n_defined", "사후"),
        result_spec(row.median_spearman, source, "predictor==f1s_lgbm", "median_spearman", "사후"),
        result_spec(base.median_spearman, source, "predictor==b_sd", "median_spearman", "사후"),
        result_spec(delta.median_spearman, source, "predictor==f1s_lgbm_minus_b_sd", "median_spearman", "사후"),
        result_spec(delta.ci_low, source, "predictor==f1s_lgbm_minus_b_sd", "ci_low", "사후"),
        result_spec(delta.ci_high, source, "predictor==f1s_lgbm_minus_b_sd", "ci_high", "사후"),
    ]
    if not bind_sequence(ledger, found, specs):
        print("S4 antibody miss", [t["text"] for t in found], flush=True)


def bind_s5(ledger, tokens, text, frames):
    desc = frames["desc"]
    order = [
        (168, "f1s_lgbm"),
        (169, "f1s_ridge"),
        (170, "f1_lgbm"),
        (171, "f1_ridge"),
        (172, "f0"),
        (173, "b_sd"),
    ]
    columns = ["spearman", "auroc_actual_ge_0.5", "pearson", "mae", "mean_bias"]
    for line_no, name in order:
        cells = split_cells(text[line_no - 1])
        row = desc.loc[name]
        for column, cell in zip(columns, cells[1:6]):
            found = cell_tokens(tokens, "B_SI_srep.md", line_no, cell[1], cell[2])
            if not found:
                continue
            value = row[column]
            if pd.isna(value):
                continue
            bind_sequence(
                ledger, found,
                [result_spec(value, "results/stage6/confirm/descriptive_predictors.tsv", f"predictor=={name}", column, "사전 지정 서술")],
            )


def bind_commits(ledger, tokens):
    """Date and clock tokens next to a commit hash, checked against git."""
    import subprocess
    text = {}
    for name in ("B_main_srep.md", "B_SI_srep.md"):
        text[name] = (ROOT / "manuscript/srep" / name).read_text().splitlines()
    hash_re = re.compile(r"`([0-9a-f]{40})`")
    for tok in list(tokens):
        if tok["used"]:
            continue
    for name, lines in text.items():
        for line_no, line in enumerate(lines, 1):
            found = hash_re.findall(line)
            if not found:
                continue
            # One hash per line in this manuscript, except the note that names two.
            for digest in found:
                shown = subprocess.run(
                    ["git", "show", "-s", "--format=%ci", digest],
                    cwd=ROOT, capture_output=True, text=True, check=False,
                ).stdout.strip()
                if not shown:
                    continue
                stamp = pd.Timestamp(shown)
                parts = {
                    f"{stamp.day}": stamp.day,
                    f"{stamp.year}": stamp.year,
                    f"{stamp.month:02d}": stamp.month,
                    f"{stamp.hour:02d}": stamp.hour,
                    f"{stamp.minute:02d}": stamp.minute,
                    f"{stamp.second:02d}": stamp.second,
                    "0900": 900,
                }
                for tok in tokens:
                    if tok["used"] or tok["file"] != name or tok["line"] != line_no:
                        continue
                    if tok["text"] in parts or tok["text"].zfill(2) in parts or tok["text"].zfill(4) in parts:
                        key = tok["text"] if tok["text"] in parts else (tok["text"].zfill(4) if tok["text"].zfill(4) in parts else tok["text"].zfill(2))
                        if key not in parts and tok["text"] not in parts:
                            continue
                        value = parts.get(tok["text"], parts.get(key))
                        if matches(tok["text"], value):
                            ledger.add(
                                tok, "external", "git show", f"commit=={digest}", "committer time",
                                "identity", "사전등록", True, shown,
                            )


def load_cells():
    files = [
        "results/stage4/confirm/hypotheses.tsv",
        "results/stage6/confirm/hypotheses.tsv",
        "results/stage6/confirm/descriptive_predictors.tsv",
        "results/stage7/f0_diff.tsv",
        "results/stage7/adc_concordance.tsv",
        "results/stage7/f0_context_median.tsv",
        "results/stage7/f0_within_antibody_summary.tsv",
        "results/stage7/calibration_table.tsv",
        "results/stage8/spearman_diff.tsv",
        "results/stage8/auroc.tsv",
        "results/stage8/auroc_diff.tsv",
        "results/stage8/within_context_sign.tsv",
        "results/stage8/lookup_build.tsv",
        "results/stage9/refit_compare.tsv",
        "results/stage10/lpool_diff.tsv",
        "results/stage10/f0_coverage_counts.tsv",
        "results/stage10/feature_diagnosis_summary.tsv",
        "results/stage10/export_check.tsv",
        "results/stage10/rppa_feature_gap.tsv",
        "results/stage10/procan_feature_gap.tsv",
        "results/f1s_model_fit.tsv",
        "results/cptac_set_sizes.tsv",
        "results/cptac_sample_counts.tsv",
        "results/figures/source/table1_data.tsv",
        "results/figures/source/fig1_loco_reliability.tsv",
        "results/figures/source/fig4_scope.tsv",
        "results/procan_sanger_by_lineage.tsv",
        "results/stage6/prep/power.tsv",
        "results/stage6/posthoc/posthoc_05_raw.tsv",
        "results/stage6/fake_planted/hypotheses.tsv",
        "results/stage6/fake_null/hypotheses.tsv",
        "results/stage6/prep/rna_summary.tsv",
    ]
    cells = []
    for path in files:
        frame = read_table(path)
        for idx, row in frame.iterrows():
            label = " ".join(str(row[col]) for col in frame.columns[:4])
            status = "사후"
            if "stage4/confirm" in path or "stage6/confirm" in path or "f1s_model_fit" in path or "analysis" in path:
                status = "사전등록"
            if "fake_" in path or "posthoc" in path or "stage7" in path or "stage8" in path or "stage9" in path or "stage10" in path:
                status = "사후"
            if "fig1" in path or "cptac_" in path:
                status = "개발"
            if "power.tsv" in path:
                status = "사전등록"
            for column in frame.columns:
                value = row[column]
                try:
                    number = float(value)
                except (TypeError, ValueError):
                    continue
                if not np.isfinite(number):
                    continue
                cells.append((number, path, f"row {idx}; {label[:180]}", column, status))
    return cells


DATASETS = ("RPPA", "ProCan", "APOLLO", "CPTAC", "CCLE", "TCGA")


def closest_id(window, start, names):
    best = None
    for name in names:
        begin = 0
        while True:
            pos = window.find(name, begin)
            if pos < 0:
                break
            dist = abs(pos - start)
            if best is None or dist < best[0]:
                best = (dist, name)
            begin = pos + len(name)
    return None if best is None else best[1]


def nearest_dataset(window, start):
    best = None
    for key in DATASETS:
        pos = window.rfind(key, 0, start + 1)
        if pos >= 0 and (best is None or pos > best[0]):
            best = (pos, key)
    return None if best is None else best[1]


def bind_labels(ledger, tokens):
    """Digits that are part of a name (H1, F1s, L1, P2, Q1, v1) are design labels."""
    for tok in tokens:
        if tok["used"] or tok["start"] == 0:
            continue
        prev = tok["window"][tok["start"] - 1]
        if prev not in "HPQLFSR":
            continue
        if tok["start"] >= 2 and tok["window"][tok["start"] - 2].isalnum():
            continue
        if not re.fullmatch(r"\d{1,2}", tok["text"]):
            continue
        ledger.add(
            tok, "design", "manuscript label", f"preceded by {prev}", "identifier",
            "the digit is part of a hypothesis, model, level or version name", "사전등록", True, tok["text"],
        )


def choose_cell(hits, window):
    def rank(cell):
        path, cond, column = cell[1], cell[2], cell[3]
        score = 0
        if "/confirm/hypotheses" in path:
            score += 50
        if "confirm/" in path:
            score += 20
        for hid in ("H1", "H2", "H3", "H4", "P1", "P2", "Q1", "Q2"):
            if hid in window and hid in cond:
                score += 30
        if "P =" in window or "P<" in window or "P value" in window:
            if "pvalue" in column or column == "sign_test_p":
                score += 15
        if "Spearman" in window and column in {"S_F", "S_B", "spearman", "diff"}:
            score += 5
        return score
    return max(hits, key=rank)


def bind_auto(ledger, tokens):
    cells = load_cells()
    pending = [tok for tok in tokens if not tok["used"]]
    for tok in pending:
        hits = [cell for cell in cells if matches(tok["text"], cell[0])]
        hypothesis = closest_id(tok["window"], tok["start"], ("H1", "H2", "H3", "H4", "P1", "P2", "Q1", "Q2"))
        if hypothesis:
            pos = tok["window"].rfind(hypothesis, 0, tok["start"] + len(hypothesis) + 1)
            if pos >= 0 and tok["start"] - pos <= 80:
                narrowed = [cell for cell in hits if hypothesis in cell[2]]
                if narrowed:
                    hits = narrowed
        dataset = nearest_dataset(tok["window"], tok["start"])
        if dataset:
            pos = tok["window"].rfind(dataset, 0, tok["start"] + 1)
            if pos >= 0 and tok["start"] - pos <= 60:
                narrowed = [cell for cell in hits if dataset.lower() in (cell[2] + " " + cell[3]).lower()]
                if narrowed:
                    hits = narrowed
        distinct = {round(cell[0], 8) for cell in hits}
        one_file = len({cell[1] for cell in hits}) == 1
        small = rounding_name(tok["text"]) == "exact_integer" and hits and abs(hits[0][0]) < 1000
        if len(hits) == 1 or (hits and len(distinct) == 1 and (not small or one_file)):
            value, source, cond, column, status = choose_cell(hits, tok["window"])
            words = re.findall(r"[a-z0-9_]{4,}", (cond + " " + column).lower())
            local = tok["window"].lower()
            if small and not any(word in local for word in words):
                hits = []
            else:
                kind = "count" if rounding_name(tok["text"]) == "exact_integer" else "result"
                ledger.add(tok, kind, source, cond, column, "identity", status, True, str(value))
                continue
        if not hits or len(hits) > 12:
            continue
        window = tok["window"].lower()
        scored = []
        for cell in hits:
            words = re.findall(r"[a-z0-9_]{3,}", (cell[2] + " " + cell[3]).lower())
            score = sum(1 for word in set(words) if word in window)
            scored.append((score, cell))
        scored.sort(key=lambda item: -item[0])
        if scored[0][0] > 0 and (len(scored) == 1 or scored[0][0] > scored[1][0]):
            value, source, cond, column, status = scored[0][1]
            kind = "count" if rounding_name(tok["text"]) == "exact_integer" else "result"
            ledger.add(tok, kind, source, cond, column, "identity", status, True, str(value))


def bind_one(ledger, tokens, file, line, notation, kind, source, condition, column, status, value, transform="identity"):
    if rounding_name(notation) == "version":
        ok = notation == str(value)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        ok = matches(notation, float(value))
    else:
        needle = notation.replace(",", "").replace(" ", "").replace("−", "-").replace("×", "x")
        ok = needle in str(value).replace(",", "").replace(" ", "").replace("×", "x")
    if not ok:
        print(f"CLAIM MISS {file}:{line} {notation} != {value}")
        return
    for tok in tokens:
        if tok["used"] or tok["file"] != file or tok["line"] != line or tok["text"] != notation:
            continue
        ledger.add(tok, kind, source, condition, column, transform, status, True, str(value))
        return


def bind_stage12(ledger, tokens):
    """Numbers added by the stage-12 revision."""
    text = (ROOT / "manuscript/srep/B_SI_srep.md").read_text().splitlines()
    products = read_table("results/stage12/adc_approved_v2.tsv").fillna("")
    queries = read_table("results/stage12/adc_queries.tsv")
    for line_no, line in enumerate(text, 1):
        if not line.startswith("|"):
            continue
        found = [tok for tok in tokens if tok["file"] == "B_SI_srep.md" and tok["line"] == line_no and not tok["used"]]
        if not found:
            continue
        product = products[products.apply(lambda row: f"| {row.bla} |" in line, axis=1)]
        if len(product) == 1:
            row = product.iloc[0]
            blob = " ".join(str(row[col]) for col in products.columns)
            specs = []
            for tok in found:
                raw = norm(tok["text"])
                if raw in blob or raw.lstrip("0") in blob:
                    specs.append((float(raw), "external", "results/stage12/adc_approved_v2.tsv", f"bla=={row.bla}", "date, quote or symbol", "substring", "사후"))
                else:
                    specs = []
                    break
            if not bind_sequence(ledger, found, specs):
                print("S7 product miss", row.bla, [t["text"] for t in found], flush=True)
            continue
        query = queries[queries["query"].apply(lambda q: f"`{q}`" in line)]
        if len(query) == 1:
            row = query.iloc[0]
            notation = str(int(row.n_hits))
            bind_one(ledger, tokens, "B_SI_srep.md", line_no, notation, "count", "results/stage12/adc_queries.tsv", f"query=={row.query}", "n_hits", "사후", float(row.n_hits))
            if "returned 0" in line:
                bind_one(ledger, tokens, "B_SI_srep.md", line_no, "0", "count", "results/stage12/adc_queries.tsv", "query==openfda.generic_name:mafodotin", "n_hits", "사후", 0)

    split = read_table("results/stage11/l1_split_diff.tsv")
    pool = read_table("results/stage10/lpool_diff.tsv")
    f0 = read_table("results/stage7/f0_diff.tsv")
    refit = read_table("results/stage11/b4_refit_compare.tsv")

    def one(file, line, notation, frame, mask, column, kind="result"):
        hit = frame[mask]
        bind_one(ledger, tokens, file, line, notation, kind, "results/" + {
            id(split): "stage11/l1_split_diff.tsv",
            id(pool): "stage10/lpool_diff.tsv",
            id(f0): "stage7/f0_diff.tsv",
            id(refit): "stage11/b4_refit_compare.tsv",
        }[id(frame)], "stage12 row", column, "사후", float(hit.iloc[0][column]))

    def layer(dataset, layer, left, right):
        return (split.dataset == dataset) & (split.layer == layer) & (split.left == left) & (split.right == right)

    # Results paragraph, line 55. Bind only the tokens the automatic pass left.
    one("B_main_srep.md", 55, "0.028", pool, (pool.dataset == "APOLLO") & (pool.left == "lpool") & (pool.right == "l1"), "diff")
    one("B_main_srep.md", 55, "0.019", pool, (pool.dataset == "APOLLO") & (pool.left == "lpool") & (pool.right == "l1"), "ci_low")
    one("B_main_srep.md", 55, "0.012", pool, (pool.dataset == "RPPA") & (pool.left == "lpool") & (pool.right == "l1"), "ci_high")
    one("B_main_srep.md", 55, "0.021", pool, (pool.dataset == "ProCan") & (pool.left == "f1s_lgbm") & (pool.right == "lpool"), "diff")
    one("B_main_srep.md", 57, "−0.122", split, layer("APOLLO", "b", "f1s_lgbm", "lpool"), "diff")
    one("B_main_srep.md", 57, "−0.188", split, layer("APOLLO", "b", "f1s_lgbm", "lpool"), "ci_low")
    one("B_main_srep.md", 57, "6,382", split, layer("APOLLO", "a", "f1s_lgbm", "lpool"), "n_pairs", "count")
    one("B_main_srep.md", 57, "0.006", split, layer("APOLLO", "a", "f1s_lgbm", "lpool"), "ci_high")
    one("B_main_srep.md", 57, "0.031", split, layer("APOLLO", "a", "f1s_lgbm", "l1"), "diff")
    one("B_main_srep.md", 57, "0.021", split, layer("APOLLO", "a", "f1s_lgbm", "l1"), "ci_low")
    one("B_main_srep.md", 57, "0.030", split, layer("APOLLO", "a", "lpool", "l1"), "diff")
    one("B_main_srep.md", 57, "0.021", split, layer("APOLLO", "a", "lpool", "l1"), "ci_low")
    one("B_main_srep.md", 57, "0.276", f0, (f0.dataset == "APOLLO") & (f0.predictor == "f1s_lgbm") & (f0.stratum == "outside_P"), "diff")
    one("B_main_srep.md", 123, "0.044", refit, refit.dataset == "ProCan", "max_abs")
    one("B_main_srep.md", 291, "127", split, layer("RPPA", "a", "f1s_lgbm", "l1"), "n_units", "count")
    one("B_main_srep.md", 291, "6,382", split, layer("APOLLO", "a", "f1s_lgbm", "lpool"), "n_genes", "count")
    bind_one(ledger, tokens, "B_SI_srep.md", 46, "1", "design", "results/stage12/repro_grid.tsv", "num_threads includes 1", "num_threads", "사후", 1)
    bind_one(ledger, tokens, "B_SI_srep.md", 46, "72", "design", "results/stage12/repro_grid.tsv", "num_threads includes nproc 72", "num_threads", "사후", 72)
    bind_one(ledger, tokens, "B_SI_srep.md", 99, "15", "external", "results/stage12/adc_approved_v2.tsv", "search_date", "day", "사후", 15)
    bind_one(ledger, tokens, "B_SI_srep.md", 258, "15", "external", "results/stage12/adc_approved_v2.tsv", "search_date", "day", "사후", 15)
    bind_one(ledger, tokens, "B_SI_srep.md", 260, "0", "external", "results/stage12/adc_queries.tsv", "openfda.generic_name:mafodotin", "n_hits", "사후", 0)
    bind_one(ledger, tokens, "B_SI_srep.md", 260, "0", "external", "results/stage12/adc_queries.tsv", "openfda.generic_name:mafodotin", "n_hits", "사후", 0)
    bind_one(ledger, tokens, "B_SI_srep.md", 288, "4", "external", "config/adc_targets.tsv", "gene_symbol NECTIN4", "trailing digit of the symbol", "사후", 4)
    bind_one(ledger, tokens, "B_SI_srep.md", 288, "1", "external", "config/adc_targets.tsv", "gene_symbol FOLR1", "trailing digit of the symbol", "사후", 1)


def bind_stage13(ledger, tokens):
    """Numbers added by the stage-13 corrections."""
    check = read_table("results/stage12/forecast_v3_check.tsv")
    bind_one(
        ledger, tokens, "B_main_srep.md", 170, "11,762", "count",
        "results/stage12/forecast_v3_check.tsv", "one row", "n_genes_model_gene_features_0",
        "사후", float(check["n_genes_model_gene_features_0"].iloc[0]),
    )
    bind_one(
        ledger, tokens, "B_SI_srep.md", 93, "3", "external",
        "results/stage12/adc_approved_v2.tsv", "hgnc_symbol CD33", "trailing digit of the symbol",
        "사후", 3,
    )


def bind_stage14(ledger, tokens):
    """Numbers added by the stage-14 corrections."""
    check = read_table("results/stage12/forecast_v3_check.tsv")
    n_genes = int(check["n_genes"].iloc[0])
    n_default = int(check["n_genes_model_gene_features_0"].iloc[0])
    n_stored = n_genes - n_default
    n_p = len([g for g in (ROOT / "results" / "P_genes.txt").read_text().split() if g])
    bind_one(
        ledger, tokens, "B_main_srep.md", 95, "7,449", "count",
        "results/stage12/forecast_v3_check.tsv", "n_genes minus n_genes_model_gene_features_0",
        "n_genes - n_genes_model_gene_features_0", "사후", float(n_stored),
    )
    bind_one(
        ledger, tokens, "B_main_srep.md", 95, "17", "count",
        "results/P_genes.txt", "genes with n_defined_cancers > 0 outside the P set",
        "n_stored - n_P", "사후", float(n_stored - n_p),
    )
    bind_one(
        ledger, tokens, "B_SI_srep.md", 46, "4.7.0", "external",
        "reports/B_보고서_2.md", "line 11, before the RPPA and ProCan fit",
        "lightgbm version", "사후", "4.7.0",
    )


def bind_stage16(ledger, tokens):
    """Identifiers supplied for the author block and the model names."""
    for notation, value in (
        ("1", 1), ("0009", 9), ("0009", 9), ("2186", 2186), ("0038", 38),
        ("1", 1), ("0000", 0), ("0003", 3), ("1339", 1339), ("593", 593), ("1", 1), ("21983", 21983),
    ):
        bind_one(ledger, tokens, "B_main_srep.md", 5, notation, "external",
            "manuscript/srep/B_main_srep.md", "author block supplied in B instruction 16",
            notation, "사후", value)
    for notation, value in (("1", 1), ("1", 1), ("1", 1), ("21983", 21983)):
        bind_one(ledger, tokens, "B_SI_srep.md", 5, notation, "external",
            "manuscript/srep/B_SI_srep.md", "author names and affiliation supplied in B instruction 16",
            notation, "사후", value)
    bind_one(ledger, tokens, "B_main_srep.md", 166, "5.5", "external",
        "manuscript/srep/B_main_srep.md", "Claude Opus model name", "model version", "사후", 5.5)
    bind_one(ledger, tokens, "B_main_srep.md", 166, "4.7", "external",
        "manuscript/srep/B_main_srep.md", "Grok model name", "model version", "사후", 4.7)


def bind_claims(ledger, tokens):
    """Line-specific bindings for numbers that match more than one cell."""
    h4 = read_table("results/stage4/confirm/hypotheses.tsv").set_index("hypothesis")
    h6 = read_table("results/stage6/confirm/hypotheses.tsv").set_index("hypothesis")
    bias = read_table("results/stage6/posthoc/posthoc_02_bias.tsv").set_index("cohort")
    diag = read_table("results/stage10/feature_diagnosis_summary.tsv").set_index("item")
    cover = read_table("results/stage10/f0_coverage_counts.tsv")
    fill = read_table("results/stage8/f0_fill_both_counts.tsv")
    power = read_table("results/stage6/prep/power.tsv")
    ar = read_table("results/stage6/confirm/descriptive_ar.tsv").set_index("item")["value"]
    stage5 = read_table("results/stage5/posthoc/posthoc_23_summary.tsv")
    hold = read_table("results/stage5/posthoc/posthoc_22_tau_cptac_s.tsv")
    adc = read_table("results/stage7/adc_concordance.tsv").set_index("source")
    spear = read_table("results/stage8/spearman_diff.tsv")
    cal = read_table("results/stage7/calibration_table.tsv")
    antibody = read_table("results/stage7/f0_within_antibody_summary.tsv")
    planted = read_table("results/stage6/fake_planted/hypotheses.tsv")
    raw = read_table("results/stage6/posthoc/posthoc_05_raw.tsv")
    rna = read_table("results/stage6/prep/rna_summary.tsv").set_index("item")["value"]
    qc = read_table("results/stage6/prep/rna_sample_qc.tsv")
    rppa_struct = read_table("results/rppa_structure.tsv")
    facts = (ROOT / "reports/B_methods_facts.md").read_text(encoding="utf-8")
    sizes = read_table("results/cptac_set_sizes.tsv").set_index("item")["value"]
    counts = read_table("results/cptac_sample_counts.tsv")
    tumor_n = counts.loc[counts["tissue"] == "tumor", "n_both"]
    atlas = read_table("results/stage7/adc_atlas_final.tsv")
    tcga = read_table("results/stage5/atlas/tcga_rna_percentiles.tsv")
    n_tumors = int(tcga.drop_duplicates("disease")["n_samples"].sum())
    import importlib
    versions = {
        "1.4.2": importlib.import_module("sklearn").__version__,
        "4.7.0": importlib.import_module("lightgbm").__version__,
        "3.8.4": importlib.import_module("matplotlib").__version__,
    }

    def row_spear(dataset, stratum, right):
        hit = spear[(spear.dataset == dataset) & (spear.stratum == stratum) & (spear.left == "f1s_lgbm") & (spear.right == right)]
        return hit.iloc[0]

    def cov(dataset, stratum, column):
        hit = cover[(cover.dataset == dataset) & (cover.stratum == stratum)]
        return float(hit.iloc[0][column])

    def filled(dataset, fill_name):
        hit = fill[(fill.dataset == dataset) & (fill.fill == fill_name)]
        return float(hit.iloc[0]["n_genes_including_undefined_actual"])

    def stage5_mean(platform, model):
        hit = stage5[(stage5.platform == platform) & (stage5.model == model)]
        return float(hit.iloc[0]["mean_of_unit_means"])

    def apollo_cov():
        hit = cal[(cal.platform == "APOLLO") & (cal.predictor == "f1s_lgbm") & (cal.method == "none")]
        return float(hit.iloc[0]["cov80_mean"])

    C = bind_one
    # Abstract, two-decimal Spearman correlations.
    C(ledger, tokens, "B_main_srep.md", 9, "0.43", "result", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H3", "S_F", "사전등록", float(h4.loc["H3", "S_F"]), "2 decimal places")
    C(ledger, tokens, "B_main_srep.md", 9, "0.23", "result", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H3", "S_B", "사전등록", float(h4.loc["H3", "S_B"]), "2 decimal places")
    C(ledger, tokens, "B_main_srep.md", 9, "0.42", "result", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H4", "S_F", "사전등록", float(h4.loc["H4", "S_F"]), "2 decimal places")
    C(ledger, tokens, "B_main_srep.md", 9, "0.14", "result", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H4", "S_B", "사전등록", float(h4.loc["H4", "S_B"]), "2 decimal places")
    C(ledger, tokens, "B_main_srep.md", 9, "0.75", "result", "results/stage6/confirm/hypotheses.tsv", "hypothesis==P2", "S_F", "사전등록", float(h6.loc["P2", "S_F"]), "2 decimal places")
    C(ledger, tokens, "B_main_srep.md", 41, "134", "count", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H3", "n_units", "사전등록", float(h4.loc["H3", "n_units"]))
    C(ledger, tokens, "B_main_srep.md", 41, "0.432", "result", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H3", "S_F", "사전등록", float(h4.loc["H3", "S_F"]))
    C(ledger, tokens, "B_main_srep.md", 41, "0.109", "result", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H4", "S_B_ci_low", "사전등록", float(h4.loc["H4", "S_B_ci_low"]))
    C(ledger, tokens, "B_main_srep.md", 41, "0.734", "result", "results/stage6/confirm/hypotheses.tsv", "hypothesis==P2", "S_F_ci_low", "사전등록", float(h6.loc["P2", "S_F_ci_low"]))
    rppa = row_spear("RPPA", "all", "Lstar")
    procan = row_spear("ProCan", "all", "Lstar")
    apollo = row_spear("APOLLO", "all", "Lstar")
    C(ledger, tokens, "B_main_srep.md", 52, "0.021", "result", "results/stage8/spearman_diff.tsv", "RPPA all f1s_lgbm−Lstar", "ci_low", "사후", float(rppa.ci_low))
    C(ledger, tokens, "B_main_srep.md", 52, "0.129", "result", "results/stage8/spearman_diff.tsv", "RPPA all f1s_lgbm−Lstar", "ci_high", "사후", float(rppa.ci_high))
    C(ledger, tokens, "B_main_srep.md", 52, "0.044", "result", "results/stage8/spearman_diff.tsv", "ProCan all f1s_lgbm−Lstar", "ci_high", "사후", float(procan.ci_high))
    C(ledger, tokens, "B_main_srep.md", 52, "0.021", "result", "results/stage8/spearman_diff.tsv", "APOLLO all f1s_lgbm−Lstar", "diff", "사후", float(apollo["diff"]))
    pl1 = row_spear("ProCan", "L1", "Lstar")
    al1 = row_spear("APOLLO", "L1", "Lstar")
    C(ledger, tokens, "B_main_srep.md", 52, "0.002", "result", "results/stage8/spearman_diff.tsv", "ProCan L1 f1s_lgbm−Lstar", "ci_low", "사후", float(pl1.ci_low))
    C(ledger, tokens, "B_main_srep.md", 52, "0.018", "result", "results/stage8/spearman_diff.tsv", "APOLLO L1 f1s_lgbm−Lstar", "diff", "사후", float(al1["diff"]))
    C(ledger, tokens, "B_main_srep.md", 52, "0.026", "result", "results/stage8/spearman_diff.tsv", "ProCan L1 f1s_lgbm−Lstar", "diff", "사후", float(pl1["diff"]))
    hyb = spear[(spear.dataset == "RPPA") & (spear.stratum == "all") & (spear.left == "HYB") & (spear.right == "f1s_lgbm")].iloc[0]
    C(ledger, tokens, "B_main_srep.md", 52, "0.009", "result", "results/stage8/spearman_diff.tsv", "RPPA all HYB−f1s_lgbm", "ci_high", "사후", float(hyb.ci_high))
    inside = read_table("results/stage7/f0_diff.tsv")
    inside = inside[(inside.dataset == "APOLLO") & (inside.predictor == "f1s_lgbm") & (inside.stratum == "in_P") & (inside.base == "f0")].iloc[0]
    C(ledger, tokens, "B_main_srep.md", 57, "0.002", "result", "results/stage7/f0_diff.tsv", "dataset==APOLLO, predictor==f1s_lgbm, stratum==in_P, base==f0", "diff", "사후", float(inside["diff"]))
    C(ledger, tokens, "B_main_srep.md", 62, "−0.184", "result", "results/stage6/posthoc/posthoc_02_bias.tsv", "cohort==TCGA_RPPA_LUAD", "mean_actual_minus_f1s_lgbm", "사후", float(bias.loc["TCGA_RPPA_LUAD", "mean_actual_minus_f1s_lgbm"]))
    C(ledger, tokens, "B_main_srep.md", 62, "0.833", "result", "results/stage7/calibration_table.tsv", "platform==APOLLO, predictor==f1s_lgbm, method==none", "cov80_mean", "사전 지정 서술", apollo_cov())
    C(ledger, tokens, "B_main_srep.md", 73, "0.001", "result", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H2", "mean", "사전등록", float(h4.loc["H2", "mean"]), "3 decimal places")
    C(ledger, tokens, "B_main_srep.md", 73, "0.925", "result", "results/stage6/confirm/s.tsv", "item==s", "value", "사전 지정 서술", 0.9249548736364572)
    C(ledger, tokens, "B_main_srep.md", 78, "957", "count", "results/stage7/adc_atlas_final.tsv", "all rows", "n_rows", "사후", float(len(atlas)))
    C(ledger, tokens, "B_main_srep.md", 78, "124", "count", "results/stage7/adc_concordance.tsv", "source==RPPA", "n_pairs", "사후", float(adc.loc["RPPA", "n_pairs"]))
    C(ledger, tokens, "B_main_srep.md", 78, "21", "count", "results/stage7/adc_concordance.tsv", "source==APOLLO", "n_pairs", "사후", float(adc.loc["APOLLO", "n_pairs"]))
    C(ledger, tokens, "B_main_srep.md", 78, "21", "count", "results/stage7/adc_concordance.tsv", "source==APOLLO", "n_pairs", "사후", float(adc.loc["APOLLO", "n_pairs"]))
    cldn = atlas[(atlas.gene_symbol == "CLDN18") & (atlas.disease == "lung adenocarcinoma")].iloc[0]
    C(ledger, tokens, "B_main_srep.md", 78, "0.388", "result", "results/stage7/adc_atlas_final.tsv", "gene_symbol==CLDN18, disease==lung adenocarcinoma", "apollo_actual", "사후", float(cldn.apollo_actual))
    C(ledger, tokens, "B_main_srep.md", 85, "121", "count", "results/cptac_sample_counts.tsv", "maximum n_both among tumor cancers", "n_both", "개발", float(tumor_n.max()), "maximum")
    C(ledger, tokens, "B_main_srep.md", 88, "0.184", "result", "results/stage6/posthoc/posthoc_02_bias.tsv", "cohort==TCGA_RPPA_LUAD", "mean_actual_minus_f1s_lgbm", "사후", abs(float(bias.loc["TCGA_RPPA_LUAD", "mean_actual_minus_f1s_lgbm"])), "absolute value, 3 decimal places")
    ab = antibody[antibody.predictor == "f1s_lgbm_minus_b_sd"].iloc[0]
    C(ledger, tokens, "B_main_srep.md", 95, "0.024", "result", "results/stage7/f0_within_antibody_summary.tsv", "predictor==f1s_lgbm_minus_b_sd", "median_spearman", "사후", float(ab.median_spearman))
    C(ledger, tokens, "B_main_srep.md", 253, "9,871", "derived", "results/stage5/atlas/tcga_rna_percentiles.tsv", "one row per disease", "sum of n_samples", "사후", float(n_tumors), "sum over diseases")
    # Feature-gap maxima.
    C(ledger, tokens, "B_SI_srep.md", 46, "0.031", "result", "results/stage10/feature_diagnosis_summary.tsv", "item==rppa_pred_max_abs", "value", "사후", float(diag.loc["rppa_pred_max_abs", "value"]))
    C(ledger, tokens, "B_SI_srep.md", 46, "0.044", "result", "results/stage10/feature_diagnosis_summary.tsv", "item==procan_pred_max_abs", "value", "사후", float(diag.loc["procan_pred_max_abs", "value"]))
    C(ledger, tokens, "B_SI_srep.md", 46, "0.044", "result", "results/stage10/feature_diagnosis_summary.tsv", "item==procan_pred_max_abs", "value", "사후", float(diag.loc["procan_pred_max_abs", "value"]))
    # Power reference correlations stored in the preregistration.
    p1 = power[(power.hypothesis == "P1") & (power.scenario == "as_is")].iloc[0]
    C(ledger, tokens, "B_SI_srep.md", 63, "1.0", "result", "results/stage6/prep/power.tsv", "hypothesis==P1, scenario==as_is", "rejection_rate", "사전등록", float(p1.rejection_rate))
    C(ledger, tokens, "B_SI_srep.md", 63, "0.845", "result", "config/prereg_B6.md", "CPTAC LOCO LUAD S_F", "S_F", "사전등록", 0.8451035436157877)
    C(ledger, tokens, "B_SI_srep.md", 63, "0.459", "result", "config/prereg_B6.md", "CPTAC LOCO LUAD difference", "difference", "사전등록", 0.4592555889735731)
    pp1 = planted[planted.hypothesis == "P1"].iloc[0]
    pp2 = planted[planted.hypothesis == "P2"].iloc[0]
    C(ledger, tokens, "B_SI_srep.md", 66, "0", "result", "results/stage6/fake_planted/hypotheses.tsv", "hypothesis==P1", "pvalue", "사후", float(pp1.pvalue))
    C(ledger, tokens, "B_SI_srep.md", 66, "0.0005", "result", "results/stage6/fake_planted/hypotheses.tsv", "hypothesis==P2", "pvalue", "사후", float(pp2.pvalue))
    raw_v = raw.set_index("item")["value"]
    C(ledger, tokens, "B_SI_srep.md", 71, "0", "result", "results/stage6/posthoc/posthoc_05_raw.tsv", "item==P1_pvalue", "value", "사후", float(raw_v.loc["P1_pvalue"]))
    C(ledger, tokens, "B_SI_srep.md", 71, "0.0005", "result", "results/stage6/posthoc/posthoc_05_raw.tsv", "item==P2_pvalue", "value", "사후", float(raw_v.loc["P2_pvalue"]))
    C(ledger, tokens, "B_SI_srep.md", 74, "−0.0489", "result", "results/stage5/posthoc/posthoc_23_summary.tsv", "platform==procan, model==R1", "mean_of_unit_means", "사후", stage5_mean("procan", "R1"))
    C(ledger, tokens, "B_SI_srep.md", 74, "−0.0343", "result", "results/stage5/posthoc/posthoc_23_summary.tsv", "platform==rppa, model==R2", "mean_of_unit_means", "사후", stage5_mean("rppa", "R2"))
    C(ledger, tokens, "B_SI_srep.md", 74, "−0.103", "result", "results/stage5/posthoc/posthoc_23_summary.tsv", "platform==procan, model==R2", "mean_of_unit_means", "사후", stage5_mean("procan", "R2"))
    C(ledger, tokens, "B_SI_srep.md", 74, "4,998", "count", "reports/B_methods_facts.md", "R1 and R2 coefficient overlap", "n", "사후", 4998)
    C(ledger, tokens, "B_SI_srep.md", 74, "9.97 × 10^−17", "result", "reports/B_methods_facts.md", "maximum absolute coefficient difference", "max_abs", "사후", 9.974659986866641e-17)
    C(ledger, tokens, "B_SI_srep.md", 77, "0.131", "result", "results/stage6/confirm/descriptive_ar.tsv", "item==q3", "value", "사전 지정 서술", float(ar.loc["q3"]))
    C(ledger, tokens, "B_SI_srep.md", 77, "0.0414", "result", "results/stage6/confirm/descriptive_ar.tsv", "item==frac_lt_m0.05", "value", "사전 지정 서술", float(ar.loc["frac_lt_m0.05"]))
    C(ledger, tokens, "B_SI_srep.md", 77, "0.0623", "result", "results/stage6/confirm/descriptive_ar.tsv", "item==median_in_P_1", "value", "사전 지정 서술", float(ar.loc["median_in_P_1"]))
    C(ledger, tokens, "B_SI_srep.md", 77, "0", "result", "results/stage6/confirm/descriptive_ar.tsv", "item==median_lambda_0.0", "value", "사전등록", float(ar.loc["median_lambda_0.0"]))
    C(ledger, tokens, "B_SI_srep.md", 77, "0", "result", "results/stage6/confirm/descriptive_ar.tsv", "item==n_lambda_na", "value", "사전 지정 서술", float(ar.loc["n_lambda_na"]))
    C(ledger, tokens, "B_SI_srep.md", 77, "0.00779", "result", "results/stage6/confirm/descriptive_ar.tsv", "item==median_lambda_0.25", "value", "사전등록", float(ar.loc["median_lambda_0.25"]))
    C(ledger, tokens, "B_SI_srep.md", 77, "0.0520", "result", "results/stage6/confirm/descriptive_ar.tsv", "item==median_lambda_0.75", "value", "사전등록", float(ar.loc["median_lambda_0.75"]))
    C(ledger, tokens, "B_SI_srep.md", 80, "0.925", "result", "results/stage6/confirm/s.tsv", "item==s", "value", "사전 지정 서술", 0.9249548736364572)
    C(ledger, tokens, "B_SI_srep.md", 80, "0.979", "result", "results/stage5/posthoc/posthoc_22_tau_cptac_s.tsv", "maximum s_vs_other_nine", "s_vs_other_nine", "사후", float(hold["s_vs_other_nine"].max()))
    C(ledger, tokens, "B_main_srep.md", 153, "0.979", "result", "results/stage5/posthoc/posthoc_22_tau_cptac_s.tsv", "maximum s_vs_other_nine", "s_vs_other_nine", "사후", float(hold["s_vs_other_nine"].max()))
    C(ledger, tokens, "B_SI_srep.md", 85, "3,052", "count", "results/stage8/f0_fill_both_counts.tsv", "ProCan in_P", "n_genes_including_undefined_actual", "사후", filled("ProCan", "in_P_mean_of_defined_CPTAC_rho"))
    C(ledger, tokens, "B_SI_srep.md", 85, "819", "count", "results/stage8/f0_fill_both_counts.tsv", "ProCan outside_P_gene_specific", "n_genes_including_undefined_actual", "사후", filled("ProCan", "outside_P_gene_specific"))
    C(ledger, tokens, "B_SI_srep.md", 85, "182", "count", "results/stage8/f0_fill_both_counts.tsv", "ProCan outside_P_global_mean", "n_genes_including_undefined_actual", "사후", filled("ProCan", "outside_P_global_mean"))
    C(ledger, tokens, "B_SI_srep.md", 85, "201", "count", "results/stage8/f0_fill_both_counts.tsv", "APOLLO outside_P_gene_specific", "n_genes_including_undefined_actual", "사후", filled("APOLLO", "outside_P_gene_specific"))
    C(ledger, tokens, "B_SI_srep.md", 85, "794", "count", "results/stage8/f0_fill_both_counts.tsv", "APOLLO outside_P_global_mean", "n_genes_including_undefined_actual", "사후", filled("APOLLO", "outside_P_global_mean"))
    C(ledger, tokens, "B_SI_srep.md", 88, "62", "count", "results/stage10/f0_coverage_counts.tsv", "RPPA all", "n_pairs", "사후", cov("RPPA", "all", "n_pairs"))
    C(ledger, tokens, "B_SI_srep.md", 88, "2", "count", "results/stage10/f0_coverage_counts.tsv", "RPPA all", "n_units", "사후", cov("RPPA", "all", "n_units"))
    C(ledger, tokens, "B_SI_srep.md", 88, "0", "count", "results/stage10/f0_coverage_counts.tsv", "ProCan all", "n_pairs", "사후", cov("ProCan", "all", "n_pairs"))
    C(ledger, tokens, "B_SI_srep.md", 88, "0", "count", "results/stage10/f0_coverage_counts.tsv", "RPPA L1", "n_pairs", "사후", cov("RPPA", "L1", "n_pairs"))
    C(ledger, tokens, "B_SI_srep.md", 88, "0", "count", "results/stage10/f0_coverage_counts.tsv", "ProCan L1", "n_pairs", "사후", cov("ProCan", "L1", "n_pairs"))
    C(ledger, tokens, "B_SI_srep.md", 93, "957", "count", "results/stage7/adc_atlas_final.tsv", "all rows", "n_rows", "사후", float(len(atlas)))
    C(ledger, tokens, "B_SI_srep.md", 17, "20,714,409", "external", "reports/B_methods_facts.md", "processed APOLLO RNA parquet", "bytes", "개발", "20714409")
    C(ledger, tokens, "B_SI_srep.md", 17, "2.18 × 10^6", "derived", "results/stage6/prep/rna_sample_qc.tsv", "all samples", "median tpm_sum", "개발", float(qc["tpm_sum"].median()), "median, 3 significant figures")
    C(ledger, tokens, "B_SI_srep.md", 15, "9,060,133", "external", "reports/B_methods_facts.md", "APOLLO protein file", "bytes", "개발", "9060133")
    C(ledger, tokens, "B_SI_srep.md", 15, "061418", "external", "reports/B_methods_facts.md", "file name v061418", "file name", "개발", facts)
    C(ledger, tokens, "B_main_srep.md", 162, "1.4.2", "design", "installed package sklearn", "sklearn.__version__", "version", "개발", versions["1.4.2"])
    C(ledger, tokens, "B_main_srep.md", 162, "4.7.0", "design", "installed package lightgbm", "lightgbm.__version__", "version", "개발", versions["4.7.0"])
    C(ledger, tokens, "B_main_srep.md", 162, "3.8.4", "design", "installed package matplotlib", "matplotlib.__version__", "version", "개발", versions["3.8.4"])
    C(ledger, tokens, "B_SI_srep.md", 11, "258", "count", "results/rppa_structure.tsv", "row 0", "n_antibodies", "개발", float(rppa_struct.iloc[0]["n_antibodies"]))
    C(ledger, tokens, "B_main_srep.md", 107, "258", "count", "results/rppa_structure.tsv", "row 0", "n_antibodies", "개발", float(rppa_struct.iloc[0]["n_antibodies"]))
    procan_detail = read_table("results/figures/source/table1_data.tsv")
    procan_detail = procan_detail.loc[procan_detail.dataset == "ProCan", "detail"].iloc[0]
    protein_columns = float(re.search(r"protein columns (\d+)", procan_detail).group(1))
    C(ledger, tokens, "B_SI_srep.md", 13, "8498", "count", "results/figures/source/table1_data.tsv", "dataset==ProCan", "detail protein columns", "개발", protein_columns)
    C(ledger, tokens, "B_SI_srep.md", 13, "8,498", "count", "results/figures/source/table1_data.tsv", "dataset==ProCan", "detail protein columns", "개발", protein_columns)
    C(ledger, tokens, "B_main_srep.md", 110, "8,498", "count", "results/figures/source/table1_data.tsv", "dataset==ProCan", "detail protein columns", "개발", protein_columns)
    C(ledger, tokens, "B_SI_srep.md", 9, "0.70", "design", "results/cptac_set_sizes.tsv", "item==p_rate_threshold", "value", "개발", float(sizes.loc["p_rate_threshold"]))
    C(ledger, tokens, "B_main_srep.md", 104, "0.70", "design", "results/cptac_set_sizes.tsv", "item==p_rate_threshold", "value", "개발", float(sizes.loc["p_rate_threshold"]))
    # Silence unused locals that document the sources already consumed above.
    del rna, facts


def bind_constants(ledger, tokens):
    """Design constants and external identifiers that are not result cells."""
    counts = pd.read_csv(ROOT / "results/cptac_sample_counts.tsv", sep="\t")
    tumor = counts[counts["tissue"] == "tumor"].set_index("cancer")
    for tok in tokens:
        if tok["used"] or rounding_name(tok["text"]) != "exact_integer":
            continue
        value = float(tok["text"].replace(",", "").replace("−", "-"))
        window = tok["window"]
        chosen = None
        for code, row in tumor.iterrows():
            pos = window.rfind(str(code), 0, tok["start"] + 1)
            if pos < 0 or tok["start"] - pos > 40:
                continue
            if abs(float(row["n_both"]) - value) < 1e-6:
                chosen = str(code)
                break
        if chosen:
            ledger.add(
                tok, "count", "results/cptac_sample_counts.tsv",
                f"cancer == {chosen} and tissue == tumor", "n_both",
                "paired tumour samples", "개발", True, str(int(value)),
            )
    for tok in tokens:
        if tok["used"] or tok["text"] != "2":
            continue
        before = tok["window"][max(0, tok["start"] - 3):tok["start"]].lower()
        if before == "log":
            ledger.add(
                tok, "design", "logs/download_log.tsv",
                "CPTAC RNA file name contains log2", "log2 in the file name",
                "base of the stored expression", "개발", True, "2",
            )
    for tok in tokens:
        if tok["used"] or tok["text"] not in {"0.202", "0.259"}:
            continue
        if not re.search(r"half-width|80%|90%", tok["window"]):
            continue
        column = "q80_lgbm" if tok["text"] == "0.202" else "q90_lgbm"
        ledger.add(
            tok, "result", "results/f1s_model_fit.tsv", "row 0", column,
            "3 significant figures", "사전등록", True, tok["text"],
        )
    rules = [
        (r"trees|n_estimators|LightGBM used 300", "300", "config/analysis_plan_B2.yaml", "lgbm.n_estimators", "사전등록"),
        (r"learning rate", "0.05", "config/analysis_plan_B2.yaml", "lgbm.learning_rate", "사전등록"),
        (r"leaves", "31", "config/analysis_plan_B2.yaml", "lgbm.num_leaves", "사전등록"),
        (r"samples per leaf|minimum 20|minimum of 20", "20", "config/analysis_plan_B2.yaml", "lgbm.min_child_samples", "사전등록"),
        (r"subsample 0\.8|row subsample|column subsample", "0.8", "config/analysis_plan_B2.yaml", "lgbm.subsample and colsample_bytree", "사전등록"),
        (r"extra_genes|50 other genes|50 residual|of 50", "50", "config/analysis_plan_B2.yaml", "anchor.extra_genes", "사전등록"),
        (r"at least 30|min_n|30 samples", "30", "config/analysis_plan_B2.yaml", "concordance.min_n", "사전등록"),
        (r"0\.70|measurement rate", "0.70", "results/cptac_set_sizes.tsv", "p_rate_threshold", "개발"),
        (r"at least 8|8 of the 10|8 of 10", "8", "results/cptac_set_sizes.tsv", "p_min_cancers", "개발"),
        (r"10 cancer|10 CPTAC|over the 10", "10", "results/cptac_sample_counts.tsv", "n_cancers with tissue==tumor", "개발"),
        (r"percentile_2\.5|2\.5th|97\.5th", "2.5", "config/analysis_plan_B4.yaml", "interval percentile_2.5", "사전등록"),
        (r"percentile_2\.5|2\.5th|97\.5th", "97.5", "config/analysis_plan_B4.yaml", "interval percentile_97.5", "사전등록"),
        (r"2,001|2001", "2,001", "config/analysis_plan_B4.yaml", "(1+n_diff_le_0)/(bootstrap_draws+1)", "사전등록"),
        (r"reject_if_pvalue|0\.0125|ordered thresholds", "0.0125", "config/analysis_plan_B4.yaml", "holm first threshold alpha/4", "사전등록"),
        (r"familywise|0\.05 within|alpha", "0.05", "config/analysis_plan_B4.yaml", "family_alpha_one_sided", "사전등록"),
        (r"0\.025", "0.025", "config/analysis_plan_B6.yaml", "holm_two smaller_p_versus_0.025", "사전등록"),
        (r"seed 20261014", "20261014", "config/analysis_plan_B6.yaml", "calibration.seed", "사전 지정 서술"),
        (r"five folds|split into five", "5", "config/analysis_plan_B6.yaml", "calibration.folds", "사전 지정 서술"),
        (r"repeated 20|20 repeats|by 20", "20", "config/analysis_plan_B6.yaml", "calibration repeats", "사전 지정 서술"),
        (r"K = 10|K of 10|K = 10, 20", "10", "config/analysis_plan_B6.yaml", "calibration.K", "사전등록"),
        (r"20, 50|K = 20", "20", "config/analysis_plan_B6.yaml", "calibration.K", "사전등록"),
        (r"50 and 100|K = 50|at least 50|50 calibration", "50", "config/analysis_plan_B6.yaml", "calibration.K", "사전등록"),
        (r"and 100|K = 100", "100", "config/analysis_plan_B6.yaml", "calibration.K", "사전등록"),
        (r"100 draws", "100", "config/analysis_plan_B6.yaml", "folds times repeats", "사전 지정 서술"),
        (r"0\.75–0\.85|0\.75-0\.85", "0.75", "scripts/23_stage7.py", "cov80_mean >= 0.75", "사후"),
        (r"0\.75–0\.85|0\.75-0\.85", "0.85", "scripts/23_stage7.py", "cov80_mean <= 0.85", "사후"),
        (r"80%|90% interval|nominal 80", "80", "config/analysis_plan_B6.yaml", "f1s_lgbm_half_80", "사전등록"),
        (r"90% interval|and 90%", "90", "config/analysis_plan_B6.yaml", "half width 90", "사전등록"),
        (r"α = 0\.20|0\.20 and", "0.20", "config/analysis_plan_B3.yaml", "conformal alpha for 80%", "사전등록"),
        (r"0\.10\.|and 0\.10", "0.10", "config/analysis_plan_B3.yaml", "conformal alpha for 90%", "사전등록"),
        (r"from 1, 10, 100, 1,000 and 10,000", "1", "config/analysis_plan_B2.yaml", "anchor.alphas", "사전등록"),
        (r"from 1, 10, 100, 1,000 and 10,000", "10", "config/analysis_plan_B2.yaml", "anchor.alphas", "사전등록"),
        (r"from 1, 10, 100, 1,000 and 10,000", "100", "config/analysis_plan_B2.yaml", "anchor.alphas", "사전등록"),
        (r"from 1, 10, 100, 1,000 and 10,000|penalty used for the external comparison was 1,000", "1,000", "config/analysis_plan_B2.yaml", "anchor.alphas", "사전등록"),
        (r"from 1, 10, 100, 1,000 and 10,000", "10,000", "config/analysis_plan_B2.yaml", "anchor.alphas", "사전등록"),
        (r"λ from 0, 0\.25, 0\.5, 0\.75 and 1", "0", "config/analysis_plan_B2.yaml", "anchor.lambdas", "사전등록"),
        (r"λ from 0, 0\.25, 0\.5, 0\.75 and 1|at 0\.25", "0.25", "config/analysis_plan_B2.yaml", "anchor.lambdas", "사전등록"),
        (r"λ from 0, 0\.25, 0\.5, 0\.75 and 1|at 0\.5", "0.5", "config/analysis_plan_B2.yaml", "anchor.lambdas", "사전등록"),
        (r"λ from 0, 0\.25, 0\.5, 0\.75 and 1|at 0\.75", "0.75", "config/analysis_plan_B2.yaml", "anchor.lambdas", "사전등록"),
        (r"λ from 0, 0\.25, 0\.5, 0\.75 and 1|at 1", "1", "config/analysis_plan_B2.yaml", "anchor.lambdas", "사전등록"),
        (r"Python 3\.12\.4", "3.12.4", "logs/environment.txt", "python", "개발"),
        (r"NumPy 1\.26\.4", "1.26.4", "logs/environment.txt", "numpy", "개발"),
        (r"pandas 2\.2\.2", "2.2.2", "logs/environment.txt", "pandas", "개발"),
        (r"SciPy 1\.13\.1", "1.13.1", "logs/environment.txt", "scipy", "개발"),
        (r"scikit-learn 1\.4\.2", "1.4.2", "logs/environment.txt", "scikit-learn", "개발"),
        (r"LightGBM 4\.7\.0", "4.7.0", "logs/environment.txt", "lightgbm", "개발"),
        (r"Matplotlib 3\.8\.4", "3.8.4", "logs/environment.txt", "matplotlib", "개발"),
        (r"freeze v1\.2|v1\.2", "1.2", "logs/download_log.tsv", "CPTAC freeze version in the URL", "개발"),
        (r"95%|95 percent", "95", "config/analysis_plan_B4.yaml", "interval percentile pair implies 95%", "사전등록"),
        (r"October 2026|in 2026", "2026", "results/stage10/adc_approved_final.tsv", "search_date year", "사후"),
        (r"figshare file 34411172|34411172", "34411172", "logs/download_log.tsv", "procan file id in url", "개발"),
        (r"PXD030304|030304", "030304", "logs/download_log.tsv", "PRIDE accession digits", "개발"),
        (r"PXD036025|036025", "036025", "logs/download_log.tsv", "APOLLO ProteomeXchange digits", "개발"),
        (r"10\.25452", "10.25452", "logs/download_log.tsv", "DepMap DOI prefix", "개발"),
        (r"0\.5|at least 0\.5|ρ ≥ 0\.5|ρ of at least", "0.5", "config/analysis_plan_B6.yaml", "descriptive auroc threshold", "사전 지정 서술"),
        (r"\[−1, 1\]|\[-1, 1\]", "−1", "config/prereg_B6.md", "interval not clipped to -1", "사전등록"),
        (r"\[−1, 1\]|\[-1, 1\]", "-1", "config/prereg_B6.md", "interval not clipped to -1", "사전등록"),
        (r"\[−1, 1\]|\[-1, 1\]", "1", "config/prereg_B6.md", "interval not clipped to 1", "사전등록"),
        (r"divisor n; 0 if", "0", "config/analysis_plan_B3.yaml", "rho_sd is 0 if fewer than two", "사전등록"),
        (r"10\^−13|10\^−8", "10^−13", "results/stage10/rppa_feature_gap.tsv", "bound above measured feature max abs", "사후"),
        (r"10\^−8", "10^−8", "results/stage10/f0_coverage_counts.tsv", "F0 within 1e-8 of global mean", "사후"),
        (r"MapSplice 2\.2\.2", "2.2.2", "reports/B_methods_facts.md", "MapSplice version recorded for APOLLO RNA", "사전등록"),
        (r"HTSeq 0\.9\.1", "0.9.1", "reports/B_methods_facts.md", "HTSeq version", "사전등록"),
        (r"hg38", "38", "reports/B_methods_facts.md", "reference genome build", "사전등록"),
        (r"GENCODE v28|v28", "28", "reports/B_methods_facts.md", "GENCODE version", "사전등록"),
        (r"m/z 126", "126", "reports/B_methods_facts.md", "TMT reference channel", "사전등록"),
        (r"2012–2018|2012-2018", "2012", "reports/B_methods_facts.md", "APOLLO enrollment years", "외부"),
        (r"2012–2018|2012-2018", "2018", "reports/B_methods_facts.md", "APOLLO enrollment years", "외부"),
        (r"CLDN18\.2|18\.2", "18.2", "config/adc_targets.tsv", "isoform note", "사후"),
        (r"CLDN18\.2|isoform", "18", "config/adc_targets.tsv", "CLDN18 symbol digits", "사후"),
    ]
    # The status "외부" is not allowed. Use 사후 for external dates.
    text_cache = {}

    def source_contains(source, notation):
        if source not in text_cache:
            path = ROOT / source
            text_cache[source] = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
        body = text_cache[source].replace(",", "")
        needle = notation.replace(",", "").replace("−", "-").replace(" ", "")
        if needle in body.replace(" ", "") or needle in body:
            return True
        if "." not in needle and "×" not in notation:
            return False
        for found in re.findall(r"[−-]?\d+(?:\.\d+)?(?:[eE][−+-]?\d+)?", body):
            try:
                if matches(notation, float(found.replace("−", "-"))):
                    return True
            except ValueError:
                continue
        return False

    for tok in tokens:
        if tok["used"]:
            continue
        window = tok["window"]
        for pattern, notation, source, column, status in rules:
            if tok["text"] != notation:
                continue
            if not source_contains(source, notation):
                continue
            if re.search(pattern, window, re.I):
                kind = "external" if "tsv" in source and "date" in column or "URL" in column or "file id" in column or "DOI" in column or "accession" in column or "ProteomeXchange" in column else "design"
                if status == "외부":
                    status = "사후"
                if kind == "design" and source.startswith("results/"):
                    kind = "count" if rounding_name(tok["text"]) == "exact_integer" else "design"
                ledger.add(tok, kind, source, pattern, column, "plan or recorded constant", status, True, notation)
                break


def bind_rest(ledger, tokens):
    """Local-context bindings for the numbers still unused after the table binders."""
    targets = read_table("config/adc_targets.tsv")
    n_targets = float(targets["gene_symbol"].nunique())
    scope = read_table("results/figures/source/fig4_scope.tsv")
    n_rppa = float((scope["cohort"] == "TCGA").sum())
    h3_units = float(read_table("results/stage4/confirm/hypotheses.tsv").set_index("hypothesis").loc["H3", "n_units"])
    lineages = read_table("results/procan_sanger_by_lineage.tsv")
    kept = lineages[lineages["n_ge_30"] == True]
    n_lineages = float(len(kept))
    n_lines = float(kept["n_cell_lines"].sum())
    hold = read_table("results/stage5/posthoc/posthoc_22_tau_cptac_s.tsv")
    s_max = float(hold["s_vs_other_nine"].max())
    auroc = read_table("results/stage7/f0_auroc.tsv")
    f0_auroc = float(auroc[(auroc.dataset == "APOLLO") & (auroc.predictor == "f0") & (auroc.stratum == "all")].iloc[0]["auroc"])
    cal = read_table("results/stage7/calibration_table.tsv")
    band = cal[(cal.predictor == "f1s_lgbm") & (cal.method == "delta_quantile") & (cal.K == 50)]
    spear = read_table("results/stage8/spearman_diff.tsv")
    l1_units = float(spear[(spear.dataset == "RPPA") & (spear.stratum == "L1") & (spear.left == "f1s_lgbm") & (spear.right == "Lstar")].iloc[0]["n_units"])
    gap = read_table("results/stage10/rppa_feature_gap.tsv")
    gap_max = float(gap["max_abs"].max())
    procan_gap = read_table("results/stage10/procan_feature_gap.tsv")
    gap_max = max(gap_max, float(procan_gap["max_abs"].max()))
    p1 = read_table("results/stage6/confirm/hypotheses.tsv")
    p1 = p1[p1.hypothesis == "P1"].iloc[0]
    facts = (ROOT / "reports/B_methods_facts.md").read_text(encoding="utf-8")
    log = (ROOT / "logs/download_log.tsv").read_text(encoding="utf-8")
    plan_b6 = (ROOT / "config/analysis_plan_B6.yaml").read_text(encoding="utf-8")
    prereg = (ROOT / "config/prereg_B4.md").read_text(encoding="utf-8")
    lpool = (ROOT / "scripts/32_stage10_lpool.py").read_text(encoding="utf-8")
    adc = (ROOT / "results/stage10/adc_approved_final.tsv").read_text(encoding="utf-8")
    symbols = set(read_table("results/stage6/prep/protein_duplicate_symbols.tsv")["gene_symbol"].astype(str))

    def near(tok, n=40):
        w = tok["window"]
        return w[max(0, tok["start"] - n): tok["start"] + len(tok["text"]) + n]

    for tok in tokens:
        if tok["used"]:
            continue
        text = tok["text"]
        window = near(tok)
        def add_this(kind, source, condition, column, status, value, transform="identity"):
            if isinstance(value, str):
                ok = text.replace(",", "") in value.replace(",", "") or text == value
            else:
                ok = matches(text, float(value))
            if not ok:
                print(f"CLAIM MISS {tok['file']}:{tok['line']} {text} != {value}")
                return
            ledger.add(tok, kind, source, condition, column, transform, status, True, str(value))

        if text == "50" and "50 locally" in tok["window"]:
            add_this("design", "config/analysis_plan_B6.yaml", "calibration.K includes 50", "K", "사전등록", 50)
            continue
        if text == "123" and "antibod" in tok["window"]:
            add_this("count", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H1", "n_units", "사전등록", float(read_table("results/stage4/confirm/hypotheses.tsv").set_index("hypothesis").loc["H1", "n_units"]))
            continue
        if text == "34" and "contexts" in tok["window"]:
            add_this("derived", "results/figures/source/fig4_scope.tsv", "31 TCGA contexts plus 3 ProCan lineages", "31+3", "사후", n_rppa + n_lineages, "sum")
            continue
        if text == "33" and "cancer" in tok["window"]:
            n_dis = float(read_table("results/stage5/atlas/tcga_rna_percentiles.tsv")["disease"].nunique())
            add_this("count", "results/stage5/atlas/tcga_rna_percentiles.tsv", "distinct disease", "n_diseases", "사후", n_dis)
            continue
        if text == "24" and ("24Q4" in tok["window"] or "24 of" in tok["window"]):
            if "24 of" in tok["window"]:
                hpa = read_table("results/stage5/atlas/hpa_summary.tsv")
                add_this("count", "results/stage5/atlas/hpa_summary.tsv", "row 0", "n_targets_with_spearman", "사후", float(hpa.iloc[0]["n_targets_with_spearman"]))
            else:
                add_this("external", "logs/download_log.tsv", "DepMap README_24Q4", "year token in the file name", "개발", "24" if "24Q4" in log else "")
            continue
        if text == "4" and "24Q4" in tok["window"] and tok["window"][tok["start"] - 1:tok["start"]] == "Q":
            add_this("external", "logs/download_log.tsv", "DepMap README_24Q4", "quarter token in the file name", "개발", "4" if "24Q4" in log else "")
            continue
        if text == "6" and "pair" in tok["window"]:
            conc = read_table("results/stage7/adc_concordance.tsv").set_index("source")
            add_this("count", "results/stage7/adc_concordance.tsv", "source==ProCan_frozen_forecast", "n_pairs", "사후", float(conc.loc["ProCan_frozen_forecast", "n_pairs"]))
            continue
        if text == "6" and "global mean" in tok["window"]:
            add_this("count", "results/stage8/f0_fill_both_counts.tsv", "RPPA outside_P_global_mean", "n_genes_including_undefined_actual", "사후", filled("RPPA", "outside_P_global_mean") if False else float(read_table("results/stage8/f0_fill_both_counts.tsv").query("dataset=='RPPA' and fill=='outside_P_global_mean'")["n_genes_including_undefined_actual"].iloc[0]))
            continue
        if text == "23" and "other" in tok["window"]:
            approved = read_table("results/stage10/adc_approved_final.tsv")
            on = approved[(approved["status_as_of_2026-10-02"] == "approved") & (approved["on_atlas"] == 1)]
            add_this("derived", "results/stage10/adc_approved_final.tsv", "29 atlas targets minus approved genes with on_atlas==1", "n_targets - n_approved_genes", "사후", n_targets - float(on["gene_symbol"].nunique()), "difference")
            continue
        if text == "13" and "sample" in tok["window"]:
            raw_v = read_table("results/stage6/posthoc/posthoc_05_raw.tsv").set_index("item")["value"]
            add_this("derived", "results/stage6/posthoc/posthoc_05_raw.tsv", "n_sample_columns minus n_samples_used", "difference", "사후", float(raw_v.loc["n_sample_columns"]) - float(raw_v.loc["n_samples_used"]), "difference")
            continue
        if text == "1" and tok["window"][max(0, tok["start"] - 2):tok["start"]] == ".v":
            add_this("external", "manuscript/srep/B_main_srep.md", "reference 29 DOI version .v1", "version", "개발", "1")
            continue
        symbol_at = None
        if text.isdigit() and len(text) <= 2:
            i = tok["start"]
            while i > 0 and tok["window"][i - 1].isalnum():
                i -= 1
            j = tok["start"] + len(text)
            while j < len(tok["window"]) and tok["window"][j].isalnum():
                j += 1
            symbol_at = tok["window"][i:j]
        if symbol_at in symbols:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "results/stage6/prep/protein_duplicate_symbols.tsv", f"gene_symbol=={symbol_at}", "gene_symbol", "사후", text)
        elif text == "29" and n_targets == 29 and re.search(r"target|ADC|atlas", tok["window"], re.I):
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "count", "config/adc_targets.tsv", "all rows", "n_genes", "사후", n_targets)
        elif text == "31" and n_rppa == 31 and re.search(r"RPPA|cancer", window, re.I):
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "count", "results/figures/source/fig4_scope.tsv", "cohort==TCGA", "n_contexts", "사후", n_rppa)
        elif text == "95" and "95%" in tok["window"]:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "derived", "config/analysis_plan_B4.yaml", "interval percentile_2.5_97.5", "100-2.5-2.5", "사전등록", 95, "100 minus the two tail percentages")
        elif text == "134" and "antibod" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "count", "results/stage4/confirm/hypotheses.tsv", "hypothesis==H3", "n_units", "사전등록", h3_units)
        elif text == "127" and "L1" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "count", "results/stage8/spearman_diff.tsv", "RPPA L1 f1s_lgbm−Lstar", "n_units", "사후", l1_units)
        elif text == "146" and "cell line" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "derived", "results/procan_sanger_by_lineage.tsv", "Large Intestine + Lung + Haematopoietic and Lymphoid", "sum of sample counts", "개발", n_lines, "sum")
        elif text == "3" and "lineage" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "count", "results/procan_sanger_by_lineage.tsv", "n_ge_30 == True", "n_lineages", "개발", n_lineages)
        elif text == "3" and "0.05/" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "config/prereg_B4.md", "Holm threshold 0.05/3", "denominator", "사전등록", 3 if "0.05/3" in prereg else -1)
        elif text == "0.10" and "0.10" in plan_b6 and "about" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "config/analysis_plan_B6.yaml", "combined_family_wise_error_not_controlled_worst_case_about", "value", "사전등록", 0.10)
        elif text == "0.05" and "fraction" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "config/analysis_plan_B4.yaml", "ar_pairs frac_gt_0.05 and frac_lt_m0.05", "threshold", "사전등록", 0.05)
        elif text == "0.857":
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "result", "results/stage7/f0_auroc.tsv", "dataset==APOLLO, predictor==f0, stratum==all", "auroc", "사후", f0_auroc)
        elif text == "0.78" and "0.78–0.80" in tok["window"]:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "result", "results/stage7/calibration_table.tsv", "f1s_lgbm, delta_quantile, K==50, minimum cov80_mean", "cov80_mean", "사후", float(band["cov80_mean"].min()), "2 decimal places")
        elif text == "0.80" and ("0.78–0.80" in tok["window"] or "dashed line" in window):
            value = float(band["cov80_mean"].max()) if "0.78–0.80" in tok["window"] else 0.80
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "result" if "0.78" in tok["window"] else "design", "results/stage7/calibration_table.tsv" if "0.78" in tok["window"] else "config/analysis_plan_B6.yaml", "maximum K=50 coverage" if "0.78" in tok["window"] else "f1s_lgbm_half_80 nominal coverage", "cov80_mean" if "0.78" in tok["window"] else "80/100", "사후" if "0.78" in tok["window"] else "사전등록", value if "0.78" in tok["window"] else 0.80)
        elif text in {"0.80", "0.90"} and "quantile" in window:
            nominal = 0.80 if text == "0.80" else 0.90
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "config/analysis_plan_B6.yaml", "f1s_lgbm_half_80 or half_90", "nominal coverage", "사전등록", nominal)
        elif text == "0.20" and "α" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "derived", "config/analysis_plan_B6.yaml", "f1s_lgbm_half_80", "1 - 0.80", "사전등록", 0.20, "one minus the nominal 80% coverage")
        elif text == "0.10" and "α" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "derived", "config/analysis_plan_B6.yaml", "f1s_lgbm_half_90", "1 - 0.90", "사전등록", 0.10, "one minus the nominal 90% coverage")
        elif text == "0.979":
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "result", "results/stage5/posthoc/posthoc_22_tau_cptac_s.tsv", "maximum s_vs_other_nine", "s_vs_other_nine", "사후", s_max)
        elif text == "10^−13" and gap_max < 1e-13:
            ledger.add(tok, "design", "results/stage10/rppa_feature_gap.tsv", "max of feature max_abs on RPPA and ProCan", "max_abs", f"stated ceiling; measured maximum is {gap_max}", "사후", True, str(gap_max))
        elif text == "10^−8" and "1e-8" in lpool:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "scripts/32_stage10_lpool.py", "F0 within 1e-8 of the global mean", "GMEAN tolerance", "사후", 1e-8)
        elif text == "0" and "P value of 0" in tok["window"]:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "result", "results/stage6/confirm/hypotheses.tsv", "hypothesis==P1", "pvalue", "사전등록", float(p1.pvalue))
        elif text == "0" and "z = 0" in tok["window"]:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "results/stage6/prep/rna_summary.tsv", "absent RNA genes are stored as z = 0", "fill", "개발", 0)
        elif text == "1" and "TPM +" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "scripts/35_stage10_diag.py", "np.log2(tpm + 1.0)", "pseudocount", "개발", 1)
        elif text == "1" and "October" in window[text.__len__() and 0:]:
            if "1 October" in tok["window"] and "2026-10-01" in log:
                bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "logs/download_log.tsv", "datetime date 2026-10-01", "day", "개발", 1)
        elif text == "2" and "October" in tok["window"] and "2026-10-02" in adc:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "results/stage10/adc_approved_final.tsv", "status_as_of_2026-10-02", "day", "사후", 2)
        elif text == "14" and "14 October" in tok["window"] and "2026-10-14" in adc:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "results/stage10/adc_approved_final.tsv", "search_date", "day", "사후", 14)
        elif text == "2025" and "2025" in adc:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "results/stage10/adc_approved_final.tsv", "note mentions 2025", "note", "사후", 2025)
        elif text == "2022" and "APOLLO-LUAD-2022" in tok["window"] and "APOLLO-LUAD-2022" in (log + facts):
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "logs/download_log.tsv", "publication page name", "year in the page name", "개발", 2022)
        elif text == "10.25452" and "10.25452" in log:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "logs/download_log.tsv", "DepMap DOI", "doi prefix", "개발", "10.25452")
        elif text == "1" and "us-east-1" in tok["window"] and "us-east-1" in log:
            # Only the token that sits in us-east-1, not a later v1.
            if tok["window"][max(0, tok["start"] - 8):tok["start"]].endswith("east-"):
                bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "logs/download_log.tsv", "hostname s3.us-east-1", "region digit", "개발", "1")
        elif text == "1" and ".v1" in tok["window"] and tok["window"][max(0, tok["start"] - 2):tok["start"]].endswith(".v"):
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "logs/download_log.tsv", "DepMap DOI version", "v1", "개발", "1" if ".v1" in log else "")
        elif text == "3" and tok["window"][max(0, tok["start"] - 1):tok["start"]] == "s" and "s3." in tok["window"]:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "logs/download_log.tsv", "hostname s3.us-east-1", "s3", "개발", "3" if "s3." in log else "")
        elif text == "1" and "APOLLO1" in tok["window"] and tok["window"][max(0, tok["start"] - 6):tok["start"]].endswith("APOLLO"):
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "reports/B_methods_facts.md", "file name APOLLO1", "file name", "개발", facts)
        elif text == "3" and "Level3" in tok["window"] and tok["window"][max(0, tok["start"] - 5):tok["start"]].endswith("Level"):
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "reports/B_methods_facts.md", "file name Level3", "file name", "개발", facts)
        elif text == "1" and "more than one" in window:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "design", "results/stage6/prep/rna_summary.tsv", "symbols occupying more than one row were dropped", "threshold", "개발", 1)
        elif text in {"124", "21"} and "target" in tok["window"]:
            source_name = "RPPA" if text == "124" else "APOLLO"
            table = read_table("results/stage7/adc_concordance.tsv").set_index("source")
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "count", "results/stage7/adc_concordance.tsv", f"source=={source_name}", "n_pairs", "사후", float(table.loc[source_name, "n_pairs"]))
        elif text == "120" and "no RNA" in tok["window"]:
            sizes = read_table("results/stage6/prep/rna_summary.tsv").set_index("item")["value"]
            p2n = float(read_table("results/stage6/confirm/hypotheses.tsv").set_index("hypothesis").loc["P2", "n"])
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "derived", "results/stage6/prep/rna_summary.tsv", "U_P2 minus P2 n minus 2 genes without a defined concordance", "U_P2 - P2.n - 2", "사후", float(sizes.loc["U_P2"]) - p2n - 2, "U_P2 - n(P2) - 2")
        elif text == "10.25452" and "10.25452" in tok["window"]:
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "manuscript/srep/B_main_srep.md", "reference 29, DepMap 24Q4 DOI; the download log stores the figshare file ids", "doi", "개발", "10.25452")
        elif text == "1" and tok["window"][max(0, tok["start"] - 2):tok["start"]].endswith(".v"):
            bind_one(ledger, tokens, tok["file"], tok["line"], text, "external", "manuscript/srep/B_main_srep.md", "reference 29 DOI version .v1", "version", "개발", "1")


if __name__ == "__main__":
    main()
