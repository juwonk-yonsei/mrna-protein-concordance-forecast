"""Post-hoc openFDA and Drugs@FDA search for approved antibody-drug conjugates.

The cutoff is 2026-10-02. A 404 from the API is a count of zero. Any other
failure is retried twice and then stops the search.
"""

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage12"
CUTOFF = date(2026, 10, 2)
SEARCH_DATE = date(2026, 10, 15)
STEMS = [
    "vedotin", "mafodotin", "deruxtecan", "govitecan", "emtansine",
    "soravtansine", "tesirine", "ozogamicin", "duocarmazine", "tirumotecan",
]
REQUIRED = ["ENHERTU", "KADCYLA", "TRODELVY", "DATROWAY", "PADCEV", "TIVDAK", "ELAHERE", "EMRELIS", "BLENREP"]
# Longest alias first. Each alias was read in a mechanism-of-action sentence.
ALIASES = [
    ("folate receptor alpha", "FOLR1"),
    ("B-cell maturation antigen", "TNFRSF17"),
    ("tissue factor", "F3"),
    ("Nectin-4", "NECTIN4"),
    ("nectin-4", "NECTIN4"),
    ("Trop-2", "TACSTD2"),
    ("c-Met", "MET"),
    ("CD79b", "CD79B"),
    ("CD123", "IL3RA"),
    ("CD30", "TNFRSF8"),
    ("CD33", "CD33"),
    ("CD22", "CD22"),
    ("CD19", "CD19"),
    ("HER2", "ERBB2"),
    ("BCMA", "TNFRSF17"),
    ("FRα", "FOLR1"),
]


def fetch(url):
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as handle:
                return json.load(handle)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return {"meta": {"results": {"total": 0}}, "results": []}
            last = exc
        except urllib.error.URLError as exc:
            last = exc
        if attempt < 2:
            time.sleep(2)
    raise RuntimeError(f"openFDA failed after 2 retries: {last}")


def search(endpoint, query, note):
    url = f"https://api.fda.gov/drug/{endpoint}.json?limit=100&search=" + urllib.parse.quote(query)
    data = fetch(url)
    total = int(data["meta"]["results"]["total"])
    rows = list(data.get("results", []))
    skip = 100
    while len(rows) < total:
        page = fetch(url + f"&skip={skip}")
        rows.extend(page.get("results", []))
        skip += 100
    return {"endpoint": endpoint, "query": query, "n_hits": total, "note": note, "results": rows}


def ymd(text):
    text = str(text)
    return date(int(text[0:4]), int(text[4:6]), int(text[6:8]))


def first_approval(record):
    dates = []
    for sub in record.get("submissions", []):
        if sub.get("submission_type") == "ORIG" and sub.get("submission_status") == "AP" and sub.get("submission_status_date"):
            dates.append(ymd(sub["submission_status_date"]))
    return min(dates) if dates else None


def sentence(text):
    text = " ".join(str(text).split())
    body = text.replace("12.1 Mechanism of Action ", "")
    pieces = [piece.strip() for piece in body.split(". ") if piece.strip()]
    target = next((piece for piece in pieces if symbol_of(piece)), "")
    adc = next((piece for piece in pieces if ("antibody" in piece.lower() and "conjugate" in piece.lower()) or " is an adc" in piece.lower()), "")
    chosen = [piece.rstrip(".") + "." for piece in (target, adc) if piece]
    unique = list(dict.fromkeys(chosen))
    if unique:
        return " ".join(unique)
    return pieces[0].rstrip(".") + "." if pieces else ""


def symbol_of(quote):
    found = []
    for alias, symbol in ALIASES:
        if alias in quote and symbol not in found:
            found.append(symbol)
    return ";".join(found)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    logged = []
    logged.append(search("label", '"antibody-drug conjugate"', "label text contains the phrase"))
    for stem in STEMS:
        logged.append(search("drugsfda", f"openfda.generic_name:{stem}", "planned generic-name stem query"))
    # generic_name:mafodotin returned no record. Blenrep is indexed on the active-ingredient name.
    for stem in STEMS:
        logged.append(search(
            "drugsfda",
            f"products.active_ingredients.name:{stem}",
            "added after openfda.generic_name:mafodotin returned 0 and missed Blenrep",
        ))
    logged.append(search("label", '"belantamab"', "added because the Blenrep application number has no openFDA label record"))
    logged.append(search("drugsfda", "products.active_ingredients.name:sunirine", "phrase hit whose payload stem is not in the planned list"))

    products = {}
    for item in logged:
        if item["endpoint"] != "drugsfda" or "active_ingredients.name" not in item["query"]:
            continue
        for rec in item["results"]:
            app = rec.get("application_number")
            if not app or app in products:
                continue
            brand = rec["products"][0]["brand_name"]
            generic = rec["products"][0]["active_ingredients"][0]["name"]
            products[app] = {"brand": brand, "generic": generic, "bla": app, "first_approval": first_approval(rec)}

    labels = {}
    for rec in logged[0]["results"]:
        of = rec.get("openfda") or {}
        apps = of.get("application_number") or []
        moa = (rec.get("mechanism_of_action") or [""])[0]
        quote = sentence(moa)
        is_adc = "is an antibody" in quote.lower() or "is a " in quote.lower() and "conjugate" in quote.lower()
        for app in apps:
            labels.setdefault(app, []).append((rec.get("effective_time"), quote, is_adc, of.get("brand_name")))
    for rec in next(item for item in logged if item["query"] == '"belantamab"')["results"]:
        moa = (rec.get("mechanism_of_action") or [""])[0]
        labels.setdefault("BLA761440", []).append((rec.get("effective_time"), sentence(moa), True, ["BLENREP"]))

    rows = []
    for app, prod in sorted(products.items()):
        choices = labels.get(app, [])
        if not choices:
            donor = None
            for other, info in products.items():
                if other != app and info["brand"] == prod["brand"] and labels.get(other):
                    donor = other
            note = f"no label record for {app}"
            if donor:
                choices = labels[donor]
                note += f"; mechanism quoted from {donor}, the other application of the same brand"
            else:
                note += "; mechanism not quoted"
        else:
            note = ""
        if choices:
            dated = [(ymd(eff), quote) for eff, quote, _is, _brand in choices if eff]
            dated.sort()
            revision, quote = dated[-1]
        else:
            revision, quote = None, ""
        symbol = symbol_of(quote)
        if quote and not symbol:
            raise RuntimeError(f"no HGNC symbol for {app}: {quote}")
        approved = prod["first_approval"] is not None and prod["first_approval"] <= CUTOFF
        rows.append({
            "brand": prod["brand"],
            "generic": prod["generic"],
            "bla": prod["bla"],
            "first_approval": prod["first_approval"].isoformat() if prod["first_approval"] else "",
            "label_revision": revision.isoformat() if revision else "",
            "moa_quote": quote,
            "hgnc_symbol": symbol,
            "approved_by_cutoff": int(approved),
            "search_date": SEARCH_DATE.isoformat(),
            "note": note if note else "",
        })
    brands = {row["brand"].upper() for row in rows}
    missing = [name for name in REQUIRED if name not in brands]
    if missing:
        raise RuntimeError(f"required products missing: {missing}")
    pd.DataFrame(rows).to_csv(OUT / "adc_approved_v2.tsv", sep="\t", index=False)
    pd.DataFrame([{k: item[k] for k in ("endpoint", "query", "n_hits", "note")} for item in logged]).to_csv(
        OUT / "adc_queries.tsv", sep="\t", index=False
    )
    targets = pd.read_csv(ROOT / "config" / "adc_targets.tsv", sep="\t")["gene_symbol"].tolist()
    acting = {}
    for row in rows:
        if row["approved_by_cutoff"] and row["hgnc_symbol"]:
            for symbol in row["hgnc_symbol"].split(";"):
                acting.setdefault(symbol, []).append(row["brand"])
    calls = []
    for gene in targets:
        if gene in acting:
            calls.append({"gene_symbol": gene, "call": "승인 ADC 있음", "products": ";".join(acting[gene]), "note": ""})
        else:
            calls.append({
                "gene_symbol": gene,
                "call": "none found",
                "products": "",
                "note": "no product in this search targets the gene",
            })
    pd.DataFrame(calls).to_csv(OUT / "adc_gene_calls.tsv", sep="\t", index=False)
    n_prod = len({row["brand"] for row in rows if row["approved_by_cutoff"] and row["hgnc_symbol"] in set(targets)})
    n_gene = sum(1 for rec in calls if rec["call"] == "승인 ADC 있음")
    print("products", len(rows), "on29", n_prod, "genes", n_gene, "search", SEARCH_DATE.isoformat(), flush=True)
    print("stage12 adc done", flush=True)


if __name__ == "__main__":
    main()
