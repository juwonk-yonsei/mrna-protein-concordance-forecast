"""Record Crossref checks. Writes manuscript/checks only after the queries return."""
import json
import os
import urllib.request
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"

ROOT = Path(__file__).resolve().parents[1]
DOIS = [
    "10.1016/j.cels.2020.06.013",
    "10.1016/j.crmeth.2022.100288",
    "10.1016/j.xcrm.2022.100819",
    "10.1016/j.cell.2024.05.039",
    "10.1016/j.ejca.2023.113379",
    "10.1016/j.cels.2023.07.007",
]


def main():
    rows = ["doi\tyear\tfirst_author\ttitle\tstatus"]
    for doi in DOIS:
        with urllib.request.urlopen("https://api.crossref.org/works/" + doi, timeout=60) as handle:
            msg = json.load(handle)["message"]
        author = (msg.get("author") or [{}])[0]
        year = (msg.get("issued", {}).get("date-parts") or [[""]])[0][0]
        title = " ".join(msg.get("title") or [])
        rows.append(f"{doi}\t{year}\t{author.get('family','')}\t{title}\tok")
    out = ROOT / "manuscript" / "checks"
    out.mkdir(parents=True, exist_ok=True)
    (out / "crossref.tsv").write_text("\n".join(rows) + "\n")
    print("crossref", len(DOIS), flush=True)


if __name__ == "__main__":
    main()
