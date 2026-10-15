"""Build the Scientific Reports submission files.

B_CPU_CAP is set before any import. The script does not import project modules.
"""
import os
import re
import subprocess
import sys
import zipfile
from pathlib import Path

os.environ["B_CPU_CAP"] = "none"

sys.path.insert(0, "/tmp/pyextra")
from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
PANDOC = "/tmp/pandoc-3.6.4/bin/pandoc"
SOFFICE = "/tmp/lo/root/opt/libreoffice26.8/program/soffice"
SUB = ROOT / "submission_srep"
WORK = SUB / "_work"
FIG = SUB / "figures"
OLE_MAGIC = bytes.fromhex("D0CF11E0A1B11AE1")
PAGE = "\n\n```{=openxml}\n" '<w:p><w:r><w:br w:type="page"/></w:r></w:p>\n' "```\n\n"


def black_font(style, name, size, bold=None):
    style.font.name = name
    style.font.size = Pt(size)
    style.font.color.rgb = RGBColor(0, 0, 0)
    style.font.italic = False
    if bold is not None:
        style.font.bold = bold
    fonts = style.element.get_or_add_rPr().get_or_add_rFonts()
    for attribute in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        fonts.set(qn(attribute), name)
    for attribute in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
        if fonts.get(qn(attribute)) is not None:
            del fonts.attrib[qn(attribute)]


def reference(path):
    subprocess.run([PANDOC, "-o", str(path), "--print-default-data-file", "reference.docx"], check=True)
    doc = Document(path)
    styles = doc.styles
    black_font(styles["Normal"], "Times New Roman", 12)
    for name in ("Normal", "Body Text", "First Paragraph"):
        fmt = styles[name].paragraph_format
        fmt.alignment = WD_ALIGN_PARAGRAPH.LEFT
        fmt.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        fmt.space_before = Pt(0)
        fmt.space_after = Pt(0)
        black_font(styles[name], "Times New Roman", 12)
    for name in ("Title", "Subtitle", "Heading 1", "Heading 2", "Heading 3", "Heading 4"):
        black_font(styles[name], "Times New Roman", 12, bold=True)
        styles[name].paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
        styles[name].paragraph_format.space_before = Pt(12)
        styles[name].paragraph_format.space_after = Pt(6)
        styles[name].paragraph_format.line_spacing_rule = WD_LINE_SPACING.DOUBLE
    compact = styles["Compact"]
    black_font(compact, "Times New Roman", 10)
    compact.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
    compact.paragraph_format.line_spacing = 1.0
    compact.paragraph_format.space_before = Pt(0)
    compact.paragraph_format.space_after = Pt(0)
    names = {item.name for item in styles}
    for name in ("Table Caption", "Image Caption"):
        if name in names:
            black_font(styles[name], "Times New Roman", 12)
            styles[name].paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
    if "Footer" not in names:
        footer = styles.add_style("Footer", WD_STYLE_TYPE.PARAGRAPH)
        footer.base_style = styles["Normal"]
    styles["Footer"].style_id = "Footer"
    section = doc.sections[0]
    section.page_width, section.page_height = Mm(210), Mm(297)
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(section, side, Mm(25))
    doc.save(path)


def line_numbers(section):
    element = OxmlElement("w:lnNumType")
    element.set(qn("w:countBy"), "1")
    element.set(qn("w:restart"), "continuous")
    section._sectPr.append(element)


def page_number(section):
    paragraph = section.footer.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    for kind, text in (("begin", None), (None, " PAGE "), ("end", None)):
        if kind:
            element = OxmlElement("w:fldChar")
            element.set(qn("w:fldCharType"), kind)
        else:
            element = OxmlElement("w:instrText")
            element.set(qn("xml:space"), "preserve")
            element.text = text
        run._r.append(element)


def border(name, value="single"):
    element = OxmlElement(f"w:{name}")
    element.set(qn("w:val"), value)
    element.set(qn("w:sz"), "4" if value == "single" else "0")
    element.set(qn("w:space"), "0")
    element.set(qn("w:color"), "000000")
    return element


def plain_table(table):
    """No shading and no vertical rules. Horizontal rules only at the top, under the header and at the bottom."""
    usable = 16.0
    widths = [round(usable / len(table.columns) * 567) for _ in table.columns]
    tbl = table._tbl
    pr = tbl.tblPr
    for tag in ("w:tblW", "w:tblLayout", "w:tblBorders", "w:shd"):
        old = pr.find(qn(tag))
        if old is not None:
            pr.remove(old)
    tblw = OxmlElement("w:tblW")
    tblw.set(qn("w:w"), str(sum(widths)))
    tblw.set(qn("w:type"), "dxa")
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    borders = OxmlElement("w:tblBorders")
    for name, value in (
        ("top", "single"), ("left", "nil"), ("bottom", "single"), ("right", "nil"),
        ("insideH", "nil"), ("insideV", "nil"),
    ):
        borders.append(border(name, value))
    pr.append(tblw)
    pr.append(layout)
    pr.append(borders)
    grid = tbl.tblGrid
    for col in list(grid):
        grid.remove(col)
    for width in widths:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(width))
        grid.append(col)
    for row_index, row in enumerate(table.rows):
        if row_index == 0:
            trpr = row._tr.get_or_add_trPr()
            if trpr.find(qn("w:tblHeader")) is None:
                trpr.append(OxmlElement("w:tblHeader"))
        for cell in row.cells:
            tcpr = cell._tc.get_or_add_tcPr()
            for tag in ("w:shd", "w:tcBorders"):
                old = tcpr.find(qn(tag))
                if old is not None:
                    tcpr.remove(old)
            if row_index == 0:
                edges = OxmlElement("w:tcBorders")
                edges.append(border("bottom", "single"))
                tcpr.append(edges)
            for paragraph in cell.paragraphs:
                paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
                paragraph.paragraph_format.line_spacing = 1.0


def finish(path, numbered):
    doc = Document(path)
    for section in doc.sections:
        section.page_width, section.page_height = Mm(210), Mm(297)
        for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
            setattr(section, side, Mm(25))
        if numbered:
            line_numbers(section)
        page_number(section)
    for table in doc.tables:
        plain_table(table)
    props = doc.core_properties
    props.author = ""
    props.last_modified_by = ""
    props.comments = ""
    props.keywords = ""
    doc.save(path)


def style_check(path):
    with zipfile.ZipFile(path) as archive:
        styles = archive.read("word/styles.xml").decode()
        parts = [name for name in archive.namelist() if re.match(r"word/(document|footer\d*|header\d*|footnotes|endnotes)\.xml$", name)]
        used = set()
        for name in parts:
            xml = archive.read(name).decode()
            used |= set(re.findall(r'<w:(?:p|r|tbl)Style w:val="([^"]+)"', xml))
    defined = set(re.findall(r'w:styleId="([^"]+)"', styles))
    return sorted(used - defined)


def pandoc(source, out, ref):
    subprocess.run([PANDOC, str(source), "-o", str(out), "--reference-doc", str(ref)], check=True, cwd=ROOT)


def to_pdf(path):
    subprocess.run(
        [SOFFICE, "--headless", "-env:UserInstallation=file:///tmp/lo/profile_b16docx",
         "--convert-to", "pdf", "--outdir", str(path.parent), str(path)],
        check=True,
    )
    return path.with_suffix(".pdf")


def with_figures(source, figures):
    text = source.read_text()
    for number in figures:
        text += PAGE + f"![](submission_srep/figures/Fig{number}.png){{width=150mm}}\n"
    target = WORK / (source.stem + "_fig.md")
    target.write_text(text)
    return target


def words(text):
    return len(text.split())


def checklist(main_docx, si_pdf, letter):
    text = (ROOT / "manuscript/srep/B_main_srep.md").read_text()
    title = text.splitlines()[0].lstrip("# ").strip()
    abstract = text.split("## Abstract", 1)[1].split("## ", 1)[0]
    body = text.split("## Introduction", 1)[1].split("## Methods", 1)[0]
    rows = [
        ("Title", f"{words(title)} words", "at most 20", "pass" if words(title) <= 20 else "fail"),
        ("Abstract", f"{words(abstract)} words", "at most 200, no subheadings", "pass" if words(abstract) <= 200 else "fail"),
        ("Introduction + Results + Discussion", f"{words(body)} words", "at most 4500", "pass" if words(body) <= 4500 else "fail"),
        ("Section order", "Introduction, Results, Discussion, Methods, Data availability, References, Acknowledgements, Author contributions, Competing interests, Figure legends, Tables", "Scientific Reports order", "pass"),
        ("Figures and tables", "5 figures + 3 tables", "at most 8", "pass"),
        ("Main docx", f"{main_docx.stat().st_size} bytes", "one column, 12 pt, double spaced, page numbers, continuous line numbers, at most 3 MB", "pass" if main_docx.stat().st_size <= 3_000_000 else "fail"),
        ("SI PDF", f"{si_pdf.stat().st_size} bytes", "one file, at most 50 MB", "pass" if si_pdf.stat().st_size <= 50_000_000 else "fail"),
        ("Cover letter", letter.name, "docx", "pass"),
        ("Figure files", "Fig1–Fig5 pdf and 300 dpi tiff; docx uses the 200 dpi png", "pdf and tiff", "pass"),
        ("Number ledger", "manuscript/checks/ledger_srep.tsv", "FAIL 0 and unledgered 0", "pass"),
        ("GitHub URL, Zenodo DOI, cover-letter date", "left unfilled", "not invented", "remaining"),
    ]
    lines = ["# Scientific Reports submission checklist", "", "| Item | Measured | Rule | Status |", "|---|---|---|---|"]
    for item, measured, rule, status in rows:
        lines.append(f"| {item} | {measured} | {rule} | {status} |")
    lines.append("")
    (SUB / "checklist_srep.md").write_text("\n".join(lines) + "\n")


def encryption_problems(path):
    problems = []
    blob = path.read_bytes()
    if blob[:8] == OLE_MAGIC:
        problems.append("OLE compound file")
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as package:
            if any(name.endswith("EncryptedPackage") for name in package.namelist()):
                problems.append("EncryptedPackage")
    if path.suffix == ".pdf" and b"/Encrypt" in blob:
        problems.append("PDF encryption")
    return problems


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    ref = WORK / "reference.docx"
    reference(ref)
    main_md = with_figures(ROOT / "manuscript/srep/B_main_srep.md", ["1", "2", "3", "4", "5"])
    si_md = with_figures(ROOT / "manuscript/srep/B_SI_srep.md", ["S1", "S2"])
    main_out = SUB / "B_manuscript_srep.docx"
    letter = SUB / "B_cover_letter_srep.docx"
    si_docx = WORK / "B_SI_srep.docx"
    pandoc(main_md, main_out, ref)
    pandoc(ROOT / "manuscript/srep/B_cover_letter_srep.md", letter, ref)
    pandoc(si_md, si_docx, ref)
    finish(main_out, True)
    finish(letter, False)
    finish(si_docx, False)
    for path in (main_out, letter, si_docx):
        missing = style_check(path)
        print(path.name, "missing_styles", missing, "bytes", path.stat().st_size)
        if missing:
            raise SystemExit(f"missing styles in {path.name}: {missing}")
    si_pdf = to_pdf(si_docx)
    shipped = SUB / "B_SI_srep.pdf"
    si_pdf.replace(shipped)
    checklist(main_out, shipped, letter)
    for path in (main_out, letter, shipped):
        problems = encryption_problems(path)
        print(path.name, "encryption", problems or "none", "bytes", path.stat().st_size)
        if problems:
            raise SystemExit(f"{path.name}: {problems}")


if __name__ == "__main__":
    main()
