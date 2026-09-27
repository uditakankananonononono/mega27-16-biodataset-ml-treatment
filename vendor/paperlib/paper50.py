"""paperlib.paper50 - docx builder for MEGA27 papers (reconstructed 2026-09-27).
Times New Roman body text; headings bold; equations italic-centered.
"""
import docx
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.shared import Pt, Inches


def _style_run(r, size=12, bold=False, italic=False, font="Times New Roman"):
    r.font.name = font
    r.font.size = Pt(size)
    r.font.bold = bold
    r.font.italic = italic


def new_doc():
    doc = docx.Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    cp = doc.core_properties
    cp.author = ""
    cp.last_modified_by = ""
    return doc


def title_block(doc, title, subtitle):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _style_run(p.add_run(title), size=20, bold=True)
    p2 = doc.add_paragraph()
    p2.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _style_run(p2.add_run(subtitle), size=12, italic=True)


def h1(doc, text):
    p = doc.add_paragraph()
    _style_run(p.add_run(text), size=16, bold=True)


def h2(doc, text):
    p = doc.add_paragraph()
    _style_run(p.add_run(text), size=13, bold=True)


def para(doc, text):
    p = doc.add_paragraph()
    _style_run(p.add_run(text))


def eq(doc, num, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _style_run(p.add_run(f"({num})   {text}"), italic=True)


def table(doc, caption, header, rows, size=10):
    p = doc.add_paragraph()
    _style_run(p.add_run(caption), bold=True)
    t = doc.add_table(rows=1 + len(rows), cols=len(header))
    t.style = "Table Grid"
    t.autofit = False
    w = Inches(6.5 / len(header))
    for r_ in t.rows:
        for c_ in r_.cells:
            c_.width = w
    for j, h in enumerate(header):
        c = t.rows[0].cells[j]
        _style_run(c.paragraphs[0].add_run(str(h)), size=size, bold=True)
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            c = t.rows[i + 1].cells[j]
            _style_run(c.paragraphs[0].add_run(str(v)), size=size)


def figure(doc, path, caption):
    try:
        doc.add_picture(path, width=Inches(5.5))
    except Exception:
        para(doc, f"[figure missing: {path}]")
    p = doc.add_paragraph()
    _style_run(p.add_run(caption), size=10, italic=True)


def page_break(doc):
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)


def save(doc, path):
    doc.core_properties.author = ""
    doc.core_properties.last_modified_by = ""
    doc.save(path)
