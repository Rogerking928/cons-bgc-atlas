"""Project-local Word builder, replacing ecgfig.build_doc for this figure book.

The shared builder is fine for its own papers and is left alone. This one fixes
the four things that made the v4 document unsubmittable, and it fixes them
structurally rather than by tuning font sizes until a particular page happened
to fit.

**Tables that span pages repeat their header.** A reader who turns to the second
page of Table 1 currently sees six unlabelled columns of numbers. The header row
is marked with w:tblHeader so Word repeats it.

**Notes do not strand.** A table that nearly fills a page pushes its footnote
alone onto the next one, which reads as a blank page with two lines on it. The
last row of every table and every caption carries w:keepNext, so the note
travels with the object it belongs to.

**Nothing is wider than the text column.** Images are capped well inside the
column width rather than at it, and tables are given an explicit fixed layout at
the column width, so no object can push a line past the margin.

**Height is checked, not hoped for.** Every item's height is estimated before it
is written — caption lines, note lines, image aspect ratio, table rows — and
`build` returns the estimates so the caller can see what will overflow instead
of discovering it in the PDF.

Geometry: landscape Letter, 1.0 inch left and right, 0.9 top and bottom. That
gives a 9.0 by 6.7 inch text block, half an inch more height than the shared
builder allowed, which is most of what the figures needed.
"""
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

PAGE_W, PAGE_H = 11.0, 8.5
MARGIN_LR, MARGIN_TB = 1.0, 0.9
TEXT_W = PAGE_W - 2 * MARGIN_LR          # 9.0
TEXT_H = PAGE_H - 2 * MARGIN_TB          # 6.7
IMG_MAX_W = 8.2                          # kept clear of the column edge
CHARS_PER_INCH = 15.5                    # 10 pt Times, rough but stable


# ECMA-376 fixes the order of w:pPr children. Only the ones this module writes
# need to be listed; anything unknown sorts to the end, which is where
# python-docx's own late additions belong anyway.
_PPR_ORDER = ("w:pStyle", "w:keepNext", "w:keepLines", "w:pageBreakBefore",
              "w:framePr", "w:widowControl", "w:numPr", "w:pBdr", "w:shd",
              "w:tabs", "w:spacing", "w:ind", "w:contextualSpacing", "w:jc",
              "w:textAlignment", "w:outlineLvl", "w:rPr", "w:sectPr")


def _prepend_pPr(par, tag):
    """Insert a paragraph property and leave w:pPr in schema order.

    Inserting at the front is right for keepNext but wrong for
    pageBreakBefore, which must follow keepLines; a caption that got both ended
    up with pageBreakBefore first. Renderers are inconsistent about invalid
    ordering, and the symptom here was the one the ordering rule exists to
    prevent: on exactly the pages that carry a page break, Word dropped the
    paragraph indentation and clipped the first characters of the caption
    ("Table 2" printed as "e 2"). Sorting the whole element after every
    insertion makes the order correct however the callers are sequenced.
    """
    pPr = par._p.get_or_add_pPr()
    el = OxmlElement(tag)
    pPr.append(el)
    rank = {qn(t): i for i, t in enumerate(_PPR_ORDER)}
    for child in sorted(list(pPr), key=lambda c: rank.get(c.tag, len(_PPR_ORDER))):
        pPr.append(child)          # re-appending an existing child moves it
    return el


def _keep_next(par):
    """Ask Word to keep this paragraph on the same page as the next one."""
    for tag in ("w:keepLines", "w:keepNext"):   # inserted at 0, so reverse order
        _prepend_pPr(par, tag)


def _page_break_before(par):
    """Start this paragraph on a new page, as a property rather than a break.

    An inserted page-break paragraph is its own paragraph: it ends the previous
    page and then occupies a line at the top of the next one. Worse, when the
    preceding item already ran to the bottom of its page, the break lands on a
    page that is already new and produces an empty one. `w:pageBreakBefore` is
    a property of the paragraph that follows, so Word starts a new page only if
    it is not already at the top of one, and no empty paragraph exists to leave
    a stray line behind.
    """
    _prepend_pPr(par, "w:pageBreakBefore")


def _strip_empty_paragraphs(doc):
    """Remove paragraphs with neither text nor a picture.

    python-docx's default template opens with an empty paragraph, and an empty
    paragraph at a page boundary is exactly what a blank page is made of. A
    paragraph holding a figure has no text, so the drawing check has to come
    before the text check.
    """
    for p in list(doc.paragraphs):
        if p.text.strip():
            continue
        if p._p.findall(".//" + qn("w:drawing")):
            continue
        if p._p.findall(".//" + qn("w:sectPr")):
            continue        # 這一段空白段落帶著前一節的版面設定，刪掉整節就沒了
        parent = p._element.getparent()
        if parent is not None:
            parent.remove(p._element)


def _repeat_header(row):
    """Mark a table row as a header that repeats on every continuation page."""
    trPr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    trPr.append(el)


def _col_widths(item, ncol):
    """Fractions of the text column each table column gets.

    Equal columns waste space: the label column carries sentences and the rest
    carry six characters. Cramming the labels makes every row wrap to two lines,
    which is what turns a twenty-row table into a two-page one.
    """
    cw = item.get("col_widths")
    if cw:
        # A caller-supplied list that does not match the table is a bug in the
        # caller, but it must not take the build down: pad or trim, then
        # renormalise, so the document still renders and the mismatch shows up
        # as a slightly odd column rather than a traceback.
        cw = list(cw)[:ncol]
        if len(cw) < ncol:
            cw += [cw[-1] if cw else 1.0] * (ncol - len(cw))
        tot = float(sum(cw)) or 1.0
        return [c / tot for c in cw]
    if ncol <= 1:
        # a one-column table is a section divider, not a table; it takes the
        # whole width rather than dividing by a column count of zero
        return [1.0]
    first = 0.36 if ncol >= 5 else 0.30
    return [first] + [(1 - first) / (ncol - 1)] * (ncol - 1)


def _fixed_layout(table, width_in):
    """Fix the table width so a long cell cannot push it past the margin."""
    table.autofit = False
    tblPr = table._tbl.tblPr
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    tblPr.append(layout)
    w = OxmlElement("w:tblW")
    w.set(qn("w:w"), str(int(width_in * 1440)))
    w.set(qn("w:type"), "dxa")
    tblPr.append(w)


def _wrapped_lines(text, size_pt, width_in):
    """How many lines a run of text takes at this size in this column."""
    if not text:
        return 0
    per_line = max(20, int(width_in * CHARS_PER_INCH * (10.0 / size_pt)))
    return max(1, -(-len(str(text)) // per_line))


def estimate_height(item):
    """Inches the item will occupy, caption and note included."""
    cap = _wrapped_lines(item.get("caption", ""), 12, TEXT_W) * 0.22
    note = _wrapped_lines(item.get("note", ""), 9.5, TEXT_W) * 0.17
    if item["kind"] == "figure":
        from PIL import Image
        with Image.open(item["path"]) as im:
            px_w, px_h = im.size
        body = fit_figure_width(item) * px_h / px_w
    else:
        rows = item["rows"]
        size = item.get("body_size", 9.5)
        per_row = size / 72.0 * 1.45 + 0.03
        body = 0.0
        fr = _col_widths(item, len(rows[0]))
        for r in rows:
            lines = 1
            for ci, c in enumerate(r):
                lines = max(lines, str(c).count("\n") + 1,
                            _wrapped_lines(c, size, fr[ci] * TEXT_W))
            body += lines * per_row
    return cap + note + body + 0.12


def fit_figure_width(item):
    """Largest width that keeps caption, note and image on one page."""
    from PIL import Image
    with Image.open(item["path"]) as im:
        px_w, px_h = im.size
    ratio = px_h / px_w
    cap = _wrapped_lines(item.get("caption", ""), 12, TEXT_W) * 0.22
    note = _wrapped_lines(item.get("note", ""), 9.5, TEXT_W) * 0.17
    avail = TEXT_H - cap - note - 0.16
    return max(3.0, min(item.get("width", IMG_MAX_W), IMG_MAX_W, avail / ratio))


def build(items, out_path, font="Times New Roman", doc=None):
    """items: {"kind": "figure"|"table", "caption", "note", ...} as before.

    `doc` appends the items to an existing document instead of starting a new
    one, in a landscape section of their own, so a portrait manuscript and its
    landscape figure book can be delivered as a single file. The paragraphs
    already in that document must carry explicit run sizes, because the Normal
    style is set here for the items.
    """
    if doc is None:
        doc = Document()
        sec = doc.sections[0]
    else:
        sec = doc.add_section(WD_SECTION.NEW_PAGE)
    n = doc.styles["Normal"]
    n.font.name = font
    n.font.size = Pt(11)
    n.element.rPr.rFonts.set(qn("w:eastAsia"), font)
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = Inches(PAGE_W), Inches(PAGE_H)
    sec.left_margin = sec.right_margin = Inches(MARGIN_LR)
    sec.top_margin = sec.bottom_margin = Inches(MARGIN_TB)

    def run(r, size=11, bold=False, italic=False):
        r.font.name = font
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.italic = italic
        r._element.rPr.rFonts.set(qn("w:eastAsia"), font)
        return r

    def para(align=WD_ALIGN_PARAGRAPH.LEFT, sa=4, sb=0):
        p = doc.add_paragraph()
        p.alignment = align
        pf = p.paragraph_format
        pf.space_after = Pt(sa)
        pf.space_before = Pt(sb)
        pf.left_indent = Inches(0)        # explicit: nothing hangs off the margin
        pf.right_indent = Inches(0)
        pf.first_line_indent = Inches(0)
        return p

    report = []
    for i, it in enumerate(items):
        cap = para(sa=4)
        run(cap.add_run(it["caption"]), size=11.5, bold=True)
        _keep_next(cap)
        if i:
            _page_break_before(cap)

        if it["kind"] == "figure":
            if it.get("note"):
                p = para(sa=4)
                run(p.add_run(it["note"]), size=9.5, italic=True)
                _keep_next(p)
            w = fit_figure_width(it)
            p = para(align=WD_ALIGN_PARAGRAPH.CENTER, sa=0, sb=2)
            p.add_run().add_picture(str(it["path"]), width=Inches(w))
            report.append((it["caption"][:38], round(w, 2),
                           round(estimate_height(it), 2)))
        else:
            rows = it["rows"]
            size = it.get("body_size", 9.5)
            # A table taller than the page will continue overleaf. Word can
            # repeat header rows but has no notion of "Table 1 continued", so
            # the caption is moved into the table as a merged first row and
            # marked as a header alongside the column labels. It then repeats
            # verbatim on every page the table occupies, and the standalone
            # caption paragraph above it is dropped to avoid printing it twice.
            spans = estimate_height(it) > TEXT_H
            if spans:
                cap._p.getparent().remove(cap._p)
            bold_rows = set(it.get("bold_rows", [0]))
            n_head = 1 + int(spans)
            t = doc.add_table(rows=len(rows) + int(spans), cols=len(rows[0]))
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            _fixed_layout(t, TEXT_W)
            if spans:
                top = t.rows[0].cells[0]
                for c in t.rows[0].cells[1:]:
                    top = top.merge(c)
                cp = top.paragraphs[0]
                cp.alignment = WD_ALIGN_PARAGRAPH.LEFT
                run(cp.add_run(it["caption"]), size=11.5, bold=True)
                rows = [list(r) for r in rows]
                _repeat_header(t.rows[0])
            for r_ in t.rows[:n_head]:
                _repeat_header(r_)
            fr = _col_widths(it, len(rows[0]))
            for ci, f_ in enumerate(fr):
                for row in t.rows[int(spans):]:
                    row.cells[ci].width = Inches(f_ * TEXT_W)
            for ri, r_ in enumerate(rows):
                for ci, val in enumerate(r_):
                    cp = t.cell(ri + int(spans), ci).paragraphs[0]
                    cp.alignment = (WD_ALIGN_PARAGRAPH.LEFT if ci == 0
                                    else WD_ALIGN_PARAGRAPH.CENTER)
                    cp.paragraph_format.space_after = Pt(0)
                    cp.paragraph_format.space_before = Pt(0)
                    lines = str(val).split("\n")
                    rr = cp.add_run(lines[0])
                    for extra in lines[1:]:
                        rr.add_break()
                        rr.add_text(extra)
                    run(rr, size=size, bold=(ri in bold_rows))
                    if ri == len(rows) - 1 and it.get("note"):
                        _keep_next(cp)      # the note follows the last row
            if it.get("note"):
                p = para(sb=3, sa=0)
                run(p.add_run(it["note"]), size=9.0, italic=True)
            report.append((it["caption"][:38], "-",
                           round(estimate_height(it), 2)))

    _strip_empty_paragraphs(doc)
    out_path = Path(out_path)
    doc.save(out_path)
    ntab = sum(1 for it in items if it["kind"] == "table")
    print(f"saved {out_path}  ({ntab} tables + {len(items)-ntab} figures)",
          flush=True)
    print(f"\n{'item':<40}{'width':>7}{'est. height':>13}   "
          f"(text block {TEXT_W} x {TEXT_H} in)", flush=True)
    for name, w, h in report:
        flag = "  OVERFLOWS" if h > TEXT_H else ""
        print(f"{name:<40}{str(w):>7}{h:>13.2f}{flag}", flush=True)
    return out_path
