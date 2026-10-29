#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""What is set in italic, in one place: gene symbols and organism names.

Nothing in this pipeline was italicised until now - not the binomials and not
the gene symbols - because the document builder writes plain runs. The rule is
the ordinary one:

  * a GENE symbol is italic (mecA, erm(A), ausA, the dlt operon);
  * the PROTEIN of the same name is roman (MecA, AusA, Lnu(A), YcaO, LanB),
    which is why the vocabulary is an explicit list and not a pattern: a rule
    like "lower case then a capital" cannot tell lugA from LugA reliably once
    a sentence starts with the word;
  * an organism name is italic, genus alone or genus with epithet, spelled out
    or abbreviated.

GENES is closed on purpose. Adding a gene to the paper means adding it here,
and the build fails loudly rather than shipping one roman symbol among fifty
italic ones: check_italics() is the gate.
"""
import re

# Gene symbols used anywhere in this paper: text, legends, tables, figures.
# Order does not matter here; the pattern is built longest-first below.
GENES = [
    # biosynthesis
    "ausA", "ausB", "lugA", "lugB", "lugC", "lugD",
    "lanA", "lanB", "lanC", "lanM", "lanKC", "ycaO",
    # resistance: every determinant symbol the AMRFinderPlus reference set
    # returns on this panel (out/amr_located.csv), whether or not it is printed
    # today, so that adding a row to a table cannot ship a roman symbol.
    # "str" is deliberately absent: it occurs on no genome here and the
    # three-letter word is too easy to match inside ordinary prose.
    "mecA", "mecC", "blaZ", "blaI", "blaR1",
    "erm(A)", "erm(C)", "erm(33)", "erm(44)", "erm(44)v", "erm(45)",
    "lnu(A)", "lnu(B)", "lsa(B)", "lsa(E)", "msr(A)", "mph(C)",
    "tet(K)", "tet(L)", "tet(M)", "tet(S)", "vat(B)", "vga(A)", "vga(B)",
    "aac(6')-Ie", "aac(6')-aph(2'')", "aph(2'')-Ia", "aph(3')-Ia",
    "aph(3')-IIIa", "ant(4')-Ia", "ant(6)-Ia", "ant(9)-Ia",
    # the gene maps label an aminoglycoside gene whose allele the annotation
    # does not give with the bare family symbol
    "aph", "ant",
    "catA7", "catA8", "cfr", "dfrG", "dfrK", "fexA", "fusB", "fusD",
    "mupA", "optrA", "qacA", "qacB", "qacC", "qacG", "qacJ",
    "cadA", "arsC", "chrA",
    # other loci named in the text
    "saeS", "dlt",
]

# Species epithets standing on their own. A figure has no room for the genus
# on every tick, so "F002 (capitis)" and a tree leaf reading "epidermidis" are
# the organism name just as much as the binomial is, and are set the same way.
EPITHETS = [
    "epidermidis", "haemolyticus", "capitis", "lugdunensis", "hominis",
    "saprophyticus", "simulans", "warneri", "pasteuri", "cohnii", "caprae",
    "pettenkoferi", "aureus", "sciuri", "xylosus", "succinus",
]

# Genera that appear in this paper. A binomial is the genus plus a lower-case
# epithet; the abbreviated form is the initial, a full stop and the epithet.
GENERA = ["Staphylococcus", "Mammaliicoccus", "Enterococcus", "Streptococcus",
          "Lactococcus", "Lactobacillus", "Limosilactobacillus",
          "Cutibacterium", "Micrococcus", "Escherichia", "Listeria"]

_ABBR = sorted({g[0] for g in GENERA})

# An axis with twelve species on it has room for "epid.", not for
# "S. epidermidis"; a truncation of a name is still the name.
_TRUNC = sorted({e[:4] for e in EPITHETS})

# A gene symbol must stand alone: without these guards "ant" matches inside
# "significant" and "aph" inside "graph".
_EDGE_L, _EDGE_R = r"(?<![A-Za-z0-9_])", r"(?![A-Za-z0-9_])"

_PARTS = (
    # binomials first, so "Staphylococcus epidermidis" is one run, not two,
    # and the bare epithet below never fires inside one
    [rf"\b{g}\s+[a-z]{{3,}}\b" for g in GENERA]
    + [rf"\b[{''.join(_ABBR)}]\.\s[a-z]{{3,}}\b"]
    + [rf"\b{g}\b" for g in GENERA]
    + [_EDGE_L + "(?:" + "|".join(EPITHETS) + ")" + _EDGE_R]
    + [_EDGE_L + "(?:" + "|".join(_TRUNC) + r")\."]
    # gene symbols, longest first so aac(6')-Ie wins over aac(6')
    + [_EDGE_L + "(?:"
       + "|".join(re.escape(g) for g in sorted(GENES, key=len, reverse=True))
       + ")" + _EDGE_R]
)
TOKEN = re.compile("(" + "|".join(_PARTS) + ")")

# "SCCmec" is one word whose mec part is italic; it is handled separately
# because the boundary falls inside the word.
SCCMEC = re.compile(r"\bSCC(mec)\b")


def split(text):
    """[(chunk, italic), ...] covering text exactly once, in order."""
    out = []

    def emit(s, ital):
        if not s:
            return
        if out and out[-1][1] == ital:
            out[-1] = (out[-1][0] + s, ital)
        else:
            out.append((s, ital))

    def roman(s):
        # SCCmec is the only case where the italic part is inside a word
        last = 0
        for m in SCCMEC.finditer(s):
            emit(s[last:m.start()], False)
            emit("SCC", False)
            emit("mec", True)
            last = m.end()
        emit(s[last:], False)

    last = 0
    for m in TOKEN.finditer(text):
        roman(text[last:m.start()])
        emit(m.group(0), True)
        last = m.end()
    roman(text[last:])
    return out


# ------------------------------------------------------------------ docx --
def italicise_docx(path):
    """Re-run every paragraph and table cell of a .docx through split().

    Post-processing the built file rather than teaching each builder about
    italics keeps one implementation for the manuscript, the figure book and
    the submission package, and keeps the shared docbuild module untouched.
    Returns the number of runs that were split.
    """
    import copy
    from docx import Document
    from docx.oxml.ns import qn

    doc = Document(str(path))

    def paragraphs(container):
        yield from container.paragraphs
        for t in container.tables:
            for row in t.rows:
                for cell in row.cells:
                    yield from paragraphs(cell)

    n = 0
    for p in list(paragraphs(doc)):
        for run in list(p.runs):
            parts = split(run.text)
            if len(parts) < 2 and not (parts and parts[0][1]):
                continue
            n += 1
            el = run._r
            anchor = el
            for i, (chunk, ital) in enumerate(parts):
                if i == 0:
                    run.text = chunk
                    run.font.italic = True if ital else run.font.italic
                    if not ital:
                        run.font.italic = False
                    continue
                new = copy.deepcopy(el)
                for t in new.findall(qn("w:t")):
                    new.remove(t)
                t = new.makeelement(qn("w:t"), {})
                t.text = chunk
                t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
                new.append(t)
                rPr = new.find(qn("w:rPr"))
                if rPr is None:
                    rPr = new.makeelement(qn("w:rPr"), {})
                    new.insert(0, rPr)
                for tag in ("w:i", "w:iCs"):
                    for old in rPr.findall(qn(tag)):
                        rPr.remove(old)
                if ital:
                    for tag in ("w:i", "w:iCs"):
                        rPr.append(rPr.makeelement(qn(tag), {}))
                anchor.addnext(new)
                anchor = new
    doc.save(str(path))
    return n


def roman_everything(path):
    """Clear every italic run, so italicise_docx() decides italics alone.

    docbuild sets a figure or table note entirely in italic. That is the house
    look for the internal figure book, but it hides the thing this module
    exists for: inside an all-italic legend a gene symbol is not distinguishable
    from the words around it, and the same legend appears roman in the
    manuscript's own Figure legends section, so one file carried it both ways.
    Archives of Microbiology prints roman legends. Nothing else in these
    documents is italic by design, so clearing and re-deciding is safe.
    """
    from docx import Document
    doc = Document(str(path))

    def paragraphs(container):
        yield from container.paragraphs
        for t in container.tables:
            for row in t.rows:
                for cell in row.cells:
                    yield from paragraphs(cell)

    n = 0
    for p in paragraphs(doc):
        for run in p.runs:
            if run.font.italic:
                run.font.italic = False
                n += 1
    doc.save(str(path))
    return n


def check_italics(path, allow=()):
    """Fail if any gene symbol or binomial in the file is still roman.

    This is what stops a new gene symbol from shipping unitalicised: it is not
    enough to have italicised the ones that were there when the pass was
    written.
    """
    from docx import Document
    doc = Document(str(path))

    def paragraphs(container):
        yield from container.paragraphs
        for t in container.tables:
            for row in t.rows:
                for cell in row.cells:
                    yield from paragraphs(cell)

    bad = []
    for p in paragraphs(doc):
        for run in p.runs:
            if run.font.italic:
                continue
            for chunk, ital in split(run.text):
                if ital and chunk not in allow:
                    bad.append(chunk)
    return sorted(set(bad))


# --------------------------------------------------------------- figures --
def _math_escape(chunk):
    out = chunk.replace("\\", r"\\")
    for ch in "$^_{}#&%":
        out = out.replace(ch, "\\" + ch)
    return out


def mpl(text):
    """(label, base_italic) for matplotlib: the gene and organism parts italic.

    Two ways round, because mathtext will not italicise a bracket or a full
    stop: inside `$\\it{...}$` the letters come from the italic font but the
    punctuation always comes from the roman one, so `erm(A)` would print with
    upright brackets and `epid.` with an upright stop while the same names are
    fully italic in the text and the tables.

      * ordinary case - the italic pieces are plain words - each goes in its
        own math group and the label is drawn in the normal upright font;
      * a name carrying punctuation - erm(A), aac(6')-Ie, epid., S. capitis -
        the whole label is drawn in the *italic* font and the roman pieces go
        in math groups instead, so the punctuation is italic with the name.

    The caller sets the font style from the second return value.  A label with
    more than one line is handled a line at a time: a math group cannot
    contain a newline.
    """
    parts = split(text)
    if not any(ital for _, ital in parts):
        return text, False
    flip = any(ital and not re.fullmatch(r"[A-Za-z0-9 ]+", chunk)
               for chunk, ital in parts)

    # 進 math group 的只有英數字與括號類：mathtext 會照數學規則排版，
    # "n=293" 變成 "n = 293"、連字號變成減號，那是把作者寫的數字改掉。
    # 其餘標點（空白、句點、逗號、冒號、等號、連字號）留在底層字型：
    # 斜體的句點與等號跟正體看不出差別，位置與字元則一個都沒動。
    MATHY = re.compile(r"[A-Za-z0-9()\[\]{}/%<>]+")

    def groups(chunk, cmd):
        """把 chunk 切成 (字串, 要不要進 math group) 的序列。"""
        out, last = [], 0
        for m in MATHY.finditer(chunk):
            if m.start() > last:
                out.append((chunk[last:m.start()], False))
            out.append((m.group(0), True))
            last = m.end()
        if last < len(chunk):
            out.append((chunk[last:], False))
        return "".join(rf"${cmd}{{{_math_escape(part)}}}$" if math else part
                       for part, math in out)

    def render(line):
        out = []
        for chunk, ital in split(line):
            if ital != flip:          # the piece that goes in math groups
                out.append(groups(chunk, r"\mathrm" if flip else r"\it"))
            else:                     # the piece the base font already covers
                out.append(chunk.replace("\\", r"\\").replace("$", r"\$"))
        return "".join(out)

    out = "\n".join(render(line) for line in text.split("\n"))
    # 這一段自己驗一次：把 markup 拆掉要還原成原字串。少一個字或多一個句點
    # （"S. epid." 曾經印成 "S. epid.."）在圖上不會有人看出來，讓它在這裡爆。
    assert plain(out) == text, (out, text)
    return out, flip


_MARKUP = re.compile(r"\$\\(?:it|mathrm)\{(.*?)\}\$")


def plain(label):
    """把 mpl() 產生的 mathtext markup 拆回原字串。"""
    out = _MARKUP.sub(lambda m: m.group(1), label)
    for ch in "$^_{}#&%":
        out = out.replace("\\" + ch, ch)
    return out.replace(r"\\", "\\")


def force(text):
    """Mathtext for a string that is italic in full, vocabulary or not.

    An axis with twelve species on it is labelled "epid." rather than
    "S. epidermidis", and a truncation is not a name the vocabulary can list.
    The caller states that this piece is an organism name; everything it is
    concatenated with stays roman.
    """
    esc = text.replace("\\", r"\\")
    for ch in "$^_{}#&%":
        esc = esc.replace(ch, "\\" + ch)
    return " ".join(rf"$\it{{{w}}}$" for w in esc.split(" ") if w)


if __name__ == "__main__":
    for s in ("Carriage of mecA in S. epidermidis and erm(A) in ST23",
              "the dlt and biotin operons", "SCCmec type IV carries mecA",
              "AusA is the protein, ausA the gene",
              "aac(6')-Ie and ant(9)-Ia", "Staphylococcus lugdunensis makes lugdunin"):
        print(f"{s}\n  -> {split(s)}\n  -> {mpl(s)}\n")
