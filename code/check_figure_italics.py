#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gate: every gene symbol and organism name in a figure is set in italic.

The .docx files have had such a gate since the italic pass was written
(italics.check_italics), but the figures had none: the vocabulary was applied
on the way out of matplotlib and nobody read the result back. Two faults hid
there until this was written.

  * a tick label is re-written from the axis formatter at every draw, so the
    italic text set on the label object was thrown away and Figure 4 printed
    F001-ausA in roman throughout;
  * a label set italic in full - style="italic" on a whole tick label - made
    "n=293" and "other (4 spp.)" italic as well, which says they are names.

Reading the rendered PDF back is the only check that sees both, because it
asks what the reader will see rather than what the script asked for.

Run: python3 code/check_figure_italics.py [pdf ...]
"""
import glob
import os
import re
import sys

import fitz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import italics

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "out")

WORD = re.compile(r"[A-Za-z0-9]")


def lines(path):
    """(text, [italic per character]) for every rendered line of a PDF."""
    for page in fitz.open(path):
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                text, flags = "", []
                for span in line["spans"]:
                    italic = "Italic" in span["font"] or "Oblique" in span["font"]
                    text += span["text"]
                    flags += [italic] * len(span["text"])
                if text.strip():
                    yield text, flags


def audit(path):
    roman_names, italic_others = [], []
    for text, flags in lines(path):
        i = 0
        for chunk, want_italic in italics.split(text):
            got = flags[i:i + len(chunk)]
            i += len(chunk)
            letters = [f for f, c in zip(got, chunk) if WORD.match(c)]
            if not letters:
                continue
            if want_italic and not all(letters):
                roman_names.append((text, chunk))
            elif not want_italic and any(letters):
                # mathtext puts the spaces between roman groups in the base
                # font, which is italic in the flipped mode; only the words
                # themselves carry the claim, so only they are checked
                italic_others.append((text, "".join(
                    c for f, c in zip(got, chunk) if f and WORD.match(c))))
    return roman_names, italic_others


def main(paths):
    bad = 0
    for path in paths:
        roman_names, italic_others = audit(path)
        name = os.path.basename(path)
        for line, chunk in roman_names:
            print(f"  {name}: {chunk!r} is roman, in {line.strip()!r}")
        for line, chunk in italic_others:
            print(f"  {name}: {chunk!r} is italic but is not a name, "
                  f"in {line.strip()!r}")
        bad += len(roman_names) + len(italic_others)
    print(f"check_figure_italics: {len(paths)} figures, {bad} problem(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    args = sys.argv[1:] or sorted(glob.glob(os.path.join(OUT, "Figure*.pdf")))
    sys.exit(main(args))
