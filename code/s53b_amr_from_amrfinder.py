#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Turn the AMRFinderPlus output of s53 into the two files the pipeline reads.

寫出的檔名與 s8／s8b 完全一樣（out/amr_per_genome.csv、out/amr_located.csv），
所以下游 13 支腳本不用改一行。原本的 containment 結果改名保存為
out/amr_per_genome_containment.csv、out/amr_located_containment.csv，
兩套呼叫的差異本身要寫進稿件（審稿人 2 問的就是這件事）。

為了讓範圍與原稿一致：
- Element type = AMR 全收（未給 --organism，所以不會有點突變呼叫）。
- 另收 STRESS/BIOCIDE 裡的四級銨（qacA/qacC 那組）——原稿的 CLASSES 本來就含
  QUATERNARY_AMMONIUM。VIRULENCE 不收。
- Method 一併記下來（EXACTX／ALLELEX／BLASTX／HMMX／PARTIALX…），PARTIALX 這類
  部分比對留著但標記，敏感度分析可以整組拿掉。

AMRFinderPlus 已經套用每個基因策劃好的門檻，所以「有回報」就等於「有呼叫」：
containment 欄位寫 1.0，下游的 `>= CALL (0.5)` 判斷因此剛好等於它的呼叫。
真正的比對品質改記在 identity 與 coverage 兩欄。
"""
import collections, csv, glob, gzip, os, shutil, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s8_amr_screen import CLASSES, FAM_RE  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SRC = os.path.join(OUT, "amrfinder")

# v4 改了欄名（Gene symbol → Element symbol），兩種都接。
# 實測 4.2.7 的型別欄是 Type／Subtype，不是 Element type／Element subtype。
SYM = ("Element symbol", "Gene symbol")
NAME = ("Element name", "Sequence name")
ETYPE = ("Type", "Element type")
ESUB = ("Subtype", "Element subtype")


def pick(row, names, default=""):
    for n in names:
        if n in row and row[n] not in (None, ""):
            return row[n]
    return default


def wanted(row):
    """原稿的範圍，一字不差：s8 的 FAM_RE（基因符號）且 CLASSES（藥物類別）。

    2026-10-29 改。第一版只按 Type 篩（AMR 全收＋qac），結果把調節子
    （blaI／blaR1、mecR1／mecI、qacR）、fosB、bleO、dfrS1、泛稱的 `bla` 也收進
    「後天抗藥基因」，進入距離分析的染色體簇從 396 變成 500 —— 那是範圍變了，
    不是呼叫方法變了。這裡沿用原稿同一套白名單，讓兩套呼叫只差在方法。
    AMRFinderPlus 的 Class 用空格（FUSIDIC ACID），s8 的 CLASSES 用底線。
    """
    etype = pick(row, ETYPE).upper()
    if etype not in ("AMR", "STRESS"):
        return False
    if not FAM_RE.match(pick(row, SYM)):
        return False
    cls = (row.get("Class") or "").upper().replace(" ", "_")
    return bool(set(cls.split("/")) & CLASSES)


def replicon_map():
    """(assembly, replicon) -> (replicon_type, plasmid)，取自 out/cds.tsv.gz。"""
    rep = {}
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rep[(f[0], f[1])] = (f[2], f[4])
    return rep


def main():
    files = sorted(glob.glob(os.path.join(SRC, "*.tsv")))
    if not files:
        sys.exit("no AMRFinderPlus output in out/amrfinder/ - run s53 first")
    panel = sorted({r["assembly"] for r in
                    csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv")))})
    rep = replicon_map()

    hits = collections.defaultdict(list)   # assembly -> [(sym, contig, start, end, ident, cov, method)]
    skipped_scope = collections.Counter()
    for path in files:
        acc = os.path.basename(path)[:-4]
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sym = pick(row, SYM)
                if not sym:
                    continue
                if not wanted(row):
                    skipped_scope[pick(row, ETYPE) or "?"] += 1
                    continue
                contig = row.get("Contig id") or ""
                try:
                    start, end = int(row["Start"]), int(row["Stop"])
                except (KeyError, TypeError, ValueError):
                    start, end = "", ""
                hits[acc].append((sym, contig, start, end,
                                  row.get("% Identity to reference", ""),
                                  row.get("% Coverage of reference", ""),
                                  row.get("Method", ""),
                                  row.get("Protein id") or row.get("Protein identifier", "")))

    syms = sorted({h[0] for v in hits.values() for h in v})
    print(f"{len(files)} genomes parsed, {len(syms)} determinants, "
          f"{sum(len(v) for v in hits.values())} calls", file=sys.stderr)
    print(f"out of scope: {dict(skipped_scope)}", file=sys.stderr)

    # 先把舊的 containment 結果留起來，再覆蓋正規檔名
    for base in ("amr_per_genome", "amr_located"):
        src = os.path.join(OUT, base + ".csv")
        keep = os.path.join(OUT, base + "_containment.csv")
        if os.path.exists(src) and not os.path.exists(keep):
            shutil.copy2(src, keep)
            print(f"kept {os.path.basename(keep)}", file=sys.stderr)

    with open(os.path.join(OUT, "amr_per_genome.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["assembly"] + syms)
        w.writeheader()
        for acc in panel:
            row = {"assembly": acc}
            present = {h[0] for h in hits.get(acc, [])}
            for s in syms:
                row[s] = 1.0 if s in present else 0.0
            w.writerow(row)

    with open(os.path.join(OUT, "amr_located.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["assembly", "determinant", "protein", "containment",
                    "replicon", "replicon_type", "plasmid", "start", "end",
                    "identity", "coverage", "method"])
        n = 0
        for acc in panel:
            for sym, contig, start, end, ident, cov, method, pid in sorted(hits.get(acc, [])):
                rtype, plasmid = rep.get((acc, contig), ("", ""))
                w.writerow([acc, sym, pid, 1.0, contig, rtype, plasmid,
                            start, end, ident, cov, method]); n += 1
    print(f"wrote out/amr_per_genome.csv ({len(panel)} genomes) and "
          f"out/amr_located.csv ({n} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
