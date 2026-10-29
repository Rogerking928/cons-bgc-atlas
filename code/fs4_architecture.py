#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 6：合成酶的功能域結構、專一性碼，以及島與基因體的 GC 差。"""
import collections, csv, gzip, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Patch

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
DOMCOL = {"A": F.NAVY, "C": F.RED, "PCP/ACP": F.AMBER, "TE": F.GREEN,
          "MT": F.PURPLE, "KR": F.TEAL, "ER/NAD": F.BROWN, "A_C": "#7f8c99",
          "YcaO": F.PLUM, "LanB_N": F.NAVY, "LanB_C": "#56708a",
          "LanC": F.GREEN, "LanM_cyc": F.RED}
# 家族重新編號過（2026-09-20）：F005 是 hybrid、F006 才是 lugdunin
SHOW = [("F006", "lugdunin"), ("F001", "aureusimine (epidermidis)"),
        ("F003", "aryl-capped NRPS (lugdunensis)"), ("F005", "NRPS x YcaO hybrid")]
# 蛋白名在 y 軸上截斷會讀不懂（"bifunctional Gfo/Idh/MocA family o"）：
# 常見的先換成短名，其餘的在空白處折行，不要硬切。
SHORTNAME = [
    ("lugdunin non-ribosomal peptide synthetase LugA", "LugA"),
    ("lugdunin non-ribosomal peptide synthetase LugB", "LugB"),
    ("lugdunin non-ribosomal peptide synthetase LugC", "LugC"),
    ("lugdunin non-ribosomal peptide synthetase LugD", "LugD"),
    ("aureusimine non-ribosomal peptide synthetase AusA", "AusA"),
    ("(2,3-dihydroxybenzoyl)adenylate synthase", "DHB-AMP synthase"),
    ("bifunctional Gfo/Idh/MocA family oxidoreductase", "Gfo/Idh/MocA oxidoreductase"),
    ("condensation domain-containing protein", "condensation protein"),
    ("non-ribosomal peptide synthetase", "NRPS"),
    ("YcaO-like family protein", "YcaO protein"),
    ("thioesterase II family protein", "thioesterase II"),
]


def prettify(name, width=30):
    for long, short in SHORTNAME:
        if long.lower() in (name or "").lower():
            return short
    name = (name or "").split(" domain-containing")[0]
    if len(name) <= width:
        return name
    cut = name.rfind(" ", 0, width)
    return name[:cut] + "\n" + name[cut + 1:][:width] if cut > 0 else name[:width]


def main():
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cds = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"):
        cds[r["cluster"]].append(r)
    dom = collections.defaultdict(list)
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            a, prot, pf, lab, af, at, ev, sc, pl = line.rstrip("\n").split("\t")
            dom[(a, prot)].append((int(af), int(at), lab, int(pl)))
    prof = json.load(open(os.path.join(OUT, "gc_profiles.json")))

    fig = plt.figure(figsize=(F.WIDTH_IN, 6.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.6, 1.15, 0.95], hspace=0.85, wspace=0.42)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, :])
    axD = fig.add_subplot(gs[2, :])

    # A 功能域結構
    rows, labels = [], []
    for fam, title in SHOW:
        cs = [c for c, f in fam_of.items() if f == fam]
        rep = cs[0]
        prots = [(g, dom.get((g["assembly"], g["protein"]), []))
                 for g in cds[rep] if int(g["aa"] or 0) >= 400 and g["domains"]]
        prots = [p for p in prots if p[1]]
        prots.sort(key=lambda p: -int(p[0]["aa"]))
        for g, ds in prots[:4]:
            rows.append((g, sorted(ds)))
            labels.append(f"{fam}  {prettify(g['product'])}")
    maxlen = max(int(g["aa"]) for g, _ in rows)
    y = np.arange(len(rows))[::-1]
    for yi, (g, ds) in zip(y, rows):
        L = int(g["aa"])
        axA.add_patch(Rectangle((0, yi - 0.17), L, 0.34, facecolor=F.FAINT,
                                edgecolor="none"))
        for af, at, lab, pl in ds:
            c = DOMCOL.get(lab)
            if c is None:
                continue
            axA.add_patch(Rectangle((af, yi - 0.3), at - af, 0.6, facecolor=c,
                                    edgecolor="white", linewidth=0.4))
        axA.text(L + 40, yi, f"{L:,} aa", va="center", fontsize=6, color=F.GREY)
    axA.set_xlim(-60, maxlen * 1.16)
    axA.set_xticks([0, 500, 1000, 1500, 2000, 2500])
    axA.set_ylim(-0.8, len(rows) - 0.2)
    axA.set_yticks(y); axA.set_yticklabels(labels, fontsize=5.6)
    axA.set_xlabel("residue")
    handles = [Patch(facecolor=DOMCOL[k], label=k) for k in
               ("A", "C", "PCP/ACP", "TE", "MT", "KR", "ER/NAD", "YcaO")]
    axA.legend(handles=handles, frameon=False, fontsize=6.4, ncol=8,
               loc="upper center", bbox_to_anchor=(0.5, 1.13), handlelength=1.1,
               columnspacing=1.0)
    F.bare(axA); F.panel(axA, "A", dx=-150, dy=14)

    # B 專一性碼
    if os.path.exists(os.path.join(OUT, "adomain_codes.csv")):
        ad = list(csv.DictReader(open(os.path.join(OUT, "adomain_codes.csv"))))
        byfam = collections.defaultdict(collections.Counter)
        SHORTP = [("lugdunin non-ribosomal peptide synthetase Lug", "Lug"),
                  ("aureusimine non-ribosomal peptide synthetase AusA", "AusA"),
                  ("(2,3-dihydroxybenzoyl)adenylate synthase", "DHB-AMP"),
                  ("non-ribosomal peptide synthetase", "NRPS"),
                  ("condensation domain-containing protein", "C-protein"),
                  ("AMP-binding protein", "AMP-binding")]

        def pshort(p):
            for long, short in SHORTP:
                if long.lower() in (p or "").lower():
                    return short + (p.strip()[-1] if short == "Lug" else "")
            return (p or "")[:10]

        for r in ad:
            byfam[r["family"]][(pshort(r["product"]) + " " + r["domain"],
                                r["signature"], r["closest"], r["identity"])] += 1
        fams = [f for f, _ in SHOW if f in byfam]
        rowlab, cells = [], []
        for fam in fams:
            for (pos, code, closest, ident), n in byfam[fam].most_common(3):
                rowlab.append(f"{fam} {pos}  n={n}  ~{closest} {ident}/10")
                cells.append(code)
        y = np.arange(len(cells))[::-1]
        # 不畫格子（格子是 patch，figstyle.check() 會把格子裡的字判成壓字）：
        # 一條細底線就夠分行了
        for yi, code in zip(y, cells):
            axB.plot([-0.1, 10.1], [yi - 0.4, yi - 0.4], color=F.FAINT, lw=0.6,
                     zorder=1)
            for j, ch in enumerate(code):
                axB.text(j + 0.46, yi, ch, ha="center", va="center",
                         fontsize=6.6, family="monospace", color=F.INK, zorder=3)
        axB.set_xlim(-0.2, 10.4); axB.set_ylim(-0.8, len(cells) - 0.2)
        axB.set_yticks(y); axB.set_yticklabels(rowlab, fontsize=5.2)
        axB.set_xticks([]); axB.set_xlabel("Stachelhaus specificity code (10 positions)")
        F.bare(axB, keep=()); F.panel(axB, "B", dx=-150)

    # C 島與基因體的 GC 差
    # C 代表座的 GC 視窗曲線
    pick = None
    for cid, d in prof.items():
        if d["family"] == "F006":
            pick = (cid, d); break
    if pick is None:
        pick = next(iter(prof.items()))
    cid, d = pick
    # gc 是整條 replicon 的視窗序列，start/end 是基因簇本身的座標
    xs = np.arange(len(d["gc"])) * d["window"] / 1000
    axD.plot(xs, d["gc"], color=F.INK, lw=0.9)
    axD.axhline(d["genome_gc"], color=F.RED, lw=0.9, ls="--",
                label=f"genome mean {d['genome_gc']:.1f}%")
    cl = {r["cluster"]: r for r in csv.DictReader(
        open(os.path.join(OUT, "clusters.tsv")), delimiter="\t")}
    c = cl[cid]
    axD.axvspan(int(c["start"]) / 1000, int(c["end"]) / 1000, color=F.AMBER,
                alpha=0.25, label=f"{d['family']} cluster")
    # 只看島的前後 40 kb：整條 replicon 的曲線擠在一起看不出島在哪
    lo = (int(c["start"]) - 40000) / 1000
    hi = (int(c["end"]) + 40000) / 1000
    axD.set_xlim(max(0, lo), min(xs[-1], hi))
    axD.set_xlabel(f"position on {d['replicon']} (kb)", fontsize=7.5)
    axD.set_ylabel("GC in 2 kb windows (%)", fontsize=7.5)
    # 角落那一格：y 的最低刻度會和 x 的 0 撞在一起，砍掉最低那一格
    axD.yaxis.set_major_locator(plt.MaxNLocator(nbins=4, prune="lower"))
    axD.legend(frameon=False, fontsize=5.8, loc="upper center",
               bbox_to_anchor=(0.5, 1.22), ncol=2, handlelength=1.3)
    F.bare(axD); F.panel(axD, "C", dx=-52)

    F.save(fig, "FigureS4_architecture")


if __name__ == "__main__":
    main()
