#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure S1：面板描述（物種與分離部位）。攜帶率那一格已經移到主圖 Figure 1A，
釋出年份那一格不支撐任何主張，直接拿掉。"""
import collections, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F
import names

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out")


def load():
    md = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "metadata.csv")))}
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    return md, mat, cl


def main():
    md, mat, cl = load()
    sp_of = {a: r["species"] for a, r in mat.items()}
    order = [s for s in F.SPECIES_ORDER if s != "other"]
    counts = collections.Counter(sp_of.values())
    rare = sorted(s for s in counts if s not in order)
    kinds = collections.defaultdict(lambda: collections.defaultdict(set))
    for c in cl:
        kinds[sp_of[c["assembly"]]][c["kind"]].add(c["assembly"])

    fig, (axA, axC) = plt.subplots(1, 2, figsize=(F.WIDTH_IN, 2.9))

    # A 每個物種幾株
    labs = order + ["other (%d spp.)" % len(rare)]
    vals = [counts[s] for s in order] + [sum(counts[s] for s in rare)]
    y = np.arange(len(labs))[::-1]
    for yi, s, v in zip(y, order + ["other"], vals):
        axA.barh(yi, v, color=F.sp_color(s), hatch=F.sp_hatch(s),
                 edgecolor="white", linewidth=0.4, height=0.72)
        axA.text(v + 6, yi, str(v), va="center", fontsize=8, color=F.INK)
    axA.set_yticks(y)
    # style="italic" 會連 "other (4 spp.)" 一起斜；交給 figstyle 的斜體掃描
    axA.set_yticklabels([("S. " + l) if not l.startswith("other") else l for l in labs],
                        fontsize=8)
    axA.set_xlim(0, max(vals) * 1.18)
    axA.set_xticks([0, 100, 200, 300])
    axA.set_xlabel("complete genomes")
    F.bare(axA); F.panel(axA, "A", dx=-96)

    # C 分離部位
    g = collections.Counter(names.source_group(md[a]["isolation_source"]) for a in mat)
    keys = ["blood", "device / implant", "skin / nose", "urine",
            "other clinical", "other / environmental", "not stated"]
    v = [g[k] for k in keys]
    yy = np.arange(len(keys))[::-1]
    axC.barh(yy, v, color=[F.NAVY, F.RED, F.AMBER, F.TEAL, F.GREEN, F.GREY, F.FAINT],
             edgecolor="white", linewidth=0.4, height=0.72)
    for yi, vv in zip(yy, v):
        axC.text(vv + 3, yi, str(vv), va="center", fontsize=8, color=F.INK)
    axC.set_yticks(yy); axC.set_yticklabels(keys, fontsize=8)
    axC.set_xlim(0, max(v) * 1.2)
    axC.set_xticks([0, 50, 100, 150, 200])
    axC.set_xlabel("genomes")
    F.bare(axC); F.panel(axC, "B", dx=-104)

    fig.tight_layout(h_pad=1.6, w_pad=2.4)
    F.save(fig, "FigureS1_panel")


if __name__ == "__main__":
    main()
