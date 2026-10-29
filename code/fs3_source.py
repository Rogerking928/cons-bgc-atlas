#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure S3：分離部位與帶菌（描述用）。

這張留在補充：分離部位跟送件計畫綁在一起，正文不拿它做任何推論，
放在主圖只會引誘讀者誤讀。
"""
import csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
GROUPS = ["blood", "device / implant", "skin / nose", "other clinical"]


def main():
    md = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "metadata.csv")))}
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(F.WIDTH_IN, 2.8))

    for ax, sp, fam, letter in ((axA, "epidermidis", "F001", "A"),
                                (axB, "lugdunensis", "F006", "B")):
        pool = [a for a, r in mat.items() if r["species"] == sp]
        vals, ns = [], []
        for g in GROUPS:
            sub = [a for a in pool if names.source_group(md[a]["isolation_source"]) == g]
            ns.append(len(sub))
            vals.append(100 * sum(1 for a in sub if mat[a].get(fam)) / len(sub) if sub else 0)
        y = np.arange(len(GROUPS))[::-1]
        ax.barh(y, vals, color=[F.NAVY, F.RED, F.AMBER, F.GREEN], edgecolor="white",
                linewidth=0.4, height=0.7)
        for yi, v, n in zip(y, vals, ns):
            ax.text(v + 3, yi, f"{v:.0f}%  (n={n})", va="center", fontsize=6.4,
                    color=F.INK)
        ax.set_yticks(y); ax.set_yticklabels(GROUPS, fontsize=7)
        ax.set_xlim(0, 128); ax.set_xticks([0, 50, 100])
        ax.set_xlabel(f"S. {sp} genomes carrying {fam} (%)", fontsize=7.5)
        F.bare(ax); F.panel(ax, letter, dx=-86)

    fig.tight_layout(w_pad=3.0)
    F.save(fig, "FigureS3_source")


if __name__ == "__main__":
    main()
