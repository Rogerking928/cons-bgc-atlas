#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure S9：跟已發表的工具（DeepBGC）在同一段 DNA 上對照。"""
import csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")


def main():
    deep = list(csv.DictReader(open(os.path.join(OUT, "deepbgc_results.csv"))))
    loci = [r for r in deep if not r["region"].startswith("control")]
    ctrl = [r for r in deep if r["region"].startswith("control")]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(F.WIDTH_IN, 2.8),
                                   gridspec_kw={"width_ratios": [1.0, 1.25]})
    cats = ["≥1 BGC call", "called RiPP"]
    v_loci = [100 * sum(1 for r in loci if int(r["bgcs"]) > 0) / len(loci),
              100 * sum(1 for r in loci if int(r["ripp"]) > 0) / len(loci)]
    v_ctrl = [100 * sum(1 for r in ctrl if int(r["bgcs"]) > 0) / len(ctrl),
              100 * sum(1 for r in ctrl if int(r["ripp"]) > 0) / len(ctrl)]
    x = np.arange(len(cats)); w = 0.38
    axA.bar(x - w / 2, v_loci, w, color=F.NAVY, edgecolor="white", linewidth=0.4,
            label=f"our cluster loci (n={len(loci)})")
    axA.bar(x + w / 2, v_ctrl, w, color=F.FAINT, edgecolor="white", linewidth=0.4,
            label=f"control windows (n={len(ctrl)})")
    axA.set_xticks(x); axA.set_xticklabels(cats, fontsize=6.6)
    axA.set_ylabel("windows (%)", fontsize=7.2)
    axA.set_ylim(0, 115); axA.set_yticks([0, 50, 100])
    axA.legend(frameon=False, fontsize=5.6, loc="upper right", handlelength=1.1)
    F.bare(axA); F.panel(axA, "A", dx=-40)

    y = np.arange(len(loci))[::-1]
    for yi, r in zip(y, sorted(loci, key=lambda r: r["region"])):
        called = int(r["bgcs"]) > 0
        ripp = int(r["ripp"]) > 0
        axB.scatter(0, yi, s=26, marker="o" if called else "x",
                    color=F.NAVY if called else F.GREY, zorder=3)
        axB.scatter(1, yi, s=26, marker="o" if ripp else "x",
                    color=F.AMBER if ripp else F.GREY, zorder=3)
        axB.text(1.35, yi, f"{r['region']}  {names.fam_short(r['region'])}",
                 fontsize=5.6, va="center", color=F.INK)
        if r["best_score"]:
            axB.text(0.5, yi + 0.26, f"score {float(r['best_score']):.2f}",
                     fontsize=4.8, ha="center", color=F.GREY)
    axB.set_xlim(-0.4, 3.4); axB.set_ylim(-0.6, len(loci) - 0.4)
    axB.set_xticks([0, 1]); axB.set_xticklabels(["BGC called", "classed RiPP"],
                                                fontsize=6.4)
    axB.set_yticks([])
    F.bare(axB, keep=("bottom",)); F.panel(axB, "B", dx=-14)

    fig.tight_layout(w_pad=3.0)
    F.save(fig, "FigureS9_deepbgc")


if __name__ == "__main__":
    main()
