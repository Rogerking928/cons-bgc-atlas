#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure S6：分數分布與保留集（主圖 Figure 5、6 的支持材料）。"""
import csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")


def main():
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    base = json.load(open(os.path.join(OUT, "baselines.json")))
    scan = list(csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv"))))
    hold = list(csv.DictReader(open(os.path.join(OUT, "precursor_holdout.csv"))))
    thr = rep["threshold_999th_percentile_of_negatives"]

    seen = {}
    for r in scan:
        seen.setdefault(r["seq"], r)
    hs = {h["seq"] for h in hold}
    pos = [float(r["score"]) for q, r in seen.items() if r["known_domain"]]
    neg = np.array([float(r["score"]) for q, r in seen.items()
                    if not r["known_domain"] and q not in hs])
    ho = [float(r["score"]) for q, r in seen.items() if q in hs]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(F.WIDTH_IN, 2.8))
    bins = np.linspace(0, 1, 41)
    axA.hist(neg, bins=bins, color=F.FAINT, log=True, bottom=0.8,
             label=f"other distinct CDS (n={len(neg):,})")
    axA.hist(ho, bins=bins, color=F.AMBER, log=True, bottom=0.8, alpha=0.9,
             label=f"RefSeq-named, no profile (n={len(ho)})")
    axA.hist(pos, bins=bins, color=F.NAVY, log=True, bottom=0.8, alpha=0.9,
             label=f"profile precursor (n={len(pos)})")
    axA.axvline(thr, color=F.RED, lw=1.0)
    axA.text(thr - 0.04, 1.5e3, "q<0.05", fontsize=6, color=F.RED, ha="right")
    axA.set_ylim(0.8, 3e5); axA.set_yticks([1, 100, 10000])
    axA.set_xticks([0, 0.5, 1.0])
    axA.set_xlabel("detector score"); axA.set_ylabel("distinct sequences")
    axA.legend(frameon=False, fontsize=5.4, loc="upper center")
    F.bare(axA); F.panel(axA, "A", dx=-40)

    ms = ["ESM-2 + logistic", "3-mer counts", "composition", "nearest neighbour (10-mer)"]
    labs = ["ESM-2", "3-mer", "composition", "10-mer NN"]
    v = [base["methods"][m]["holdout_ap"] or 0 for m in ms]
    y = np.arange(len(labs))[::-1]
    axB.barh(y, v, color=[F.NAVY, F.AMBER, F.GREEN, F.GREY], edgecolor="white",
             linewidth=0.4, height=0.66)
    for yi, vv in zip(y, v):
        axB.text(vv + 0.003, yi, f"{vv:.3f}", va="center", fontsize=6, color=F.INK)
    axB.set_yticks(y); axB.set_yticklabels(labs, fontsize=6.6)
    axB.set_xlim(0, 0.11)
    axB.set_xlabel("average precision on the held-out\nRefSeq-named set", fontsize=7)
    axB.set_title(f"recall at the calling threshold: "
                  f"{100 * rep['holdout_annotation']['recall_at_threshold']:.0f}%",
                  fontsize=6.6, loc="left")
    F.bare(axB); F.panel(axB, "B", dx=-48)

    fig.tight_layout(w_pad=3.0)
    F.save(fig, "FigureS6_scores")


if __name__ == "__main__":
    main()
