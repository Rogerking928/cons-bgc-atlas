#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure S7：LoRA 微調 vs 凍結（同一批序列、同一組 fold 的配對比較）。
這是一個對照實驗，不是主張，所以放補充。"""
import collections, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
GROUPS = [("profile-detected precursor", F.NAVY),
          ("candidate, q < 0.05", F.RED),
          ("candidate, q < 0.20", F.AMBER),
          ("negative", F.GREY)]


def main():
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    ft = json.load(open(os.path.join(OUT, "finetune_ap.json")))
    sal = json.load(open(os.path.join(OUT, "finetune_summary.json")))
    rows = list(csv.DictReader(open(os.path.join(OUT, "saliency.csv"))))
    paired = json.load(open(os.path.join(OUT, "paired_tests.json")))

    # 歸因曲線（四組都在 motif 上有峰、單看分不開）縮成 inset，把版面讓給
    # 真正分得開的統計量（2026-09-20 審閱建議）
    fig, axA = plt.subplots(figsize=(F.WIDTH_IN * 0.5, 2.6))

    # A 凍結 vs 微調（**同一批序列、同一組 fold** 的配對比較）
    # 舊版拿全部 18,695 條負例訓練出來的 0.894 去比微調的 0.949，兩者的負例集
    # 根本不同，那不是配對比較。現在兩邊都用 s33 的子樣本（149 正 ＋ 3,000 負）
    # 與同一組 GroupKFold，逐折相減。（2026-09-20 審閱意見）
    pf = paired["finetuned_vs_frozen_same_folds"]
    fr, ftf = pf["frozen_folds"], pf["finetuned_folds"]
    axA.bar([0, 1], [pf["frozen_mean"], pf["finetuned_mean"]], 0.5,
            color=[F.GREY, F.NAVY], edgecolor="white", linewidth=0.4)
    for x, v in zip((0, 1), (fr, ftf)):
        axA.scatter([x] * len(v), v, s=11, color=F.INK, zorder=3)
    for x, v in zip((0, 1), (pf["frozen_mean"], pf["finetuned_mean"])):
        axA.text(x, 1.03, f"{v:.3f}", ha="center", fontsize=7, color=F.INK)
    for a_, b_ in zip(fr, ftf):
        axA.plot([0, 1], [a_, b_], color=F.FAINT, lw=0.7, zorder=2)
    axA.set_xticks([0, 1])
    axA.set_xticklabels(["frozen 35M\n+ logistic", "LoRA\nfine-tuned"], fontsize=7)
    axA.set_xlim(-0.6, 1.6)
    axA.set_ylim(0, 1.18); axA.set_yticks([0, 0.5, 1.0])
    axA.set_ylabel("average precision", fontsize=7.5)
    axA.set_title(f"paired difference {pf['paired_difference']:+.3f} "
                  f"({pf['ci_t'][0]:+.3f} to {pf['ci_t'][1]:+.3f})",
                  fontsize=6.4, loc="left")
    F.bare(axA); F.panel(axA, "A", dx=-40)

    F.save(fig, "FigureS7_finetune")


if __name__ == "__main__":
    main()
