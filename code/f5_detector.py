#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 5：偵測器的效能與這個任務的真實難度。"""
import csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SHORT = {"easy (far from machinery)": "easy",
         "hard: inside a cluster window": "in a\ncluster",
         "hard: composition-shuffled positives": "shuffled",
         "hard: length-matched short CDS": "length-\nmatched"}


def main():
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    base = json.load(open(os.path.join(OUT, "baselines.json")))
    paired = json.load(open(os.path.join(OUT, "paired_tests.json")))
    hard = json.load(open(os.path.join(OUT, "hard_negatives.json")))
    ladder = json.load(open(os.path.join(OUT, "ladder.json")))
    rule = json.load(open(os.path.join(OUT, "rule_baseline.json")))

    fig, axes = plt.subplots(2, 3, figsize=(F.WIDTH_IN, 4.8))
    (axA, axB, axC), (axD, axE, axF) = axes

    # A PR 曲線
    styles = {"ESM-2 + logistic": (F.NAVY, "-"), "3-mer counts": (F.AMBER, "--"),
              "composition": (F.GREEN, "-."), "nearest neighbour (10-mer)": (F.GREY, ":")}
    for m, cur in base["curves"].items():
        c, ls = styles.get(m, (F.FAINT, "-"))
        pr, rc = np.array(cur[0]), np.array(cur[1])
        axA.plot(rc, pr, color=c, ls=ls, lw=1.2,
                 label=f"{m.split(' (')[0]} {base['methods'][m]['cv_ap']:.2f}")
    axA.set_xlabel("recall"); axA.set_ylabel("precision")
    axA.set_xlim(0, 1); axA.set_ylim(0, 1.04)
    axA.set_xticks([0, 0.5, 1.0]); axA.set_yticks([0, 0.5, 1.0])
    # PR 曲線在左下角會降到 0.3 附近，圖例擺 lower left 會被 composition 與
    # nearest-neighbour 兩條穿過；移到軸下方，兩欄兩列
    axA.legend(frameon=False, fontsize=5.2, loc="upper center",
               bbox_to_anchor=(0.5, -0.30), ncol=2, handlelength=1.6,
               columnspacing=1.0, labelspacing=0.3)
    axA.set_title("grouped CV, easy negatives", fontsize=7, loc="left")
    F.bare(axA); F.panel(axA, "A", dx=-34)

    # B 配對差值（對最強古典基線）
    labs = ["ESM-2", "3-mer", "composition", "10-mer NN", "leader rules"]
    vals = [base["methods"][m]["cv_ap"] for m in
            ("ESM-2 + logistic", "3-mer counts", "composition",
             "nearest neighbour (10-mer)")] + [rule["ap"]]
    y = np.arange(len(labs))[::-1]
    axB.barh(y, vals, color=[F.NAVY, F.AMBER, F.GREEN, F.GREY, F.BROWN],
             edgecolor="white", linewidth=0.4, height=0.68)
    lo, hi = rep["bootstrap"]["ap_language_model_ci"]
    axB.plot([lo, hi], [y[0]] * 2, color=F.INK, lw=1.2)
    for yi, v in zip(y, vals):
        if v < 0.15:      # 0.02 的長條在 0–1 的軸上等於看不見，把數值寫出來
            axB.text(v + 0.02, yi, f"{v:.3f}", va="center", fontsize=5.6, color=F.INK)
    axB.set_yticks(y); axB.set_yticklabels(labs, fontsize=6.4)
    axB.set_xlim(0, 1.05); axB.set_xticks([0, 0.5, 1.0])
    axB.set_ylim(-0.75, len(labs) - 0.25)
    axB.tick_params(axis="x", pad=4)
    d3 = paired["language_model_vs_3mer"]
    axB.set_xlabel("average precision", fontsize=7.5)
    axB.set_title(f"paired difference vs 3-mer\n{d3['difference']:.2f} "
                  f"({d3['ci'][0]:.2f}–{d3['ci'][1]:.2f})", fontsize=6.2, loc="left")
    F.bare(axB); F.panel(axB, "B", dx=-44)

    # C easy vs hard 負例
    names_ = [k for k in SHORT if k in hard["negative_sets"]]
    x = np.arange(len(names_)); w = 0.38
    ap = [hard["negative_sets"][k]["ap"] for k in names_]
    chance = [hard["negative_sets"][k]["prevalence"] for k in names_]
    auroc = [hard["negative_sets"][k]["auroc"] for k in names_]
    axC.bar(x - w / 2, ap, w, color=F.NAVY, edgecolor="white", linewidth=0.4, label="AP")
    axC.bar(x + w / 2, auroc, w, color=F.TEAL, hatch="///", edgecolor="white",
            linewidth=0.4, label="AUROC")
    axC.scatter(x - w / 2, chance, marker="_", s=90, color=F.RED, zorder=4,
                label="AP of random guessing")
    # The length-only model was quoted in the legend and drawn nowhere, so a
    # reader had to take it on trust. It belongs on this panel and not on B:
    # panel B is cross-validated average precision on the training material,
    # while this is average precision on a named negative set, and the two are
    # not the same measurement. Here it is the same metric on the same set as
    # the AP bar it sits on.
    ap_len = [hard["negative_sets"][k]["ap_length_only"] for k in names_]
    axC.scatter(x - w / 2, ap_len, marker="v", s=22, color=F.BROWN, zorder=5,
                edgecolor="white", linewidth=0.3, label="AP, length alone")
    axC.set_xticks(x); axC.set_xticklabels([SHORT[k] for k in names_], fontsize=5.4)
    axC.set_ylim(0, 1.15); axC.set_yticks([0, 0.5, 1.0])
    # 這個圖例比格子寬，錨在 (0.5, 1.06) 的「upper center」等於把它的上緣釘在
    # 軸框裡面，左欄兩個標記就壓在 y 軸軸線上（紅色的 _ 被軸線穿過）。改成
    # 錨在軸框上緣的「lower center」，整塊移到軸框外，標題再往上讓位。
    axC.legend(frameon=False, fontsize=5.0, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, 1.01), handlelength=1.0, columnspacing=0.6,
               labelspacing=0.25)
    axC.set_title("negatives of different difficulty", fontsize=6.6, loc="left", pad=26)
    axC.set_xlabel("negative set", fontsize=6.4)
    F.bare(axC); F.panel(axC, "C", dx=-34, dy=14)

    # D 相似度分層
    st = hard["identity_strata"]
    keys = [k for k in ("<30%", "30-50%", "50-70%", ">=70%") if k in st]
    x = np.arange(len(keys))
    axD.bar(x, [st[k]["ap"] for k in keys], 0.6, color=F.NAVY, edgecolor="white",
            linewidth=0.4)
    for xi, k in zip(x, keys):
        ci = st[k].get("ap_ci")
        if ci and ci[0] is not None:
            axD.plot([xi, xi], ci, color=F.INK, lw=1.0)
        axD.text(xi, min(1.02, (ci[1] if ci and ci[1] else st[k]["ap"]) + 0.04),
                 f"n={st[k]['n_positives']}", ha="center", fontsize=5.4, color=F.INK)
    axD.set_xticks(x); axD.set_xticklabels(keys, fontsize=6)
    axD.set_ylim(0, 1.1); axD.set_yticks([0, 0.5, 1.0])
    axD.set_ylabel("average precision", fontsize=7)
    axD.set_xlabel("10-mer identity to the nearest\nother training peptide", fontsize=6.4)
    F.bare(axD); F.panel(axD, "D", dx=-38)

    # E 學習曲線
    lc = hard["learning_curve"]
    axE.plot([d["n_positives"] for d in lc], [d["ap"] for d in lc], marker="o",
             color=F.NAVY, lw=1.2)
    axE.set_xlabel("training peptides", fontsize=7)
    axE.set_ylabel("average precision", fontsize=7)
    axE.set_ylim(0, 1.05); axE.set_yticks([0, 0.5, 1.0])
    axE.set_title("still rising at 149", fontsize=6.6, loc="left")
    F.bare(axE); F.panel(axE, "E", dx=-38)

    # F 模型階梯
    SIZE = {"8M": 8, "35M": 35, "150M": 150}
    items = sorted((SIZE[k], k, v) for k, v in ladder.items() if k in SIZE)
    axF.plot([i[0] for i in items], [i[2]["cv_ap"] for i in items], marker="o",
             color=F.NAVY, lw=1.2)
    for s_, k, v in items:
        axF.text(s_, v["cv_ap"] + 0.04, k, fontsize=5.8, ha="center")
    axF.set_xscale("log"); axF.set_xlim(5, 400)
    axF.set_xticks([8, 35, 150])
    axF.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    axF.set_xlabel("ESM-2 parameters (M)", fontsize=7)
    axF.set_ylabel("average precision", fontsize=7)
    axF.set_ylim(0, 1.12); axF.set_yticks([0, 0.5, 1.0])
    F.bare(axF); F.panel(axF, "F", dx=-38)

    fig.tight_layout(h_pad=2.4, w_pad=2.2)
    F.save(fig, "Figure5_detector")


if __name__ == "__main__":
    main()
