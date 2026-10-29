#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 4：這套檢驗（PCR 方案）在設計集與外部資料上各值多少。

原本分成兩張（設計 vs 外部），講的是同一件事，讀者要來回翻；併成一張，
並且把預測值同時用兩種盛行率算一遍（2026-09-20 審閱意見）。"""
import collections, csv, json, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")


def wilson(k, n):
    if n == 0:
        return (np.nan, np.nan)
    p, z = k / n, 1.96
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    ext = json.load(open(os.path.join(OUT, "external_pcr.json")))
    rows = list(csv.DictReader(open(os.path.join(OUT, "external_pcr.csv"))))
    deep = list(csv.DictReader(open(os.path.join(OUT, "deepbgc_results.csv"))))
    mat = list(csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv"))))
    best = {}
    for r in csv.DictReader(open(os.path.join(OUT, "primers.csv"))):
        k = r["target"]
        sc = (float(r["sensitivity"]), float(r["specificity"]), -int(r["product_bp"]))
        if k not in best or sc > best[k][0]:
            best[k] = (sc, r)
    panel = {k: v[1] for k, v in best.items()}
    targets = sorted(ext["per_target"])

    fig = plt.figure(figsize=(F.WIDTH_IN, 6.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.0, 1.0], hspace=0.95, wspace=0.42)
    axP = fig.add_subplot(gs[0, :])
    axA = fig.add_subplot(gs[1, 0])
    axB = fig.add_subplot(gs[1, 1])
    axC = fig.add_subplot(gs[2, :])

    # P 預測值：用設計集的盛行率，以及用外部 draft 的盛行率各算一次
    mat_rows = list(csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv"))))
    rngp = np.random.default_rng(20260920)
    yp = np.arange(len(targets))[::-1]
    for yi, t in zip(yp, targets):
        fam = t.split("-")[0]
        r = panel[t]
        sp = collections.Counter(x["species"] for x in mat_rows if x.get(fam))
        main_sp = sp.most_common(1)[0][0] if sp else ""
        tot = sum(1 for x in mat_rows if x["species"] == main_sp)
        car = sum(1 for x in mat_rows if x["species"] == main_sp and x.get(fam))
        d = [x for x in rows if x["group"] == "draft CoNS" and x["species"] == main_sp]
        car_d = sum(int(x[f"{t}_carrier"]) for x in d)
        se_a = int(r["amplified_carriers"]); se_b = int(r["carriers_tested"]) - se_a
        sp_b = int(r["amplified_non_carriers"])
        sp_a = int(r["non_carriers_tested"]) - sp_b
        S = rngp.beta(se_a + 0.5, se_b + 0.5, 4000)
        C = rngp.beta(sp_a + 0.5, sp_b + 0.5, 4000)
        vals = []
        for (c_, n_), col, mark in (((car, tot), F.GREY, "o"),
                                    ((car_d, len(d)), F.NAVY, "D")):
            if n_ == 0:
                continue
            PV = rngp.beta(c_ + 0.5, n_ - c_ + 0.5, 4000)
            ppv = S * PV / np.maximum(S * PV + (1 - C) * (1 - PV), 1e-12)
            lo, hi = np.percentile(ppv, [2.5, 97.5])
            vals.append((100 * np.median(ppv), 100 * lo, 100 * hi, col, mark,
                         100 * c_ / n_))
        # 門檻原本設 5 個百分點，結果 F002（100% vs 97.5%）在 A 被寫成 unchanged
        # 而 D 畫成兩個點，同一張圖自相矛盾。只要數字不同就畫兩個點。
        moved = len(vals) == 2 and round(vals[0][5]) != round(vals[1][5])
        if moved:
            axP.plot([vals[0][0], vals[1][0]], [yi, yi], color=F.FAINT, lw=1.0, zorder=1)
        # 兩個盛行率幾乎一樣的目標，兩個點會完全疊住，看起來像只算了一次：
        # 那幾列只畫一個點並註明 unchanged（2026-09-20 審閱意見）
        show = vals if moved else vals[:1]
        for med, lo, hi, col, mark, prev in show:
            axP.plot([lo, hi], [yi] * 2, color=col, lw=1.0, zorder=2)
            axP.scatter(med, yi, s=20, color=col, marker=mark, zorder=3)
        if len(vals) == 2:
            axP.text(103, yi,
                     (f"carriage {vals[0][5]:.0f}% \u2192 {vals[1][5]:.0f}%" if moved
                      else f"carriage {vals[0][5]:.0f}%, unchanged"),
                     fontsize=5.2, va="center", color=F.GREY)
    axP.set_yticks(yp); axP.set_yticklabels(targets, fontsize=6.2)
    axP.set_xlim(0, 132); axP.set_xticks([0, 50, 100])
    axP.set_xlabel("positive predictive value (%, 95% CI)", fontsize=7.2)
    # 圖例用代理標記：第一列有可能只畫了一個點，靠 label 會漏掉另一個
    hp = [plt.Line2D([], [], ls="", marker="o", color=F.GREY,
                     label="at design-panel carriage", markersize=5),
          plt.Line2D([], [], ls="", marker="D", color=F.NAVY,
                     label="at held-out draft carriage", markersize=5)]
    axP.legend(handles=hp, frameon=False, fontsize=5.6, loc="lower left",
               handlelength=1.2)
    axP.set_title("what a positive result would be worth, at two carriage estimates",
                  fontsize=6.6, loc="left")
    F.bare(axP); F.panel(axP, "A", dx=-52)

    # A 設計集 vs draft 的敏感度
    y = np.arange(len(targets))[::-1]
    for yi, t in zip(y, targets):
        p = panel[t]
        se_p = float(p["sensitivity"])
        lo, hi = wilson(int(p["amplified_carriers"]), int(p["carriers_tested"]))
        axA.plot([100 * lo, 100 * hi], [yi + 0.16] * 2, color=F.GREY, lw=1.0, zorder=2)
        axA.scatter(100 * se_p, yi + 0.16, s=18, color=F.GREY, marker="o", zorder=3,
                    label="design panel (554 complete)" if yi == y[0] else "")
        e = ext["per_target"][t]
        if e["sensitivity"] is None:
            continue
        lo, hi = e["sensitivity_ci"]
        axA.plot([100 * lo, 100 * hi], [yi - 0.16] * 2, color=F.NAVY, lw=1.0, zorder=2)
        axA.scatter(100 * e["sensitivity"], yi - 0.16, s=18, color=F.NAVY, marker="D",
                    zorder=3, label="held-out drafts" if yi == y[0] else "")
        # 標記畫在 100 上，101 就正好壓在上面；往右挪開、軸右緣跟著放寬
        axA.text(107, yi, f"{e['draft_carriers']}", fontsize=5.2, va="center",
                 color=F.GREY)
    axA.set_yticks(y); axA.set_yticklabels(targets, fontsize=6.2)
    axA.set_xlim(40, 122); axA.set_xticks([50, 75, 100])
    axA.set_xlabel("sensitivity (%, 95% CI)", fontsize=7.2)
    axA.legend(frameon=False, fontsize=5.6, loc="upper center",
               bbox_to_anchor=(0.5, -0.30), ncol=2, handlelength=1.2)
    F.bare(axA); F.panel(axA, "B", dx=-52)

    # B 特異度：draft 非帶菌者與非目標物種
    fp = [ext["per_target"][t]["draft_non_carriers"]
          - round(ext["per_target"][t]["specificity"]
                  * ext["per_target"][t]["draft_non_carriers"]) for t in targets]
    x = np.arange(len(targets))
    # 非目標物種那一列全是 0，畫成長條就是看不見的圖元：拿掉，改成一行文字
    axB.bar(x, fp, 0.6, color=F.NAVY, edgecolor="white", linewidth=0.4,
            label=f"draft non-carriers "
                  f"(n\u2248{ext['per_target'][targets[0]]['draft_non_carriers']})")
    axB.set_xticks(x); axB.set_xticklabels([t.split("-")[0] for t in targets],
                                           fontsize=5.6, rotation=90)
    axB.tick_params(axis="x", pad=3)
    axB.set_ylabel("genomes amplified", fontsize=7.2)
    axB.set_ylim(0, 4); axB.set_yticks([0, 1, 2, 3])
    axB.legend(frameon=False, fontsize=5.2, loc="upper center", ncol=1,
               handlelength=1.1)
    axB.text(0.5, -0.42, "no target amplified any non-target genome",
             transform=axB.transAxes, ha="center", fontsize=5.4, color=F.GREY)
    F.bare(axB); F.panel(axB, "C", dx=-40)

    # C 家族盛行率：設計集 vs draft（配對點圖，比散點好讀也不會疊標籤）
    fams = sorted({t.split("-")[0] for t in targets})
    y = np.arange(len(fams))[::-1]
    for yi, fam in zip(y, fams):
        tgt = next(t for t in targets if t.startswith(fam))
        sp = collections.Counter(r["species"] for r in mat if r.get(fam))
        main_sp = sp.most_common(1)[0][0] if sp else ""
        tot = sum(1 for r in mat if r["species"] == main_sp)
        p_panel = 100 * sum(1 for r in mat if r.get(fam)) / max(1, tot)
        d = [r for r in rows if r["group"] == "draft CoNS" and r["species"] == main_sp]
        p_draft = 100 * sum(int(r[f"{tgt}_carrier"]) for r in d) / max(1, len(d))
        axC.plot([p_panel, p_draft], [yi, yi], color=F.FAINT, lw=1.0, zorder=1)
        axC.scatter(p_panel, yi, s=18, color=F.GREY, marker="o", zorder=3,
                    label="design panel" if yi == y[0] else "")
        axC.scatter(p_draft, yi, s=18, color=F.NAVY, marker="D", zorder=3,
                    label="held-out drafts" if yi == y[0] else "")
        axC.text(104, yi, f"S. {main_sp[:4]}.", fontsize=5.0, va="center", color=F.GREY)
    axC.set_yticks(y); axC.set_yticklabels(fams, fontsize=6.0)
    axC.set_xlim(0, 126); axC.set_xticks([0, 50, 100])
    axC.set_xlabel("carriage within the family's main species (%)", fontsize=7)
    axC.legend(frameon=False, fontsize=5.4, loc="upper center",
               bbox_to_anchor=(0.5, -0.30), ncol=2, handlelength=1.1)
    F.bare(axC); F.panel(axC, "D", dx=-44)

    F.save(fig, "Figure4_assay")


if __name__ == "__main__":
    main()
