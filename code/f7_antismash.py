#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 7：圖譜與 antiSMASH 的逐座對照，以及它追出來的假基因化合成酶。

審稿人（AoM #42 Reviewer 1）要求與 antiSMASH 比較。比較只有在範圍對齊時才誠實，
所以四格的分工是：

  A  antiSMASH 呼叫的 3,339 個 region 按類別分成「圖譜涵蓋」與「圖譜不收」，
     說明圖譜是什麼、不是什麼（最常見的四類全是設計上不收的）。
  B  兩個方向的召回：圖譜的 562 個簇有多少被範圍內的 region 蓋到，
     範圍內的 region 有多少被圖譜的簇蓋到。
  C  邊界：同一個座位，antiSMASH 的 region 與圖譜的核心跨距比較。
  D  對照追出來的結果：圖譜判為非帶菌、但 RefSeq 把合成酶標成假基因的基因體，
     亦即「基因座存在」與「合成酶完整」是兩個數字。
"""
import collections, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")


def load():
    comp = list(csv.DictReader(open(os.path.join(OUT, "antismash_compare.csv"))))
    regs = list(csv.DictReader(open(os.path.join(OUT, "antismash_regions.csv"))))
    summ = json.load(open(os.path.join(OUT, "antismash_compare.json")))
    pseudo = json.load(open(os.path.join(OUT, "pseudogene_summary.json")))
    return comp, regs, summ, pseudo


def panel_a(ax, regs):
    """antiSMASH 的 region 按 product 分組，標明哪些在圖譜範圍內。"""
    cnt = collections.Counter()
    scope = {}
    for r in regs:
        # 多 product 的 region 以第一個命名，範圍標記用 s55 算好的那個
        p = r["product"].split(";")[0]
        cnt[p] += 1
        scope[p] = r["scope"]
    items = cnt.most_common(11)
    y = np.arange(len(items))[::-1]
    for (p, n), yy in zip(items, y):
        inside = scope[p] == "in"
        ax.barh(yy, n, height=0.72,
                color=F.NAVY if inside else F.FAINT,
                edgecolor=F.INK if inside else F.GREY, linewidth=0.5,
                hatch="" if inside else "///")
        ax.text(n + 18, yy, str(n), va="center", fontsize=6, color=F.INK)
    ax.set_yticks(y)
    ax.set_yticklabels([p for p, _ in items], fontsize=6)
    ax.set_xlabel("antiSMASH regions in the 554 genomes", fontsize=7)
    ax.set_xlim(0, max(n for _, n in items) * 1.18)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, facecolor=F.NAVY,
                                     edgecolor=F.INK, linewidth=0.5,
                                     label="within the atlas's profile scope"),
                       plt.Rectangle((0, 0), 1, 1, facecolor=F.FAINT,
                                     edgecolor=F.GREY, linewidth=0.5,
                                     hatch="///", label="outside it")],
              fontsize=5.8, loc="lower right", frameon=False)
    F.bare(ax)


def panel_b(ax, summ):
    """兩個方向的召回，各畫成一條 100% 的棒子。"""
    a_tot = summ["atlas_clusters"]
    a_hit = summ["atlas_clusters_overlapped_by_in_scope_region"]
    r_tot = summ["antismash_regions_in_scope"]
    r_hit = summ["in_scope_regions_matched_by_atlas"]
    bars = [("atlas clusters\nfound by antiSMASH", a_hit, a_tot),
            ("in-scope antiSMASH regions\nfound by the atlas", r_hit, r_tot)]
    for i, (lab, hit, tot) in enumerate(bars):
        ax.barh(i, 100.0 * hit / tot, height=0.5, color=F.NAVY,
                edgecolor=F.INK, linewidth=0.5)
        ax.barh(i, 100.0, height=0.5, color="none", edgecolor=F.GREY,
                linewidth=0.5, zorder=3)
        ax.text(100.0 * hit / tot - 1.5, i, f"{hit}/{tot}", va="center",
                ha="right", fontsize=6.4, color="white", zorder=4)
    ax.set_yticks(range(len(bars)))
    ax.set_yticklabels([b[0] for b in bars], fontsize=6.2)
    ax.set_xlim(0, 104)
    ax.set_xlabel("per cent", fontsize=7)
    ax.invert_yaxis()
    F.bare(ax)


def panel_c(ax, comp):
    """同一座位的跨距：圖譜的核心 vs antiSMASH 的 region。"""
    pairs = [(int(r["cluster_span"]), int(r["as_span"]))
             for r in comp if r["as_scope"] == "in" and r["as_span"]]
    cl = np.array([p[0] for p in pairs]) / 1000.0
    rg = np.array([p[1] for p in pairs]) / 1000.0
    ax.scatter(cl, rg, s=7, color=F.NAVY, alpha=0.45, linewidths=0)
    hi = max(cl.max(), rg.max()) * 1.05
    ax.plot([0, hi], [0, hi], color=F.GREY, linewidth=0.7, linestyle="--")
    ax.set_xlim(0, hi)
    ax.set_ylim(0, hi)
    ax.set_xlabel("atlas core span (kb)", fontsize=7)
    ax.set_ylabel("antiSMASH region span (kb)", fontsize=7)
    # 點全在對角線上方（region 一律比核心寬），所以註記放右下角的空白處
    ax.text(0.42, 0.05, f"median {np.median(rg):.0f} kb vs {np.median(cl):.1f} kb"
                        f"\n({np.median(rg) / np.median(cl):.1f}x wider)",
            transform=ax.transAxes, fontsize=6, va="bottom", ha="left",
            color=F.INK)
    F.bare(ax)


def panel_d(ax, pseudo):
    """基因座存在 vs 合成酶完整。只畫有假基因命中的家族。"""
    fams = [(f, d) for f, d in pseudo["families"].items()
            if d["genomes_with_pseudogene_but_not_called_carrier"]]
    fams.sort(key=lambda x: -x[1]["locus_present_including_pseudogenised"])
    y = np.arange(len(fams))
    for i, (f, d) in enumerate(fams):
        intact = d["atlas_carriers"]
        extra = d["genomes_with_pseudogene_but_not_called_carrier"]
        ax.barh(i, intact, height=0.6, color=F.NAVY, edgecolor=F.INK,
                linewidth=0.5, label="intact synthetase (atlas call)" if not i else None)
        ax.barh(i, extra, left=intact, height=0.6, color=F.AMBER,
                edgecolor=F.INK, linewidth=0.5, hatch="///",
                label="locus present, gene pseudogenised" if not i else None)
        ax.text(intact + extra + 1.5, i, f"{intact}+{extra}", va="center",
                fontsize=6, color=F.INK)
    ax.set_yticks(y)
    ax.set_yticklabels([names.fam_short(f) for f, _ in fams], fontsize=6)
    ax.set_xlabel("genomes of the family's species", fontsize=7)
    ax.set_xlim(0, max(d["locus_present_including_pseudogenised"]
                       for _, d in fams) * 1.22)
    ax.invert_yaxis()
    ax.legend(fontsize=5.8, loc="lower right", frameon=False)
    F.bare(ax)


def main():
    comp, regs, summ, pseudo = load()
    fig = plt.figure(figsize=(F.WIDTH_IN, 5.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0],
                          hspace=0.52, wspace=0.42,
                          left=0.155, right=0.975, top=0.965, bottom=0.085)
    ax = fig.add_subplot(gs[0, 0]); panel_a(ax, regs); F.panel(ax, "A", dx=-96)
    ax = fig.add_subplot(gs[0, 1]); panel_b(ax, summ); F.panel(ax, "B", dx=-92)
    ax = fig.add_subplot(gs[1, 0]); panel_c(ax, comp); F.panel(ax, "C", dx=-48)
    ax = fig.add_subplot(gs[1, 1]); panel_d(ax, pseudo); F.panel(ax, "D", dx=-62)
    F.save(fig, "Figure7_antismash")


if __name__ == "__main__":
    main()
