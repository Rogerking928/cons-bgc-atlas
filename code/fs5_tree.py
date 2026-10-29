#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 8：554 株的相似度樹（全蛋白體 MinHash）＋ 基因含量圈。"""
import collections, csv, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import figstyle as F
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
GAP = 0.62            # radians left open at the top for the ring labels
RINGS = [("F001", "F001 aureusimine (epi)", "#31485d"),
         ("F004", "F004 gallidermin-family", "#b8791f"),
         ("F006", "F006 lugdunin", "#a4303f"),
         ("F009", "F009 lacticin 481 (plasmid)", "#6b5b95"),
         ("mecA", "mecA", "#4f7d5b"),
         ("cand", "candidate precursor (model)", "#1f2933")]


def main():
    D = np.load(os.path.join(OUT, "tree_D.npy"))
    accs = json.load(open(os.path.join(OUT, "tree_meta.json")))["accs"]
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    amr = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv")))}
    # 模型那一圈要等 s19 跑完；還沒有就畫成空圈，不要讓整張圖跑不出來
    cand = collections.defaultdict(int)
    p_rep = os.path.join(OUT, "precursor_report.json")
    p_scan = os.path.join(OUT, "precursor_scan.csv")
    if os.path.exists(p_rep) and os.path.exists(p_scan):
        thr = json.load(open(p_rep))["threshold_999th_percentile_of_negatives"]
        for r in csv.DictReader(open(p_scan)):
            if float(r["score"]) >= thr and not r["known_domain"]:
                cand[r["assembly"]] += 1

    d = (D + D.T) / 2
    np.fill_diagonal(d, 0)
    Z = linkage(squareform(d, checks=False), method="average")
    dd = dendrogram(Z, no_plot=True)
    order = dd["leaves"]
    n = len(order)
    pos = {leaf: i for i, leaf in enumerate(order)}

    span = 2 * np.pi - GAP
    def ang(x):                       # dendrogram x coordinate -> angle
        return GAP / 2 + (x - 5.0) / (10.0 * (n - 1)) * span
    rmax = d.max()
    def rad(y):                       # distance -> radius (root at the centre)
        return 1.0 - y / rmax

    fig, ax = plt.subplots(figsize=(F.WIDTH_IN, F.WIDTH_IN * 1.02),
                           subplot_kw={"projection": "polar"})
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    sp_of = {a: mat[a]["species"] for a in accs}
    # colour a branch by species when all its leaves agree
    leaf_species = [sp_of[accs[i]] for i in order]
    for xs, ys in zip(dd["icoord"], dd["dcoord"]):
        a0, a1 = ang(xs[0]), ang(xs[3])
        i0 = int(round((xs[0] - 5.0) / 10.0))
        i1 = int(round((xs[3] - 5.0) / 10.0))
        seg = leaf_species[max(0, min(i0, i1)):max(i0, i1) + 1] or ["?"]
        col = F.sp_color(seg[0]) if len(set(seg)) == 1 else F.GREY
        ax.plot([a0, a0], [rad(ys[0]), rad(ys[1])], color=col, lw=0.45)
        ax.plot([a1, a1], [rad(ys[3]), rad(ys[2])], color=col, lw=0.45)
        t = np.linspace(a0, a1, 24)
        ax.plot(t, [rad(ys[1])] * 24, color=col, lw=0.45)

    # rings
    r0 = 1.035
    width = 0.038
    for k, (key, lab, col) in enumerate(RINGS):
        rr = r0 + k * (width + 0.012)
        ax.bar(np.linspace(GAP / 2, GAP / 2 + span, n), [width] * n,
               width=span / n, bottom=rr, color="#eceff1", edgecolor="none")
        marks, angles = [], []
        for i, leaf in enumerate(order):
            a = accs[leaf]
            if key.startswith("F"):
                on = bool(mat[a][key])
            elif key == "cand":
                on = cand.get(a, 0) > 0
            else:
                on = float(amr.get(a, {}).get(key, 0) or 0) >= 0.5
            if on:
                marks.append(width); angles.append(GAP / 2 + (i + 0.5) * span / n)
        if marks:
            ax.bar(angles, marks, width=span / n, bottom=rr, color=col,
                   edgecolor="none")
        ax.text(GAP / 2 - 0.02, rr + width / 2, lab, fontsize=6.4, ha="right",
                va="center", color=col)

    ax.set_ylim(0, r0 + len(RINGS) * (width + 0.012) + 0.02)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines["polar"].set_visible(False)
    ax.grid(False)
    counts = collections.Counter(sp_of.values())
    ax.text(0, 0, f"{len(accs)} genomes\n{len(counts)} species",
            ha="center", va="center", fontsize=7.5, fontweight="bold")
    top = [s for s in F.SPECIES_ORDER if s in counts][:8]
    handles = [Patch(facecolor=F.sp_color(s), label=f"S. {s} ({counts[s]})")
               for s in top]
    ax.legend(handles=handles, frameon=False, fontsize=6.4, ncol=4,
              loc="lower center", bbox_to_anchor=(0.5, -0.05))
    ax.set_title("average linkage on whole-proteome MinHash distances",
                 fontsize=8.5, pad=14)
    F.save(fig, "FigureS5_tree")


if __name__ == "__main__":
    main()
