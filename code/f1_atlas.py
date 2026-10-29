#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 1：圖譜（攜帶率 ＋ 家族 × 物種熱圖 ＋ 位置 ＋ 島長，共用同一個 y 軸）。"""
import collections, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
CMAP = LinearSegmentedColormap.from_list("cons", ["#f2f5f7", F.NAVY])


def main():
    mat = list(csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv"))))
    fams = list(csv.DictReader(open(os.path.join(OUT, "families.tsv")), delimiter="\t"))
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    keep = [f["family"] for f in fams if int(f["n_assemblies"]) >= 4]
    kind = {f["family"]: f["kind"] for f in fams}
    sp_of = {r["assembly"]: r["species"] for r in mat}
    counts = collections.Counter(sp_of.values())
    carriers = {r["species"] for r in mat if any(r.get(f) for f in keep)}
    order = [s for s in F.SPECIES_ORDER if s != "other" and s in carriers]
    order += sorted(carriers - set(order))

    kinds = collections.defaultdict(lambda: collections.defaultdict(set))
    for c in cl:
        kinds[sp_of[c["assembly"]]][c["kind"]].add(c["assembly"])

    fig = plt.figure(figsize=(F.WIDTH_IN, 6.0))
    # 色階條自己一列：掛在 axB 上會把 B 縮短，C／D 的列就對不齊了
    # 原本的 chr/plasmid 長條每根都是滿格（堆疊比例），等於沒有在編碼資料：
    # 改成家族標籤旁邊一個小方塊，版面讓給熱圖與島長（2026-09-20 審閱意見）
    gs = fig.add_gridspec(3, 2, height_ratios=[0.85, 2.0, 0.10],
                          width_ratios=[2.9, 0.95], hspace=0.72, wspace=0.08)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1], sharey=axB)
    cax = fig.add_subplot(gs[2, 0])

    # ---------------------------------------------------------------- A 攜帶率
    kk = ["NRPS/PKS", "hybrid", "RiPP"]
    allsp = [s for s in F.SPECIES_ORDER if s != "other" and counts.get(s, 0) >= 20]
    x = np.arange(len(allsp))
    w = 0.26
    for i, k in enumerate(kk):
        v = [100 * len(kinds[s][k]) / counts[s] for s in allsp]
        axA.bar(x + (i - 1) * w, v, w, color=F.KIND[k],
                hatch=["", "xx", ".."][i], edgecolor="white", linewidth=0.4, label=k)
    for xi, s in zip(x, allsp):
        if not any(kinds[s][k] for k in kk):
            axA.text(xi, 3, "none detected", ha="center", va="bottom", fontsize=5.4,
                     color=F.GREY, rotation=90)
    axA.set_xticks(x)
    # 只有物種縮寫是學名，下面那行 n= 不是；整個標籤設 style="italic" 會把
    # 數字也斜掉，交給 figstyle 的斜體掃描（字彙認得 "epid." 這種縮寫）
    axA.set_xticklabels([f"{s[:4]}.\nn={counts[s]}" for s in allsp], fontsize=6.0)
    axA.set_ylabel("genomes carrying\na cluster (%)", fontsize=7.5)
    axA.set_ylim(0, 100); axA.set_yticks([0, 50, 100])
    axA.legend(frameon=False, fontsize=6.6, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), handlelength=1.3, columnspacing=1.2)
    F.bare(axA); F.panel(axA, "A", dx=-56, dy=14)

    # 家族在染色體還是質體：標籤前面一個小方塊
    loc_of = {}
    for fam in keep:
        cs = [c for c in cl if fam_of[c["cluster"]] == fam]
        p_ = sum(1 for c in cs if c["replicon_type"] != "chromosome")
        loc_of[fam] = ("plasmid" if p_ == len(cs) else
                       "chromosome" if p_ == 0 else "both")

    # -------------------------------------------------------------- B 家族熱圖
    M = np.zeros((len(keep), len(order)))
    N = np.zeros_like(M, dtype=int)
    for i, fam in enumerate(keep):
        for j, sp in enumerate(order):
            n = sum(1 for r in mat if r["species"] == sp and r.get(fam))
            N[i, j] = n
            M[i, j] = 100 * n / counts[sp] if n else np.nan
    im = axB.imshow(M, cmap=CMAP, vmin=0, vmax=100, aspect="auto")
    for i in range(len(keep)):
        for j in range(len(order)):
            if N[i, j]:
                axB.text(j, i, str(N[i, j]), ha="center", va="center", fontsize=5.6,
                         color="white" if M[i, j] > 55 else F.INK)
    axB.set_xticks(range(len(order)))
    # 旋轉標籤的外框比看起來大：用 4 字母縮寫（圖說列全名）才不會互相壓到
    axB.set_xticklabels([s[:4] + "." for s in order], fontsize=6.2,
                        rotation=45, ha="right")
    axB.set_yticks(range(len(keep)))
    axB.set_yticklabels([f"{f}  {names.fam_short(f)}" for f in keep], fontsize=6.4)
    for i, fam in enumerate(keep):
        col = {"chromosome": F.REPLICON["chromosome"], "plasmid": F.REPLICON["plasmid"],
               "both": F.GREY}[loc_of[fam]]
        axB.add_patch(plt.Rectangle((-0.92, i - 0.22), 0.34, 0.44, color=col,
                                    clip_on=False, transform=axB.transData))
    axB.set_xlim(-1.0, len(order) - 0.5)
    axB.tick_params(length=0)
    for sp in axB.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[0, 50, 100])
    cb.set_label("carried by this % of the species' genomes (number in cell = genomes)",
                 fontsize=6.0)
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_visible(False)
    hnd = [plt.Rectangle((0, 0), 1, 1, color=F.REPLICON["chromosome"]),
           plt.Rectangle((0, 0), 1, 1, color=F.REPLICON["plasmid"]),
           plt.Rectangle((0, 0), 1, 1, color=F.GREY)]
    axB.legend(hnd, ["chromosome", "plasmid", "both"], frameon=False, fontsize=5.6,
               ncol=3, loc="lower left", bbox_to_anchor=(0.0, 1.0), handlelength=0.9,
               columnspacing=0.8)
    F.panel(axB, "B", dx=-158, dy=16)

    # ------------------------------------------------------------- C 島長
    y = np.arange(len(keep))
    for yi, fam in zip(y, keep):
        v = [int(c["span"]) / 1000 for c in cl if fam_of[c["cluster"]] == fam]
        axD.plot([min(v), max(v)], [yi, yi], color=F.GREY, lw=1.0, zorder=2)
        axD.scatter(np.median(v), yi, color=F.KIND[kind[fam]], s=16, zorder=3)
    axD.set_xlabel("core span (kb)", fontsize=6.4)
    axD.set_xlim(0, 32); axD.set_xticks([0, 15, 30])
    axD.tick_params(labelsize=6, left=False)
    plt.setp(axD.get_yticklabels(), visible=False)
    F.bare(axD, keep=("bottom",)); F.panel(axD, "C", dx=-6, dy=16)
    axB.set_ylim(len(keep) - 0.5, -0.5)   # 共用 y 軸，設一次就好

    F.save(fig, "Figure1_atlas")


if __name__ == "__main__":
    main()
