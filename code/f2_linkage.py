#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 2：基因簇與後天抗藥基因的基因體鄰接（座位圖、距離對虛無、可動性、GC）。"""
import collections, csv, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrow, Patch

import figstyle as F
import names
from f3_supp_maps import AMR_RE, COLOR, klass, short   # 座位圖的畫法共用一份

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SHOW = [("C00115", "plasmid"), ("C00077", "chromosome")]   # 質體、染色體各一


def draw_map(ax, genes, kb_per_unit, letter, title):
    lo = int(genes[0]["start"])
    span = int(genes[-1]["end"]) - lo
    row, last = 0, ("", 0)
    for g in genes:
        s, e = int(g["start"]) - lo, int(g["end"]) - lo
        k = klass((g["domains"] or "").split("-"), g["product"])
        x, dx = (s, e - s) if g["strand"] == "+" else (e, s - e)
        ax.add_patch(FancyArrow(x, 0, dx, 0, width=0.34, head_width=0.52,
                                length_includes_head=True,
                                head_length=min(abs(dx) * 0.35, span * 0.02),
                                facecolor=COLOR[k], edgecolor="white", linewidth=0.5))
        lab = short(g["product"], k not in ("other", "export"))
        if lab and lab == last[0] and (s - last[1]) < span * 0.12:
            lab = ""
        elif lab:
            last = (lab, e)
        if lab:
            ax.text((s + e) / 2, 0.44 + (row % 4) * 0.82, lab, ha="center",
                    va="bottom", fontsize=5.6, color=F.INK)
            row += 1
    ax.set_xlim(-1000, kb_per_unit * 1000)
    ax.set_ylim(-0.6, 4.4)
    ax.set_yticks([])
    ticks = [0, 10000, 20000, 30000]
    ax.set_xticks([t for t in ticks if t <= kb_per_unit * 1000])
    ax.set_xticklabels([f"{t // 1000}" for t in ticks if t <= kb_per_unit * 1000],
                       fontsize=6.4)
    ax.set_xlabel("kb", fontsize=6.4)
    ax.set_title(title, loc="left", fontsize=7.4, pad=18)
    F.bare(ax, keep=("bottom",))
    F.panel(ax, letter, dx=-18, dy=2)


def main():
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = {r["cluster"]: r for r in csv.DictReader(
        open(os.path.join(OUT, "clusters.tsv")), delimiter="\t")}
    cds = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"):
        cds[r["cluster"]].append(r)
    loc = list(csv.DictReader(open(os.path.join(OUT, "colocalisation.csv"))))
    fams = list(csv.DictReader(open(os.path.join(OUT, "families.tsv")), delimiter="\t"))
    keep = [f["family"] for f in fams if int(f["n_assemblies"]) >= 4]
    kind = {f["family"]: f["kind"] for f in fams}
    nullres = json.load(open(os.path.join(OUT, "distance_null.json")))
    gc = list(csv.DictReader(open(os.path.join(OUT, "cluster_gc.csv"))))
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}

    fig = plt.figure(figsize=(F.WIDTH_IN, 8.2))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 1.5, 1.25],
                          hspace=0.95, wspace=0.34)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, :])
    axC = fig.add_subplot(gs[2, 0])
    axD = fig.add_subplot(gs[2, 1])
    axE = fig.add_subplot(gs[3, 0])
    axG = fig.add_subplot(gs[3, 1])

    # A, B 兩張座位圖，共用同一個 kb 尺度
    longest = max(int(cl[c]["window_end"]) - int(cl[c]["window_start"])
                  for c, _ in SHOW) / 1000
    for ax, (cid, lab), letter in zip((axA, axB), SHOW, "AB"):
        g = sorted(cds[cid], key=lambda r: int(r["start"]))
        c = cl[cid]
        sp = mat[c["assembly"]]["species"]
        draw_map(ax, g, longest, letter,
                 f"{fam_of[cid]}  {names.fam_short(fam_of[cid])} — S. {sp}, "
                 f"{c['replicon_type']}, core {int(c['span']) / 1000:.1f} kb")
    drawn = {klass((g["domains"] or "").split("-"), g["product"])
             for cid, _ in SHOW for g in cds[cid]}
    handles = [Patch(facecolor=COLOR[k], label=k) for k in
               ("structural", "maturation", "export", "regulation", "mobility",
                "resistance", "other") if k in drawn]
    axB.legend(handles=handles, frameon=False, fontsize=6.2, ncol=7,
               loc="upper center", bbox_to_anchor=(0.5, -0.45))

    # C 距離 vs 隨機重擺
    rng = np.random.default_rng(20260920)
    n_of = {}
    for i, g in enumerate(("chromosome", "plasmid")):
        v = np.array([max(int(r["nearest_bp"]) / 1000, 0.1) for r in loc
                      if r["nearest_bp"] != ""
                      and (r["replicon_type"] != "chromosome") == (g == "plasmid")])
        n_of[g] = len(v)
        axC.scatter(rng.normal(i, 0.07, len(v)), v, s=7, alpha=0.6,
                    color=F.REPLICON[g], edgecolor="none", zorder=3)
        axC.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, color=F.INK, lw=1.6, zorder=4)
        nl = nullres["groups"][g]["null_median_kb"]
        axC.plot([i - 0.32, i + 0.32], [nl] * 2, color=F.RED, lw=1.4, ls="--", zorder=5)
        p = nullres["groups"][g]["wilcoxon_p"]
        # 1400 還在散點雲裡（染色體那組上探到 ~2000 kb），字會壓在點上；
        # 抬到資料上方，軸上緣跟著放寬
        axC.text(i, 9000, f"obs {np.median(v):.0f}\nnull {nl:.0f} kb\nP = {p:.3g}",
                 ha="center", fontsize=5.8, color=F.INK)
    axC.set_yscale("log"); axC.set_ylim(0.06, 90000)
    axC.set_yticks([0.1, 1, 10, 100, 1000])
    axC.set_yticklabels(["≤0.1", "1", "10", "100", "1000"], fontsize=6.4)
    axC.set_xticks([0, 1])
    axC.set_xticklabels([f"chromosome\n(n={n_of['chromosome']})",
                         f"plasmid\n(n={n_of['plasmid']})"], fontsize=6.8)
    axC.set_xlim(-0.6, 1.6)
    axC.set_ylabel("distance to nearest acquired\nresistance gene (kb)", fontsize=7)
    axC.plot([], [], color=F.INK, lw=1.6, label="observed median")
    axC.plot([], [], color=F.RED, lw=1.4, ls="--",
             label="expected if placed at random")
    axC.legend(frameon=False, fontsize=5.4, loc="lower left", handlelength=1.6,
               borderpad=0.1)
    F.bare(axC); F.panel(axC, "C", dx=-46)

    # D 各家族共 replicon 比例
    frac, labs, cols = [], [], []
    for fam in keep:
        rs = [r for r in loc if r["family"] == fam]
        if not rs:
            continue
        frac.append(100 * sum(1 for r in rs if int(r["n_acquired_on_replicon"]) > 0) / len(rs))
        labs.append(fam); cols.append(F.KIND[kind[fam]])
    y = np.arange(len(labs))[::-1]
    axD.barh(y, frac, color=cols, edgecolor="white", linewidth=0.4, height=0.72)
    axD.set_yticks(y); axD.set_yticklabels([l[1:] for l in labs], fontsize=5.4)
    axD.set_xlim(0, 100); axD.set_xticks([0, 50, 100])
    axD.set_xlabel("clusters sharing a replicon with an\nacquired determinant (%), by family F0..",
                   fontsize=6.6)
    axD.tick_params(labelsize=6.4)
    F.bare(axD); F.panel(axD, "D", dx=-34)

    # E 島內可動基因（堆疊比例）
    kinds = ("NRPS/PKS", "hybrid", "RiPP")
    cats = [("0", F.FAINT, ""), ("1", F.GREY, "///"), ("≥2", F.NAVY, "xxx")]
    bottom = np.zeros(len(kinds))
    counts = {k: collections.Counter() for k in kinds}
    for c in cl.values():
        m = int(c["mobile_cds"])
        counts[c["kind"]]["0" if m == 0 else "1" if m == 1 else "≥2"] += 1
    for lab, col, hh in cats:
        v = np.array([100 * counts[k][lab] / max(1, sum(counts[k].values())) for k in kinds])
        axE.bar(np.arange(len(kinds)), v, 0.58, bottom=bottom, color=col, hatch=hh,
                edgecolor="white", linewidth=0.4, label=lab)
        bottom += v
    axE.set_xticks(range(len(kinds)))
    axE.set_xticklabels([f"{k}\n(n={sum(counts[k].values())})" for k in kinds], fontsize=6.4)
    axE.set_ylabel("clusters (%)", fontsize=7)
    axE.set_ylim(0, 100); axE.set_yticks([0, 50, 100])
    axE.legend(frameon=False, fontsize=6, ncol=3, title="mobility genes in window",
               title_fontsize=6, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    F.bare(axE); F.panel(axE, "E", dx=-40, dy=16)

    # F GC 差
    gfams = [f for f in keep if any(r["family"] == f for r in gc)]
    for i, fam in enumerate(gfams):
        v = np.array([float(r["delta"]) for r in gc if r["family"] == fam])
        axG.scatter(rng.normal(i, 0.08, len(v)), v, s=6, alpha=0.55, color=F.NAVY,
                    edgecolor="none")
        axG.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color=F.RED, lw=1.2)
    axG.axhline(0, color=F.GREY, lw=0.8, ls=":")
    axG.set_xticks(range(len(gfams)))
    axG.set_xticklabels([f[1:] for f in gfams], fontsize=5.4, rotation=90)
    axG.tick_params(axis="x", pad=3)
    axG.set_ylabel("cluster GC − genome GC\n(percentage points)", fontsize=7)
    axG.set_yticks([-10, -5, 0, 5])
    F.bare(axG); F.panel(axG, "F", dx=-44)

    F.save(fig, "Figure2_linkage")


if __name__ == "__main__":
    main()
