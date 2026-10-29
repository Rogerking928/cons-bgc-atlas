#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 3：帶菌是被「怎麼呼叫基因簇」和譜系綁住的，不是生物學訊號。"""
import collections, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
CALL = 0.5
DEFS = [("strict (>=2 domain-carrying CDS)", "calling rule v1:\n≥2 domain-carrying CDS",
         F.GREY),
        ("fixed (a lone megasynthetase counts)", "calling rule v2:\na lone synthetase counts",
         F.NAVY),
        ("gene present (RefSeq annotation)", "gene present\n(RefSeq annotation)", F.AMBER)]


def main():
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    amr = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv")))}
    mlst = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "mlst.csv")))}
    rule = json.load(open(os.path.join(OUT, "calling_rule.json")))

    fig = plt.figure(figsize=(F.WIDTH_IN, 4.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 1.0], width_ratios=[1.3, 1.0],
                          hspace=0.9, wspace=0.42)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    sts = list(rule["by_ST"][DEFS[0][0]])
    x = np.arange(len(sts)); w = 0.27
    for i, (key, lab, col) in enumerate(DEFS):
        v = [100 * rule["by_ST"][key][s]["carriers"] / rule["by_ST"][key][s]["n"]
             for s in sts]
        axA.bar(x + (i - 1) * w, v, w, color=col, hatch=["", "///", "..."][i],
                edgecolor="white", linewidth=0.4, label=lab)
        # 0 高度的長條就是不存在，而 ST23 在 v1 規則下正是 0——那是整張圖的重點
        for xi, vv in zip(x, v):
            if vv == 0:
                axA.text(xi + (i - 1) * w, 1.5, "0", ha="center", va="bottom",
                         fontsize=5.6, color=col, fontweight="bold")
    axA.set_xticks(x)
    axA.set_xticklabels([f"{s}\nn={rule['by_ST'][DEFS[0][0]][s]['n']}" for s in sts],
                        fontsize=6.2)
    axA.set_ylabel("S. epidermidis genomes\ncarrying F001 (%)", fontsize=7.5)
    axA.set_ylim(0, 100); axA.set_yticks([0, 50, 100])
    axA.legend(frameon=False, fontsize=6.0, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), handlelength=1.2, columnspacing=1.2)
    F.bare(axA); F.panel(axA, "A", dx=-48, dy=26)

    dets = ["erm(A)", "ant(9)-Ia", "mecA", "blaZ"]
    y = np.arange(len(dets))[::-1]
    for i, (key, lab, col) in enumerate(DEFS):
        off = (1 - i) * 0.22
        for yi, det in zip(y, dets):
            r = rule["association"][key][det]
            o = r["odds"] if np.isfinite(r["odds"]) and r["odds"] > 0 else r["odds_haldane"]
            lo, hi = r.get("ci", [None, None])
            if lo and hi:
                axB.plot([lo, hi], [yi + off] * 2, color=col, lw=0.9, zorder=2)
            axB.scatter(o, yi + off, s=18, color=col, marker=["o", "s", "D"][i], zorder=3,
                        label=lab.replace("\n", " ") if yi == y[0] else "")
    axB.axvline(1, color=F.GREY, lw=0.8, ls=":")
    axB.set_xscale("log"); axB.set_xlim(0.01, 100)
    axB.set_yticks(y); axB.set_yticklabels(dets, fontsize=6.6)
    axB.set_xlabel("odds ratio, F001 against the determinant (95% CI)", fontsize=7.2)
    axB.legend(frameon=False, fontsize=5.4, loc="upper center",
               bbox_to_anchor=(0.5, -0.34), ncol=2, handlelength=1.1)
    e1 = rule["association"][DEFS[0][0]]["erm(A)"]
    e2 = rule["association"][DEFS[1][0]]["erm(A)"]
    axB.set_title(f"erm(A): OR {e1['odds']:.2f} (P = {e1['p']:.0e}) under v1,\n"
                  f"{e2['odds']:.2f} (P = {e2['p']:.2f}) under v2",
                  fontsize=6.2, loc="left")
    F.bare(axB); F.panel(axB, "B", dx=-44)

    lug = [a for a, r in mat.items() if r["species"] == "lugdunensis"]
    lsts = collections.Counter(mlst.get(a, {}).get("ST", "") for a in lug)
    lkeep = [s for s, n in lsts.most_common() if s and n >= 3]
    lg = [[a for a in lug if mlst.get(a, {}).get("ST") == s] for s in lkeep]
    lg.append([a for a in lug if mlst.get(a, {}).get("ST") not in lkeep])
    llab = [f"ST{s}" for s in lkeep] + ["other/\nuntyped"]
    x = np.arange(len(lg)); w = 0.27
    for i, (fam, col, hh) in enumerate((("F003", F.NAVY, ""), ("F006", F.RED, "///"),
                                        ("F005", F.PURPLE, "..."))):
        v = [100 * sum(1 for a in g if mat[a].get(fam)) / len(g) for g in lg]
        axC.bar(x + (i - 1) * w, v, w, color=col, hatch=hh, edgecolor="white",
                linewidth=0.4, label=f"{fam} {names.fam_short(fam)}")
    axC.set_xticks(x)
    axC.set_xticklabels([f"{l}\nn={len(g)}" for l, g in zip(llab, lg)], fontsize=6.0)
    axC.set_ylabel("genomes (%)", fontsize=7.5)
    axC.set_ylim(0, 100); axC.set_yticks([0, 50, 100])
    # 只有學名要斜體，"by sequence type" 不要：交給斜體掃描
    axC.set_title("S. lugdunensis by sequence type", fontsize=7.6,
                  loc="left", pad=20)
    axC.legend(frameon=False, fontsize=5.2, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 1.0), handlelength=1.0, columnspacing=0.5)
    F.bare(axC); F.panel(axC, "C", dx=-40, dy=20)

    F.save(fig, "Figure3_lineage")


if __name__ == "__main__":
    main()
