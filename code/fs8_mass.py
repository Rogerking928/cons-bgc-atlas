#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 7：臨床實驗室接口（PCR 分型方案 ＋ MALDI-TOF 可查的理論質量）。"""
import collections, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")


def best_primers():
    best = {}
    for r in csv.DictReader(open(os.path.join(OUT, "primers.csv"))):
        k = r["target"]
        sc = (float(r["sensitivity"]), float(r["specificity"]), -int(r["product_bp"]))
        if k not in best or sc > best[k][0]:
            best[k] = (sc, r)
    return [r for _, r in sorted(best.values(), key=lambda v: v[1]["target"])]


def main():
    mass = list(csv.DictReader(open(os.path.join(OUT, "precursor_masses.csv"))))
    fig, (axC, axD) = plt.subplots(1, 2, figsize=(F.WIDTH_IN, 3.0))

    byfam = collections.defaultdict(list)
    for r in mass:
        byfam[r["family"]].append(r)
    fams = sorted(byfam, key=lambda f: -len(byfam[f]))[:8]
    y = np.arange(len(fams))[::-1]
    for yi, fam in zip(y, fams):
        rs = byfam[fam]
        modal = collections.Counter(r["core_seq"] for r in rs).most_common(1)[0][0]
        r = next(x for x in rs if x["core_seq"] == modal)
        lad = [float(v.split(":")[1]) for v in r["ladder_mono"].split(";")]
        axC.plot([min(lad), max(lad)], [yi, yi], color=F.FAINT, lw=1.0, zorder=1)
        axC.scatter(lad, [yi] * len(lad), s=10, color=F.KIND["RiPP"], zorder=3)
        if r["ladder_mono_decarb"]:
            lad2 = [float(v.split(":")[1]) for v in r["ladder_mono_decarb"].split(";")]
            axC.scatter(lad2, [yi - 0.22] * len(lad2), s=10, marker="s", color=F.NAVY,
                        zorder=3)
        if r["published"]:
            axC.scatter([float(r["published_mass"])], [yi + 0.26], s=34, marker="*",
                        color=F.RED, zorder=4)
            axC.text(float(r["published_mass"]) + 60, yi + 0.26,
                     f"{r['published']} ({float(r['best_ladder_gap']):.2f} Da)",
                     fontsize=5.6, va="center", color=F.RED)
    axC.set_yticks(y)
    axC.set_yticklabels([f"{f} (n={len(byfam[f])})" for f in fams], fontsize=6.4)
    axC.set_xlabel("predicted monoisotopic mass of the mature core (Da)", fontsize=7.5)
    handles = [plt.Line2D([], [], ls="", marker=m, color=c, label=l, markersize=5)
               for m, c, l in (("o", F.KIND["RiPP"], "dehydrations"),
                               ("s", F.NAVY, "+ C-term decarboxylation"),
                               ("*", F.RED, "published mass"))]
    # 圖例原本壓在 F004／F008 兩列上，移到圖外（2026-09-21 審閱意見）
    axC.legend(handles=handles, frameon=False, fontsize=6.0, loc="lower center",
               bbox_to_anchor=(0.5, 1.02), ncol=3, handlelength=1.1)
    F.bare(axC); F.panel(axC, "A", dx=-60)

    seen = {}
    for r in mass:
        seen.setdefault((r["family"], r["core_seq"]), []).append(r)
    for (fam, core), rs in seen.items():
        axD.scatter(len(core), int(rs[0]["n_ser_thr"]), s=12 + 6 * len(rs),
                    color=F.KIND["RiPP"], alpha=0.75, edgecolor="white", linewidth=0.4)
    axD.set_xlabel("predicted core peptide length (aa)", fontsize=7.5)
    axD.set_ylabel("Ser + Thr in the core", fontsize=7.5)
    F.bare(axD); F.panel(axD, "B", dx=-44)

    fig.tight_layout(w_pad=3.0)
    F.save(fig, "FigureS8_mass")


if __name__ == "__main__":
    main()
