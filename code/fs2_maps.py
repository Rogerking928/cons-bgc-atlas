#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure S2：其餘代表基因簇的座位圖。主圖 Figure 2A/B 已經放了質體與染色體
各一張；這裡四張共用同一個 kb 尺度，短的就讓它短，不要拉滿寬度。"""
import collections, csv, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch

import figstyle as F
from names import FAMILY

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")

SHOW = ["F001", "F004", "F005", "F011"]
CLASS = {"structural": ("Gallidermin", "L_biotic_typeA", "Lanthi_prec",
                        "Mersacidin", "Antimicrobial17", "Antimicrobial18",
                        "Bacteriocin_II", "Bacteriocin_IIc", "Bacteriocin_IId",
                        "MccV", "Lactococcin_972"),
         "maturation": ("LanM_cyc", "LanC", "YcaO", "LanB_N", "LanB_C",
                        "Peptidase_S8", "A", "C", "KS", "PCP/ACP", "TE", "MT",
                        "KR", "A_C", "ER/NAD", "Radical_SAM"),
         "export": ("C39",),
         "mobility": ("Transposase", "Integrase", "Recombinase")}
# 一定要用字界：舊版寫 "erm" 會把 "permease" 也塗成抗藥基因（2026-09-20 踩過）
AMR_RE = re.compile(r"\bbeta-lactamase|\bBlaZ\b|\bBlaR1?\b|\bBlaI\b|"
                    r"\bLnu\(|\bErm\(|\berm\(|rRNA adenine N-6-methyltransferase|"
                    r"\bmecA\b|PBP2a|methicillin|"
                    r"aminoglycoside (?:O-)?(?:phospho|nucleotidyl|acetyl)transferase|"
                    r"\btet\(|tetracycline efflux|macrolide efflux|"
                    r"streptogramin|quaternary ammonium compound|"
                    r"cadmium-translocating|arsenate reductase|arsenite efflux|"
                    r"chromate efflux", re.I)
COLOR = {"structural": F.ROLE["structural"], "maturation": F.ROLE["maturation"],
         "export": F.ROLE["export"], "immunity": F.ROLE["immunity"],
         "regulation": F.ROLE["regulation"], "mobility": F.ROLE["mobility"],
         "resistance": F.ROLE["resistance"], "other": F.ROLE["other"]}
# 基因名縮寫與標籤規則只有一份，兩張座位圖共用（重複過一次就漂移過一次）
from f3_supp_maps import SHORT, klass, short  # noqa: E402


def main():
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = {r["cluster"]: r for r in csv.DictReader(
        open(os.path.join(OUT, "clusters.tsv")), delimiter="\t")}
    cds = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"):
        cds[r["cluster"]].append(r)

    # 先算出最長的一張，四張共用同一個 kb-per-inch
    longest = 0
    for fam in SHOW:
        for c, f in fam_of.items():
            if f != fam:
                continue
            g = cds[c]
            longest = max(longest, int(g[-1]["end"]) - int(g[0]["start"]))
    fig, axes = plt.subplots(len(SHOW), 1, figsize=(F.WIDTH_IN, 7.4))
    for ax, fam, letter in zip(axes, SHOW, "ABCDEF"):
        cs = sorted([c for c, f in fam_of.items() if f == fam],
                    key=lambda c: int(cl[c]["n_cds"]))
        # 這張圖要讓讀者看到抗藥基因就在旁邊，所以家族裡若有帶抗藥基因的座，
        # 就挑帶最多的那個當代表（沒有的話回到中位大小的那個）
        def n_amr(c):
            return sum(1 for g in cds[c] if AMR_RE.search(g["product"] or ""))
        rep = max(cs, key=n_amr) if any(n_amr(c) for c in cs) else cs[len(cs) // 2]
        genes = sorted(cds[rep], key=lambda r: int(r["start"]))
        lo = int(genes[0]["start"])
        span = int(genes[-1]["end"]) - lo
        row = 0
        last = ("", 0)
        for g in genes:
            s, e = int(g["start"]) - lo, int(g["end"]) - lo
            k = klass((g["domains"] or "").split("-"), g["product"])
            x, dx = (s, e - s) if g["strand"] == "+" else (e, s - e)
            ax.add_patch(FancyArrow(
                x, 0, dx, 0, width=0.34, head_width=0.52,
                length_includes_head=True,
                head_length=min(abs(dx) * 0.35, span * 0.018),
                facecolor=COLOR[k], edgecolor="white", linewidth=0.5))
            lab = short(g["product"], k not in ("other", "export"))
            if lab and lab == last[0] and (s - last[1]) < span * 0.12:
                # 連著好幾個一模一樣的基因（class III lanthipeptide 前驅肽有四條）
                # 只標一次，否則字會疊在一起
                lab = ""
            elif lab:
                last = (lab, e)
            if lab:
                ax.text((s + e) / 2, 0.42 + (row % 4) * 0.90, lab, ha="center",
                        va="bottom", fontsize=6.0, color=F.INK)
                row += 1
        a = cl[rep]
        ax.set_xlim(-longest * 0.02, longest * 1.02)   # 共用尺度
        ax.set_ylim(-0.6, 4.6)
        ax.set_yticks([])
        ticks = [t for t in (0, 10000, 20000, 30000, 40000) if t <= longest]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t // 1000}" for t in ticks], fontsize=7)
        ax.set_xlabel("kb", fontsize=7)
        ax.set_title(f"{fam}  {FAMILY[fam][0]} — {a['assembly']}, "
                     f"{a['replicon_type']}, core {int(a['span']) / 1000:.1f} kb",
                     loc="left", fontsize=8, pad=26)
        F.bare(ax, keep=("bottom",))
        F.panel(ax, letter, dx=-22, dy=2)
    handles = [Patch(facecolor=COLOR[k], label=k) for k in
               ("structural", "maturation", "export", "regulation",
                "mobility", "resistance", "other")]
    axes[-1].legend(handles=handles, frameon=False, fontsize=7, ncol=7,
                    loc="upper center", bbox_to_anchor=(0.5, -0.5))
    fig.tight_layout(h_pad=3.6)
    F.save(fig, "FigureS2_maps")


if __name__ == "__main__":
    main()
