#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 3：四個代表基因簇的座位圖（含兩側 5 kb 的鄰居，抗藥基因塗黑）。"""
import collections, csv, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch

import figstyle as F
from names import FAMILY

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")

SHOW = ["F005", "F004", "F009", "F011"]
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
SHORT = [("lugdunin non-ribosomal peptide synthetase LugA", "lugA"),
         ("lugdunin non-ribosomal peptide synthetase LugB", "lugB"),
         ("lugdunin non-ribosomal peptide synthetase LugC", "lugC"),
         ("lugdunin non-ribosomal peptide synthetase LugD", "lugD"),
         ("thioesterase II family protein", "TE II"),
         ("aureusimine non-ribosomal peptide synthetase AusA", "ausA"),
         ("aureusimine biosynthesis 4'-phosphopantetheinyl transferase AusB", "ausB"),
         ("gallidermin family lantibiotic", "lanA"),
         ("lantibiotic dehydratase", "lanB"),
         ("lanthionine synthetase C family protein", "lanC"),
         ("type 2 lanthipeptide synthetase LanM", "lanM"),
         ("class III lanthionine synthetase LanKC", "lanKC"),
         ("class III lanthipeptide", "lanA (class III)"),
         ("lacticin 481 family lantibiotic", "lanA (lacticin 481)"),
         ("peptidase domain-containing ABC transporter", "ABC (C39)"),
         ("penicillin-hydrolyzing class A beta-lactamase BlaZ", "blaZ"),
         ("beta-lactam sensor/signal transducer BlaR1", "blaR1"),
         ("penicillinase repressor BlaI", "blaI"),
         ("lincosamide nucleotidyltransferase Lnu(A)", "lnu(A)"),
         ("IS6-like element IS257 family transposase", "IS257"),
         ("IS6 family transposase", "IS6"), ("IS3 family transposase", "IS3"),
         ("IS110-like element", "IS110"), ("ISL3 family transposase", "ISL3"),
         ("bifunctional aminoglycoside", "aac(6')-aph(2'')"),
         ("aminoglycoside O-phosphotransferase", "aph"),
         ("aminoglycoside nucleotidyltransferase", "ant"),
         ("aminoglycoside O-nucleotidyltransferase ANT(4')", "ant(4')-Ia"),
         ("Mu transposase", "Mu transposase"),
         ("YcaO-like family protein", "YcaO"),
         ("YcaO domain", "YcaO"),
         ("IS256-like element", "IS256"),
         ("IS1182-like element", "IS1182"),
         ("thiopeptide-type bacteriocin biosynthesis protein", "thiopeptide"),
         ("TnsB C-terminal", "TnsB"),
         ("cadmium-translocating P-type ATPase", "cadA"),
         ("arsenate reductase", "arsC"), ("chromate efflux transporter", "chrA"),
         ("IS1182 family transposase", "IS1182"),
         ("recombinase family protein", "recombinase"),
         ("sensor histidine kinase", "HK"),
         ("response regulator transcription factor", "RR"),
         ("ABC transporter ATP-binding protein", "ABC"),
         ("ABC transporter permease", "ABC perm."),
         ("hypothetical protein", "")]


def klass(doms, product):
    if AMR_RE.search(product or ""):
        return "resistance"
    for k, v in CLASS.items():
        if any(d in v for d in doms):
            return k
    p = (product or "").lower()
    if "histidine kinase" in p or "response regulator" in p or "repressor" in p:
        return "regulation"
    if "abc transporter" in p or "permease" in p or "efflux" in p:
        return "export"
    return "other"


def short(p, labelled):
    for long, s in SHORT:
        if long.lower() in (p or "").lower():
            return s
    if not labelled:
        return ""
    # 兩行標籤在密集的座位圖上互相壓；常見的用上面的短名，其餘單行截到 14 字
    p = (p or "").split(" domain-containing")[0]
    return p if len(p) <= 14 else p[:13] + "\u2026"


def main():
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = {r["cluster"]: r for r in csv.DictReader(
        open(os.path.join(OUT, "clusters.tsv")), delimiter="\t")}
    cds = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"):
        cds[r["cluster"]].append(r)

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
        ax.set_xlim(-span * 0.02, span * 1.02)
        ax.set_ylim(-0.6, 4.6)
        ax.set_yticks([])
        ax.set_xticks([0, span // 2, span])
        ax.set_xticklabels(["0", f"{span // 2 / 1000:.1f}", f"{span / 1000:.1f} kb"],
                           fontsize=7)
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
