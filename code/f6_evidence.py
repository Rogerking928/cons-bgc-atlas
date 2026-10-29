#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 6：偵測器的外部驗證與可解釋性（含一個失敗的因果檢定）。"""
import collections, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt
import numpy as np

import figstyle as F

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
GROUPS = [("profile-detected precursor", F.NAVY, "profile"),
          ("candidate, q < 0.05", F.RED, "cand. q<0.05"),
          ("candidate, q < 0.20", F.AMBER, "cand. q<0.20"),
          ("negative", F.GREY, "negative")]


def main():
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    sal = json.load(open(os.path.join(OUT, "finetune_summary.json")))
    mut = json.load(open(os.path.join(OUT, "mutagenesis.json")))
    mrows = list(csv.DictReader(open(os.path.join(OUT, "mutagenesis.csv"))))
    srows = list(csv.DictReader(open(os.path.join(OUT, "saliency.csv"))))
    gain = json.load(open(os.path.join(OUT, "atlas_gain_summary.json")))
    cal = json.load(open(os.path.join(OUT, "fdr_calibration.json")))
    en = rep["enrichment_near_machinery"]

    fig = plt.figure(figsize=(F.WIDTH_IN, 6.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 1.0, 1.0], hspace=0.85, wspace=0.46)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])
    axD = fig.add_subplot(gs[1, :])
    axE = fig.add_subplot(gs[2, 0])
    axG = fig.add_subplot(gs[2, 1:])

    # A 不靠標籤的外部檢定
    axA.bar([0, 1], [100 * en["high_near"] / (en["high_near"] + en["high_far"]),
                     100 * en["rest_near"] / (en["rest_near"] + en["rest_far"])],
            0.55, color=[F.NAVY, F.FAINT], edgecolor="white", linewidth=0.4)
    axA.set_xticks([0, 1]); axA.set_xlim(-0.55, 1.55)
    axA.set_xticklabels(["high\nscorers", "all other\nshort CDS"], fontsize=6.2)
    axA.set_ylabel("within 10 kb of cluster\nmachinery (%)", fontsize=7)
    np_ = en.get("n_permutations", 200)
    axA.set_title(f"OR {en['odds_ratio']:.0f}\npermutation P < {1 / np_:.4f}",
                  fontsize=6.4, loc="left")
    F.bare(axA); F.panel(axA, "A", dx=-40)

    # B 圖譜補洞（只有 RiPP 與 hybrid 家族有前驅肽可言）
    kind_of = {r["family"]: r["kind"] for r in csv.DictReader(
        open(os.path.join(OUT, "families.tsv")), delimiter="\t")}
    fams = [(f, d) for f, d in sorted(gain.items())
            if f.startswith("F") and d["clusters"] >= 4
            and kind_of.get(f) in ("RiPP", "hybrid")]
    x = np.arange(len(fams)); w = 0.27
    axB.bar(x - w, [100 * d["with_known_precursor"] / d["clusters"] for _, d in fams],
            w, color=F.NAVY, edgecolor="white", linewidth=0.4, label="profile")
    axB.bar(x, [100 * d["gaining_a_candidate_q05"] / d["clusters"] for _, d in fams],
            w, color=F.AMBER, edgecolor="white", linewidth=0.4, label="model, q<0.05")
    axB.bar(x + w, [100 * d["gaining_a_candidate"] / d["clusters"] for _, d in fams],
            w, color=F.AMBER, hatch="///", alpha=0.55, edgecolor="white",
            linewidth=0.4, label="model, q<0.20")
    axB.set_xticks(x); axB.set_xticklabels([f for f, _ in fams], fontsize=5.2)
    axB.set_ylim(0, 132); axB.set_yticks([0, 50, 100])
    axB.set_ylabel("clusters with a precursor (%)", fontsize=6.6)
    axB.set_xlabel("RiPP and hybrid families", fontsize=6.4)
    # 三欄排出去比格子還寬，會壓到隔壁 panel：改兩欄兩列，y 軸上緣讓位
    axB.legend(frameon=False, fontsize=5.0, ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, 1.06), handlelength=1.1, columnspacing=0.8,
               labelspacing=0.25)
    F.bare(axB); F.panel(axB, "B", dx=-38)

    # C 歸因統計（最高歸因殘基是否緊鄰 GG）＋ inset：對齊後的曲線
    near = sal["near_gg"]
    labs, vals, ns, cols = [], [], [], []
    for g, col, short in GROUPS:
        if g not in near:
            continue
        k, n = near[g]
        labs.append(short); vals.append(100 * k / n); ns.append(n); cols.append(col)
    y = np.arange(len(labs))[::-1]
    axC.barh(y, vals, color=cols, edgecolor="white", linewidth=0.4, height=0.66)
    for yi, v, n in zip(y, vals, ns):
        axC.text(v + 2, yi, f"{v:.0f}% ({n})", va="center", fontsize=5.2, color=F.INK)
    axC.set_yticks(y); axC.set_yticklabels(labs, fontsize=5.6)
    axC.set_xlim(0, 150); axC.set_xticks([0, 50, 100])
    axC.set_xlabel("top-attributed residue\nnext to a GG motif (%)", fontsize=6.4)
    axC.tick_params(axis="x", pad=3)
    ins = axC.inset_axes([0.66, 0.54, 0.33, 0.42])
    meta = sal.get("meta_profile", {})
    for g, col, _ in GROUPS:
        prof = meta.get(g)
        if not prof:
            continue
        xs = sorted(int(k) for k in prof)
        ins.plot(xs, [prof[str(v)][0] for v in xs], color=col, lw=0.8)
    ins.axvline(0, color=F.INK, lw=0.7, ls=":")
    ins.set_xticks([-30, 0, 30]); ins.set_yticks([0])
    ins.tick_params(labelsize=4.4, length=2, pad=1)
    ins.set_xlabel("residue relative to GG", fontsize=4.6, labelpad=1)
    ins.set_title("median attribution: every\ngroup peaks at the motif",
                  fontsize=4.6, pad=2, loc="left", color=F.GREY)
    F.bare(ins)
    F.bare(axC); F.panel(axC, "C", dx=-46)

    # D 擾動實驗：把 GG 改掉，分數有沒有掉
    groups = [g for g in ("candidate", "profile precursor", "negative") if g in mut]
    x = np.arange(len(groups))
    rng = np.random.default_rng(20260920)
    # 左兩組是「刪掉 motif」、右一組是「插入 motif」，不能放在同一條軸上直接互比
    axD.axvline(len(groups) - 1.5, color=F.FAINT, lw=0.8)
    axD.text(0.5, 1.02, "motif removed (GG \u2192 AA)", transform=axD.transAxes,
             ha="right", fontsize=6.0, color=F.GREY)
    axD.text(0.98, 1.02, "motif inserted", transform=axD.transAxes,
             ha="right", fontsize=6.0, color=F.GREY)
    for i, g in enumerate(groups):
        d = np.array([float(r["delta"]) for r in mrows if r["group"] == g])
        axD.scatter(rng.normal(i, 0.06, len(d)), d, s=7, alpha=0.6,
                    color=F.NAVY if g != "negative" else F.GREY, edgecolor="none")
        axD.plot([i - 0.24, i + 0.24], [np.median(d)] * 2, color=F.RED, lw=1.6)
        c = mut[g].get("control_delta_median")
        if c is None and g == "negative":
            # 陰性組的對照：同樣插入兩個殘基，但插的是 AA 而不是 GG
            c = mut[g].get("control_delta_median_insert", 0.0)
        if c is not None:
            axD.plot([i - 0.24, i + 0.24], [c] * 2, color=F.INK, lw=1.0, ls="--")
        axD.text(i, np.max(d) + 0.6,
                 f"{mut[g]['kind']}\nmedian {mut[g]['delta_median']:+.2f}",
                 ha="center", fontsize=5.6, color=F.INK)
    axD.axhline(0, color=F.GREY, lw=0.8, ls=":")
    axD.set_xticks(x)
    axD.set_xticklabels([f"{g}\n(n={mut[g]['n']})" for g in groups], fontsize=6.2)
    axD.set_xlim(-0.6, len(groups) - 0.4)
    axD.set_ylabel("change in detector log-odds\nafter editing the motif", fontsize=7)
    axD.plot([], [], color=F.RED, lw=1.6, label="median change")
    axD.plot([], [], color=F.INK, lw=1.0, ls="--", label="control: same edit elsewhere")
    # lower right 正好落在 negative 那團散點裡，字疊在點上；軸下方又會撞到
    # 兩行的刻度標籤。lower center 落在 profile precursor 那一組下方的空白處。
    axD.legend(frameon=False, fontsize=5.6, loc="lower center", ncol=2,
               handlelength=1.6, columnspacing=1.6)
    F.bare(axD); F.panel(axD, "D", dx=-44)

    # E, F 誘餌校準（原本在補充；它是「模型到底在讀什麼」的下半句，放這裡才連得上）
    cur = cal["curve"]
    t = [c["threshold"] for c in cur]
    series = [("real_rate", F.NAVY, "-", f"real short CDS (n={cal['n_real_distinct']:,})"),
              ("decoy_background_rate", F.RED, "--",
               f"shuffled background (n={cal['n_decoy_background']:,})"),
              ("decoy_from_positives_rate", F.AMBER, "-.",
               f"shuffled precursors (n={cal['n_decoy_from_positives']})")]
    for key, col, ls, lab in series:
        axE.plot(t, [100 * c[key] for c in cur], color=col, ls=ls, lw=1.2, label=lab)
    for name, c in cal["thresholds"].items():
        axE.axvline(c["score"], color=F.GREY, lw=0.8, ls=":")
        axE.text(c["score"], 400, name, rotation=0, fontsize=5.0,
                 color=F.GREY, ha="center", va="bottom")
    axE.set_yscale("log"); axE.set_ylim(0.01, 1200); axE.set_xlim(0.3, 1.0)
    axE.set_xticks([0.4, 0.6, 0.8, 1.0])
    axE.set_yticks([0.1, 10])
    axE.set_xlabel("detector score threshold", fontsize=7)
    axE.set_ylabel("above the calling threshold (%)", fontsize=7)
    axE.legend(frameon=False, fontsize=4.6, loc="upper center",
               bbox_to_anchor=(0.5, -0.34), ncol=1, handlelength=1.3,
               labelspacing=0.25)
    F.bare(axE); F.panel(axE, "E", dx=-40)

    names_ = list(cal["thresholds"])
    x = np.arange(len(names_)); w = 0.22
    bars = [("real_rate", F.NAVY, "o", "real"),
            ("decoy_background_rate", F.RED, "s", "shuffled background"),
            ("decoy_from_positives_rate", F.AMBER, "D", "shuffled precursors")]
    # 對數軸上的長條長度不再正比於數值，會誤導；改成點加細線
    for i, (key, col, mark, lab) in enumerate(bars):
        v = [100 * cal["thresholds"][n][key] for n in names_]
        xs = x + (i - 1) * w
        for xi, vv in zip(xs, v):
            axG.plot([xi, xi], [0.06, vv], color=F.FAINT, lw=0.8, zorder=1)
            axG.text(xi, vv * 1.5, f"{vv:.2f}" if vv < 10 else f"{vv:.0f}",
                     ha="center", fontsize=5.2, color=F.INK)
        axG.scatter(xs, v, s=22, color=col, marker=mark, zorder=3, label=lab)
    axG.set_yscale("log"); axG.set_ylim(0.06, 3000)
    axG.set_xlim(-0.5, len(names_) - 0.5)
    axG.set_xticks(x); axG.set_xticklabels(names_, fontsize=6.4)
    axG.set_ylabel("above the calling threshold (%)", fontsize=7)
    axG.set_xlabel("shuffling keeps composition and destroys order: the shuffled "
                   "precursors still pass", fontsize=5.8)
    axG.legend(frameon=False, fontsize=5.0, loc="upper center",
               bbox_to_anchor=(0.5, 1.04), ncol=3, handlelength=0.9, columnspacing=0.7)
    F.bare(axG); F.panel(axG, "F", dx=-40)

    F.save(fig, "Figure6_evidence")


if __name__ == "__main__":
    main()
