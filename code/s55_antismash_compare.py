#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The atlas against the antiSMASH regions of the same genomes.

審稿人（AoM #42 Reviewer 1，2026-10-28）：「It is unclear why antiSMASH, which is
currently the gold standard, is not compared against and only DeepBGC is.」

比較只有在**範圍對齊**時才誠實。這份圖譜的 45 個 Pfam 只涵蓋非核醣體胜肽合成酶、
模組化 PKS，與 RiPP（羊毛硫肽、thiopeptide／azole、sactipeptide、細菌素）；
antiSMASH 的規則還會呼叫 terpene、第三型 PKS、NIS 型 siderophore、staphylopine
（opine-like-metallophore）、agr 自體誘導物與 CDPS —— 這些**本圖譜設計上就不收**，
而且它們是 antiSMASH 呼叫數裡的大多數（agr 與 terpene-precursor 在 554 株每株都有）。
所以：

- 範圍內（IN_SCOPE）：兩者都該找到的類別，用來算互相的召回率與邊界差異。
- 範圍外：只統計 antiSMASH 找到什麼，說明本圖譜不是一份通用 BGC 普查。

一個 region 可以掛多個 product（例如 'NRPS;azole-containing-RiPP'），只要有一個
落在範圍內就算範圍內。

Output: out/antismash_compare.csv（每個圖譜簇一列）
        out/antismash_regions.csv（每個 antiSMASH region 一列）
        out/antismash_compare.json（摘要）
"""
import collections, csv, glob, json, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SRC = os.path.join(OUT, "antismash")

# antiSMASH 的 product 名稱 -> 是否落在本圖譜的 Pfam 範圍內
IN_SCOPE = {
    "NRPS", "NRPS-like", "transAT-PKS", "transAT-PKS-like", "T1PKS", "hglE-KS",
    "lanthipeptide-class-i", "lanthipeptide-class-ii", "lanthipeptide-class-iii",
    "lanthipeptide-class-iv", "lanthipeptide-class-v", "lipolanthine",
    "thiopeptide", "azole-containing-RiPP", "RaS-RiPP", "RiPP-like",
    "sactipeptide", "thioamitides", "epipeptide", "lassopeptide",
    "linaridin", "cyanobactin", "glycocin", "bacteriocin", "head_to_tail",
    # RRE = RiPP recognition element，是 RiPP 的指標，算範圍內
    "RRE-containing",
}
# 明確記錄為範圍外的（出現在本 panel 的），其餘未知名稱也歸範圍外並列出來
OUT_OF_SCOPE = {
    "terpene", "terpene-precursor", "T3PKS", "NI-siderophore",
    "opine-like-metallophore", "cyclic-lactone-autoinducer", "CDPS",
    "betalactone", "ectoine", "phosphonate", "resorcinol", "butyrolactone",
}


def products(cell):
    return [p for p in (cell or "").split(";") if p]


def scope_of(cell):
    ps = products(cell)
    return "in" if any(p in IN_SCOPE for p in ps) else "out"


def load_regions():
    regs = collections.defaultdict(list)       # (assembly, replicon) -> [row]
    unknown = collections.Counter()
    n = 0
    for path in sorted(glob.glob(os.path.join(SRC, "*.regions.tsv"))):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r["start"], r["end"] = int(r["start"]), int(r["end"])
                r["scope"] = scope_of(r["product"])
                for p in products(r["product"]):
                    if p not in IN_SCOPE and p not in OUT_OF_SCOPE:
                        unknown[p] += 1
                regs[(r["assembly"], r["replicon"])].append(r)
                n += 1
    return regs, n, unknown


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def main():
    regs, n_reg, unknown = load_regions()
    if not n_reg:
        sys.exit("no antiSMASH output in out/antismash/ - run s54 first")
    if unknown:
        print("product names not in either scope list (treated as out of "
              f"scope): {dict(unknown)}", file=sys.stderr)

    fam = json.load(open(os.path.join(OUT, "cluster_family.json")))
    clusters = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")),
                                   delimiter="\t"))

    rows = []
    matched_region_ids = set()
    for c in clusters:
        key = (c["assembly"], c["replicon"])
        cs, ce = int(c["start"]), int(c["end"])
        ws, we = int(c["window_start"]), int(c["window_end"])
        best = None
        for r in regs.get(key, []):
            ov = overlap(cs, ce, r["start"], r["end"])
            ovw = overlap(ws, we, r["start"], r["end"])
            if ov == 0 and ovw == 0:
                continue
            cand = (ov, ovw, r)
            if best is None or (cand[0], cand[1]) > (best[0], best[1]):
                best = cand
        if best:
            ov, ovw, r = best
            matched_region_ids.add((r["assembly"], r["replicon"],
                                    r["region_number"]))
            rows.append([c["cluster"], c["assembly"], c["replicon"],
                         fam.get(c["cluster"], ""),
                         c["kind"], c["replicon_type"], cs, ce, ce - cs + 1,
                         r["region_number"], r["start"], r["end"],
                         r["end"] - r["start"] + 1, r["product"], r["scope"],
                         ov, ovw])
        else:
            rows.append([c["cluster"], c["assembly"], c["replicon"],
                         fam.get(c["cluster"], ""),
                         c["kind"], c["replicon_type"], cs, ce, ce - cs + 1,
                         "", "", "", "", "", "none", 0, 0])

    with open(os.path.join(OUT, "antismash_compare.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster", "assembly", "replicon", "family", "kind",
                    "replicon_type", "cluster_start", "cluster_end",
                    "cluster_span", "as_region", "as_start", "as_end",
                    "as_span", "as_product", "as_scope", "overlap_bp_span",
                    "overlap_bp_window"])
        w.writerows(rows)

    # 每個 antiSMASH region 一列，標記有沒有被圖譜的簇對上
    reg_rows = []
    for key in sorted(regs):
        for r in regs[key]:
            rid = (r["assembly"], r["replicon"], r["region_number"])
            reg_rows.append([r["assembly"], r["replicon"], r["region_number"],
                             r["start"], r["end"], r["end"] - r["start"] + 1,
                             r["product"], r["scope"], r["contig_edge"],
                             int(rid in matched_region_ids)])
    with open(os.path.join(OUT, "antismash_regions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["assembly", "replicon", "region_number", "start", "end",
                    "span", "product", "scope", "contig_edge",
                    "matched_by_atlas"])
        w.writerows(reg_rows)

    # ------------------------------------------------------------- 摘要
    # 具名索引：compare 表（C_）與 region 表（R_）的欄位位置
    C_FAMILY, C_SPAN, C_AS_SPAN, C_SCOPE = 3, 8, 12, 14
    R_PRODUCT, R_SCOPE, R_MATCHED = 6, 7, 9

    def med(xs):
        xs = sorted(xs)
        return None if not xs else (xs[len(xs) // 2] if len(xs) % 2
                                    else (xs[len(xs) // 2 - 1] + xs[len(xs) // 2]) / 2)

    hit = [r for r in rows if r[C_SCOPE] != "none"]
    hit_in = [r for r in rows if r[C_SCOPE] == "in"]
    in_regions = [r for r in reg_rows if r[R_SCOPE] == "in"]
    out_regions = [r for r in reg_rows if r[R_SCOPE] == "out"]
    by_fam = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        by_fam[r[C_FAMILY]][1] += 1
        if r[C_SCOPE] == "in":
            by_fam[r[C_FAMILY]][0] += 1
    prod_unmatched = collections.Counter(r[R_PRODUCT] for r in in_regions
                                         if not r[R_MATCHED])
    summary = {
        "atlas_clusters": len(rows),
        "antismash_regions": n_reg,
        "antismash_regions_in_scope": len(in_regions),
        "antismash_regions_out_of_scope": len(out_regions),
        "atlas_clusters_overlapped_by_any_region": len(hit),
        "atlas_clusters_overlapped_by_in_scope_region": len(hit_in),
        "in_scope_regions_matched_by_atlas": sum(r[R_MATCHED]
                                                 for r in in_regions),
        "in_scope_regions_not_matched_by_atlas": (
            len(in_regions) - sum(r[R_MATCHED] for r in in_regions)),
        "median_atlas_span_bp": med([r[C_SPAN] for r in rows]),
        "median_matched_region_span_bp": med([r[C_AS_SPAN] for r in hit_in]),
        "span_ratio_region_over_cluster": (
            round(med([r[C_AS_SPAN] for r in hit_in])
                  / med([r[C_SPAN] for r in hit_in]), 2)
            if hit_in else None),
        "per_family_recall_in_scope": {k: v for k, v in sorted(by_fam.items())},
        "out_of_scope_products": dict(collections.Counter(
            r[R_PRODUCT] for r in out_regions).most_common()),
        "in_scope_products_not_matched_by_atlas": dict(prod_unmatched.most_common()),
    }
    with open(os.path.join(OUT, "antismash_compare.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    print(f"atlas clusters {len(rows)}; antiSMASH regions {n_reg} "
          f"({len(in_regions)} in scope, {len(out_regions)} out)")
    print(f"  atlas clusters overlapped by an in-scope region: {len(hit_in)}"
          f" ({100*len(hit_in)/len(rows):.1f}%)")
    print(f"  in-scope regions with no atlas cluster: "
          f"{summary['in_scope_regions_not_matched_by_atlas']}")
    print(f"  median span: atlas {summary['median_atlas_span_bp']} bp vs "
          f"matched region {summary['median_matched_region_span_bp']} bp "
          f"(x{summary['span_ratio_region_over_cluster']})")
    print("-> out/antismash_compare.csv, out/antismash_regions.csv, "
          "out/antismash_compare.json")


if __name__ == "__main__":
    main()
