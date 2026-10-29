#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MIBiG 命中併成一張表，並接回圖譜的家族。

s56 每株寫一個 .kcb.tsv；這支把它們併起來，對每個 region 補上
「哪一個圖譜的簇與它重疊、那個簇屬於哪個家族」，讓審稿人 1 問的
「已知 BGC 怎麼判定、是不是 MIBiG」可以逐列看。

本圖譜的家族名來自 RefSeq product name 與 UniProt（§2.3），**不是** MIBiG，
所以這張表是對同一批座位的**獨立**指認，不是命名的來源。

Output: out/mibig_hits.csv、out/mibig_summary.json
"""
import collections, csv, glob, json, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SRC = os.path.join(OUT, "antismash_kcb")


def main():
    files = sorted(glob.glob(os.path.join(SRC, "*.kcb.tsv")))
    if not files:
        sys.exit("no knownclusterblast output - run s56 first")

    # region -> 圖譜的簇與家族（s55 已經把 region 對上簇了）
    # region 編號每條 replicon 各自從 1 開始，所以鍵一定要含 replicon，
    # 否則染色體的 region 1 會被質體的 region 1 蓋掉。
    fam_of = {}
    for r in csv.DictReader(open(os.path.join(OUT, "antismash_compare.csv"))):
        if r["as_region"]:
            fam_of[(r["assembly"], r["replicon"], r["as_region"])] = (
                r["cluster"], r["family"])
    scope = {}
    for r in csv.DictReader(open(os.path.join(OUT, "antismash_regions.csv"))):
        scope[(r["assembly"], r["replicon"], r["region_number"])] = (
            r["scope"], r["product"], r["replicon"])

    rows, n_hit = [], 0
    for path in files:
        for r in csv.DictReader(open(path), delimiter="\t"):
            key = (r["assembly"], r["replicon"], r["region_number"])
            sc, prod, rep = scope.get(key, ("", "", ""))
            cluster, family = fam_of.get(key, ("", ""))
            if r["mibig_accession"]:
                n_hit += 1
            rows.append([r["assembly"], rep, r["region_number"], prod, sc,
                         cluster, family, r["mibig_accession"],
                         r["mibig_description"], r["mibig_type"],
                         r["similarity"], r["core_gene_hits"], r["hits"],
                         r["mibig_version"]])

    rows.sort(key=lambda x: (x[0], x[1], int(x[2] or 0)))
    with open(os.path.join(OUT, "mibig_hits.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["assembly", "replicon", "region_number", "antismash_product",
                    "scope", "atlas_cluster", "atlas_family", "mibig_accession",
                    "mibig_description", "mibig_type", "similarity",
                    "core_gene_hits", "blast_hits", "mibig_version"])
        w.writerows(rows)

    # 每個家族最常對到哪一個 MIBiG 條目
    by_fam = collections.defaultdict(collections.Counter)
    for r in rows:
        if r[6] and r[7]:
            by_fam[r[6]][f"{r[7]} {r[8]}"] += 1
    best = {}
    for fam, c in sorted(by_fam.items()):
        top, n = c.most_common(1)[0]
        sims = [float(r[10]) for r in rows
                if r[6] == fam and f"{r[7]} {r[8]}" == top and r[10]]
        best[fam] = {"top_mibig": top, "clusters_with_that_hit": n,
                     "clusters_with_any_hit": sum(c.values()),
                     "median_similarity": (sorted(sims)[len(sims) // 2]
                                           if sims else None)}
    summary = {"regions": len(rows), "regions_with_a_mibig_hit": n_hit,
               "genomes": len(files), "per_family": best}
    with open(os.path.join(OUT, "mibig_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    print(f"{len(rows)} regions in {len(files)} genomes; "
          f"{n_hit} have a MIBiG hit")
    for fam, d in best.items():
        print(f"  {fam}: {d['top_mibig']} in {d['clusters_with_that_hit']}"
              f"/{d['clusters_with_any_hit']}, median similarity "
              f"{d['median_similarity']}")
    print("-> out/mibig_hits.csv, out/mibig_summary.json")


if __name__ == "__main__":
    main()
