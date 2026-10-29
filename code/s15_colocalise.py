#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Physical linkage: which resistance determinants sit on the same replicon as
each cluster, and how far away.

Species-intrinsic determinants are excluded from the linkage tables: lsa(A) is
in essentially every E. faecalis chromosome and msr(C)/aac(6')_Entco in every
E. faecium chromosome, so their co-residence carries no information. They are
kept aside as a positive control for the detection method.
"""
import collections, csv, json, os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
# 葡萄球菌版：fusD 是 S. saprophyticus 的內生性 fusidic acid 抗性，
# mecA-ceftar 與 mecA、qacB 與 qacA 在這個資料庫裡是同一個基因的兩筆參考，
# 一起算會把同一件事數兩次。
INTRINSIC = {"fusD", "mecA-ceftar", "qacB"}


def main():
    loc = list(csv.DictReader(open(os.path.join(OUT, "amr_located.csv"))))
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    fam = json.load(open(os.path.join(OUT, "cluster_family.json")))
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    by = collections.defaultdict(list)
    for r in loc:
        if r["determinant"] not in INTRINSIC:
            by[(r["assembly"], r["replicon"])].append(r)

    rows, tally = [], collections.defaultdict(collections.Counter)
    for c in cl:
        f = fam[c["cluster"]]
        near = by.get((c["assembly"], c["replicon"]), [])
        lo, hi = int(c["start"]), int(c["end"])
        d = None
        for r in near:
            s, e = int(r["start"] or 0), int(r["end"] or 0)
            gap = 0 if (s <= hi and e >= lo) else min(abs(s - hi), abs(lo - e))
            d = gap if d is None else min(d, gap)
        rows.append(dict(
            cluster=c["cluster"], family=f, assembly=c["assembly"],
            species=mat[c["assembly"]]["species"], replicon=c["replicon"],
            replicon_type=c["replicon_type"], plasmid=c["plasmid"],
            replicon_len=c["replicon_len"], cluster_span=c["span"],
            n_acquired_on_replicon=len({r["determinant"] for r in near}),
            determinants=";".join(sorted({r["determinant"] for r in near})),
            nearest_bp="" if d is None else d))
        for r in {r["determinant"] for r in near}:
            tally[(f, c["replicon_type"])][r] += 1
    p = os.path.join(OUT, "colocalisation.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("wrote", p, len(rows), "clusters")
    for key in sorted(tally):
        n = sum(1 for r in rows if (r["family"], r["replicon_type"]) == key)
        withany = sum(1 for r in rows if (r["family"], r["replicon_type"]) == key
                      and r["n_acquired_on_replicon"])
        print(f"  {key[0]} on {key[1]:10s}: {withany}/{n} share the replicon with an "
              f"acquired determinant; top {tally[key].most_common(5)}")


if __name__ == "__main__":
    main()
