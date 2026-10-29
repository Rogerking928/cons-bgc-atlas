#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GC content of every called cluster against its own genome.

Writes the windowed GC profile of a few representative replicons for the figure,
and the cluster-minus-genome difference for every cluster in the panel.
"""
import collections, csv, gzip, json, os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FNA = os.path.join(HERE, "data", "fna")
WIN = 2000
SHOW = {"F001": None, "F004": None, "F006": None, "F009": None}


def load_fna(acc):
    seqs, n, b = {}, None, []
    with gzip.open(os.path.join(FNA, acc + ".fna.gz"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    seqs[n] = "".join(b)
                n, b = line[1:].split()[0], []
            else:
                b.append(line.strip().upper())
    if n:
        seqs[n] = "".join(b)
    return seqs


def gc(s):
    n = len(s)
    return 100.0 * (s.count("G") + s.count("C")) / n if n else float("nan")


def main():
    fam = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    by_asm = collections.defaultdict(list)
    for c in cl:
        by_asm[c["assembly"]].append(c)

    # one representative cluster per family to draw. Where another item of the
    # document already shows a particular assembly, the same one is used here,
    # so the reader is looking at one genome and not three.
    PREFER = {"F002": "GCF_030020495.1", "F001": "GCF_053154045.1",
              "F004": "GCF_052084685.1"}
    for f in SHOW:
        cand = [c for c in cl if fam[c["cluster"]] == f]
        if not cand:
            continue
        pick = [c for c in cand if c["assembly"] == PREFER.get(f)]
        if pick:
            SHOW[f] = pick[0]["cluster"]
        else:
            cand.sort(key=lambda c: int(c["span"]))
            SHOW[f] = cand[len(cand) // 2]["cluster"]

    rows, profiles = [], {}
    for n, (a, cs) in enumerate(sorted(by_asm.items()), 1):
        seqs = load_fna(a)
        gen = "".join(seqs.values())
        gg = gc(gen)
        for c in cs:
            s = seqs.get(c["replicon"], "")
            if not s:
                continue
            sub = s[int(c["start"]) - 1:int(c["end"])]
            rows.append(dict(cluster=c["cluster"], family=fam[c["cluster"]],
                             assembly=a, replicon_type=c["replicon_type"],
                             span=c["span"], cluster_gc=round(gc(sub), 2),
                             genome_gc=round(gg, 2),
                             delta=round(gc(sub) - gg, 2)))
            if c["cluster"] in SHOW.values():
                prof = [gc(s[i:i + WIN]) for i in range(0, len(s) - WIN, WIN)]
                profiles[c["cluster"]] = dict(
                    family=fam[c["cluster"]], assembly=a, replicon=c["replicon"],
                    replicon_type=c["replicon_type"], window=WIN,
                    genome_gc=round(gg, 2), replicon_len=len(s),
                    start=int(c["start"]), end=int(c["end"]),
                    gc=[round(x, 2) for x in prof])
        if n % 100 == 0:
            print(f"  {n}/{len(by_asm)} genomes", flush=True)

    with open(os.path.join(OUT, "cluster_gc.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    json.dump(profiles, open(os.path.join(OUT, "gc_profiles.json"), "w"))
    for f in sorted(SHOW):
        d = [r["delta"] for r in rows if r["family"] == f]
        if d:
            print(f"  {f}: median cluster-minus-genome GC {np.median(d):+.2f} points "
                  f"(n={len(d)})")
    print("-> out/cluster_gc.csv, out/gc_profiles.json")


if __name__ == "__main__":
    main()
