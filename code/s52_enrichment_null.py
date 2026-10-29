#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The positional null of Figure 6A, at a resolution that can report a small P.

s19 ran this null 200 times, which floors the P value at 1/201 = 0.005; the
number in the figure was therefore the smallest the test could produce rather
than a measurement. This script recomputes the same null 10,000 times from the
saved scan and writes the result back into precursor_report.json.

It rebuilds the two groups and the machinery positions exactly as s19 does and
refuses to run if the rebuilt 2x2 table differs from the stored one, so a
divergence in either definition is a failure rather than a silent difference.

  PYMEM_GB=8 python3 code/s52_enrichment_null.py
"""
import collections
import csv
import gzip
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
ANCHOR_DOMS = {"C39", "LanM_cyc", "LanC", "YcaO", "LanB_N", "LanB_C", "Peptidase_S8"}
NEAR = 10000
N_PERM = 10000
SEED = 20260918


def main():
    report = json.load(open(os.path.join(OUT, "precursor_report.json")))
    stored = report["enrichment_near_machinery"]
    thr = report["threshold"]["score"]

    trained = {r["seq"] for r in csv.DictReader(
        open(os.path.join(OUT, "precursor_trainset.csv")))}
    rows = [r for r in csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv")))
            if not r["known_domain"] and r["seq"] not in trained]

    def near(r):
        return r["nearest_anchor_bp"] != "" and int(r["nearest_anchor_bp"]) <= NEAR

    hi = [r for r in rows if float(r["score"]) >= thr]
    bg = [r for r in rows if float(r["score"]) < thr]
    a_ = sum(1 for r in hi if near(r))
    b_ = len(hi) - a_
    c_ = sum(1 for r in bg if near(r))
    d_ = len(bg) - c_
    if (a_, b_, c_, d_) != (stored["high_near"], stored["high_far"],
                            stored["rest_near"], stored["rest_far"]):
        sys.exit(f"FAIL: rebuilt table {(a_, b_, c_, d_)} differs from the stored "
                 f"{(stored['high_near'], stored['high_far'], stored['rest_near'], stored['rest_far'])}; "
                 f"the reconstruction is not faithful, so the null would not be the same test")

    dom = collections.defaultdict(dict)
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            dom[f[0]].setdefault(f[1], set()).add(f[3])
    anchors = collections.defaultdict(list)
    rep_len = {}
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, rep, prot, s, e = f[0], f[1], f[5], int(f[6]), int(f[7])
            rep_len[(a, rep)] = int(f[3])
            if dom.get(a, {}).get(prot, set()) & ANCHOR_DOMS:
                anchors[(a, rep)].append((s, e, prot))

    hi_by_rep = collections.defaultdict(list)
    for r in hi:
        hi_by_rep[(r["assembly"], r["replicon"])].append((int(r["start"]), int(r["end"])))

    rng = np.random.default_rng(SEED)
    null = np.zeros(N_PERM, dtype=int)
    work = []
    for key, spans in hi_by_rep.items():
        L = rep_len.get(key, 0)
        anc = anchors.get(key, [])
        if anc and L > 0:
            work.append((L, len(anc), np.array(spans)))
    for i in range(N_PERM):
        k = 0
        for L, n_anc, spans in work:
            fake = rng.integers(0, L, size=n_anc)
            for s, e in spans:
                if np.any(np.minimum(np.abs(s - (fake + 1000)), np.abs(fake - e)) <= NEAR):
                    k += 1
        null[i] = k
        if (i + 1) % 1000 == 0:
            print(f"  {i + 1}/{N_PERM} permutations", flush=True)

    p_perm = float((int(np.sum(null >= a_)) + 1) / (N_PERM + 1))
    stored.update(permutation_null_mean=float(null.mean()),
                  permutation_p=p_perm,
                  n_permutations=N_PERM,
                  permutation_null_max=int(null.max()))
    json.dump(report, open(os.path.join(OUT, "precursor_report.json"), "w"), indent=1)
    print(f"observed {a_} of {len(hi)} high scorers near machinery; null mean "
          f"{null.mean():.1f}, max {null.max()}, P = {p_perm:.2g} over {N_PERM} "
          f"permutations", flush=True)


if __name__ == "__main__":
    main()
