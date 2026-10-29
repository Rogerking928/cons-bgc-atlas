#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Specificity-conferring codes of the adenylation domains in the panel.

Every A domain found in an NRPS or hybrid cluster is aligned globally to the
phenylalanine-activating domain of gramicidin S synthase 1 (GrsA, UniProt
P0C061) with BLOSUM62 and affine gaps, and the ten residues that line the
substrate pocket are read off at the GrsA positions 235, 236, 239, 278, 299,
301, 322, 330, 331 and 517. The signature is then compared with the published
consensus codes (Stachelhaus, Mootz and Marahiel, Chemistry & Biology
1999;6:493-505, Table 1).

The control is GrsA itself: aligned to itself the pipeline must return
DAWTIAAICK, the published signature of PheA.
"""
import collections, csv, gzip, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FAA = os.path.join(HERE, "data", "faa")
REF = os.path.join(HERE, "data", "ref")
POSITIONS = [235, 236, 239, 278, 299, 301, 322, 330, 331, 517]
WOBBLE = {278, 299}
GAP_OPEN, GAP_EXTEND = -11, -1
PAD = 60


def read_fasta(path):
    n, b = None, []
    for line in open(path):
        if line.startswith(">"):
            if n:
                yield n, "".join(b)
            n, b = line[1:].split()[0], []
        else:
            b.append(line.strip())
    if n:
        yield n, "".join(b)


def blosum62():
    rows = [l.split() for l in open(os.path.join(REF, "BLOSUM62"))
            if l.strip() and not l.startswith("#")]
    head = rows[0]
    m = {}
    for r in rows[1:]:
        for c, v in zip(head, r[1:]):
            m[(r[0], c)] = int(v)
    return m


BL = blosum62()


def sc(a, b):
    return BL.get((a, b), BL.get(("*", "*"), -4))


def nw(a, b):
    """Needleman-Wunsch with affine gaps; returns the aligned pair of strings."""
    n, m = len(a), len(b)
    NEG = -10 ** 7
    M = np.full((n + 1, m + 1), NEG, dtype=np.int32)
    Ix = np.full((n + 1, m + 1), NEG, dtype=np.int32)   # gap in b
    Iy = np.full((n + 1, m + 1), NEG, dtype=np.int32)   # gap in a
    P = np.zeros((n + 1, m + 1, 3), dtype=np.int8)
    M[0, 0] = 0
    for i in range(1, n + 1):
        Ix[i, 0] = GAP_OPEN + (i - 1) * GAP_EXTEND
    for j in range(1, m + 1):
        Iy[0, j] = GAP_OPEN + (j - 1) * GAP_EXTEND
    for i in range(1, n + 1):
        ai = a[i - 1]
        row = np.array([sc(ai, bj) for bj in b], dtype=np.int32)
        for j in range(1, m + 1):
            best_prev = max(M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1])
            src = int(np.argmax([M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1]]))
            M[i, j] = best_prev + row[j - 1]
            P[i, j, 0] = src
            o, e = M[i - 1, j] + GAP_OPEN, Ix[i - 1, j] + GAP_EXTEND
            Ix[i, j], P[i, j, 1] = (o, 0) if o >= e else (e, 1)
            o, e = M[i, j - 1] + GAP_OPEN, Iy[i, j - 1] + GAP_EXTEND
            Iy[i, j], P[i, j, 2] = (o, 0) if o >= e else (2, 2) if False else (e, 2)
    i, j = n, m
    state = int(np.argmax([M[n, m], Ix[n, m], Iy[n, m]]))
    A, B = [], []
    while i > 0 or j > 0:
        if state == 0:
            A.append(a[i - 1]); B.append(b[j - 1])
            state = int(P[i, j, 0]); i -= 1; j -= 1
        elif state == 1:
            A.append(a[i - 1]); B.append("-")
            state = 0 if P[i, j, 1] == 0 else 1
            i -= 1
        else:
            A.append("-"); B.append(b[j - 1])
            state = 0 if P[i, j, 2] == 0 else 2
            j -= 1
        if i == 0 and j > 0:
            A.append("-" * j); B.append(b[:j][::-1]); j = 0
        elif j == 0 and i > 0:
            A.append(a[:i][::-1]); B.append("-" * i); i = 0
    return "".join(A)[::-1], "".join(B)[::-1]


def signature(query, ref):
    """Read the ten GrsA pocket positions off a global alignment."""
    qa, ra = nw(query, ref)
    qi = ri = 0
    at = {}
    for q, r in zip(qa, ra):
        if r != "-":
            ri += 1
        if q != "-":
            qi += 1
        if r != "-" and ri in POSITIONS:
            at[ri] = q if q != "-" else "."
    return "".join(at.get(p, ".") for p in POSITIONS)


def main():
    ref = dict(read_fasta(os.path.join(HERE, "data", "refs.faa")))
    grsa = ref["P0C061"]
    consensus = []
    for line in open(os.path.join(REF, "stachelhaus_table1.tsv")):
        if line.startswith("#") or line.startswith("substrate"):
            continue
        f = line.rstrip("\n").split("\t")
        consensus.append((f[0], f[1]))

    print("control: GrsA against itself ->", signature(grsa[:600], grsa[:600]), flush=True)

    fam = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = {r["cluster"]: r for r in csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t")}
    want = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"):
        if cl[r["cluster"]]["kind"] != "RiPP" and "A" in (r["domains"] or "").split("-"):
            want[r["assembly"]].append((r["cluster"], r["protein"], r["product"], r["locus_tag"]))
    dom = collections.defaultdict(list)
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3] == "A":
                dom[(f[0], f[1])].append((int(f[4]), int(f[5])))

    rows = []
    for a in sorted(want):
        seqs = dict(read_fasta(os.path.join(FAA, a + ".faa")))
        for cluster, prot, product, tag in want[a]:
            s = seqs.get(prot)
            if not s:
                continue
            for k, (lo, hi) in enumerate(sorted(dom[(a, prot)]), 1):
                sub = s[max(0, lo - PAD):min(len(s), hi + 400)]
                sig = signature(sub, grsa[:600])
                best, bid = "", -1
                for name, code in consensus:
                    idn = sum(1 for x, y in zip(sig, code) if x == y and x != ".")
                    if idn > bid:
                        best, bid = name, idn
                rows.append(dict(assembly=a, cluster=cluster, family=fam[cluster],
                                 protein=prot, locus_tag=tag, product=product,
                                 domain=f"A{k}", start=lo, end=hi,
                                 signature=sig, aligned=sum(1 for c in sig if c != "."),
                                 closest=best, identity=bid))
                print(f"  {a} {fam[cluster]} {prot} A{k}: {sig}  closest {best} {bid}/10",
                      flush=True)
    with open(os.path.join(OUT, "adomain_codes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} A domains -> out/adomain_codes.csv", flush=True)


if __name__ == "__main__":
    main()
