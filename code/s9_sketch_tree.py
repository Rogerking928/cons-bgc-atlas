#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Whole-proteome k-mer sketches -> Mash-style distance matrix for the genus.

Same bottom-k MinHash over 9-residue protein k-mers as the 443-genome version,
so the two are directly comparable; the pairwise loop is parallelised because
2,007 genomes is ~2.0M pairs rather than ~98k.
It is a similarity dendrogram, not a phylogeny, and the caption must say so.
"""
import json, os, sys, glob
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from multiprocessing import Pool

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(HERE, "data", "faa")
OUT = os.path.join(HERE, "out")
K, SKETCH = 9, 2000
rng = np.random.default_rng(20260918)
W = rng.integers(1, 2**63 - 1, size=K, dtype=np.uint64)


def sketch(path):
    hs, seqs, cur = [], [], []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if cur:
                    seqs.append(b"".join(cur)); cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append(b"".join(cur))
    # Reduce to the bottom-k every 250 proteins instead of concatenating the whole
    # proteome's k-mer hashes: the full array is ~10^6 uint64 per genome and this box
    # has 11 GB shared with everything else.
    keep = np.zeros(0, dtype=np.uint64)
    for s in seqs:
        if len(s) < K:
            continue
        a = np.frombuffer(s, dtype=np.uint8).astype(np.uint64)
        h = (sliding_window_view(a, K) * W).sum(axis=1)
        h ^= (h >> np.uint64(29))
        h *= np.uint64(0x94d049bb133111eb)
        hs.append(h)
        if len(hs) >= 250:
            keep = np.unique(np.concatenate(hs + [keep]))[:SKETCH]
            hs = []
    if hs:
        keep = np.unique(np.concatenate(hs + [keep]))[:SKETCH]
    return keep


SK = None


def _init(sk):
    global SK
    SK = sk


def _row(i):
    si = SK[i]
    n = len(SK)
    out = np.zeros(n)
    for j in range(i + 1, n):
        sj = SK[j]
        merged = np.union1d(si, sj)[:SKETCH]
        inter = np.intersect1d(merged, si, assume_unique=True)
        inter = np.intersect1d(inter, sj, assume_unique=True).size
        jac = inter / merged.size if merged.size else 0.0
        out[j] = 1.0 if jac <= 0 else min(-np.log(2 * jac / (1 + jac)) / K, 1.0)
    return i, out


def main():
    accs = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(FAA, "*.faa")))
    print(f"proteomes: {len(accs)}", flush=True)
    sk = []
    for i, a in enumerate(accs, 1):
        sk.append(sketch(os.path.join(FAA, a + ".faa")))
        if i % 200 == 0:
            print(f"  sketched {i}/{len(accs)}", flush=True)
    n = len(accs)
    D = np.zeros((n, n))
    nproc = max(1, min(os.cpu_count() or 2, 6))
    print(f"pairwise on {nproc} procs ({n*(n-1)//2} pairs)", flush=True)
    with Pool(nproc, initializer=_init, initargs=(sk,)) as pool:
        for done, (i, row) in enumerate(pool.imap_unordered(_row, range(n), chunksize=8), 1):
            D[i] = np.maximum(D[i], row)
            if done % 200 == 0:
                print(f"  rows {done}/{n}", flush=True)
    D = np.maximum(D, D.T)
    np.save(os.path.join(OUT, "tree_D.npy"), D)
    json.dump({"accs": accs}, open(os.path.join(OUT, "tree_meta.json"), "w"))
    iu = np.triu_indices(n, 1)
    print(f"done: n={n}  median distance {float(np.median(D[iu])):.4f}", flush=True)


if __name__ == "__main__":
    main()
