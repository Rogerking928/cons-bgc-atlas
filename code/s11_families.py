#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Group the called clusters into families and name them.

Each cluster is sketched as a bottom-k MinHash over the 10-residue words of
its protein products; clusters are joined by single linkage at a containment
of CUT. Families are named from the RefSeq products of their members and from
containment against the curated reference peptides in data/refs.faa.
"""
import collections, csv, json, os
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(HERE, "data", "faa")
OUT = os.path.join(HERE, "out")
K, SKETCH, CUT = 10, 600, 0.50
rng = np.random.default_rng(20260918)
W = rng.integers(1, 2**63 - 1, size=K, dtype=np.uint64)
STOP = {"hypothetical", "protein", "family", "domain-containing", "putative",
        "protein,", "and", "of", "type", "system", "abc"}


def hashes(seq):
    a = np.frombuffer(seq.encode(), dtype=np.uint8).astype(np.uint64)
    if a.size < K:
        return np.zeros(0, dtype=np.uint64)
    from numpy.lib.stride_tricks import sliding_window_view
    h = (sliding_window_view(a, K) * W).sum(axis=1)
    h ^= (h >> np.uint64(29))
    h *= np.uint64(0x94d049bb133111eb)
    return h


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    cds = list(csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"))
    by_cl = collections.defaultdict(list)
    for r in cds:
        by_cl[r["cluster"]].append(r)
    # a cluster is represented by the proteins that carry a biosynthetic
    # domain, not by everything inside its called bounds: the bounds often
    # take in neighbouring genes, and those dilute the sketch enough to split
    # a family (a plasmid-borne cytolysin locus was separated this way).
    MOBILE = {"Transposase", "Integrase", "Recombinase"}
    core_of = collections.defaultdict(list)
    for r in cds:
        ds = set(filter(None, (r["domains"] or "").split("-")))
        if ds and not ds <= MOBILE:
            core_of[r["cluster"]].append(r)
    wanted = collections.defaultdict(set)
    for c, rs in core_of.items():
        for r in rs:
            wanted[r["assembly"]].add(r["protein"])
    seqs = {}
    for a, prots in wanted.items():
        for n, s in read_fasta(os.path.join(FAA, a + ".faa")):
            if n in prots:
                seqs[(a, n)] = s
    print(f"{len(cl)} clusters, {len(seqs)} member proteins", flush=True)

    sk, ids = [], []
    for c in cl:
        hs = [hashes(seqs[(c["assembly"], r["protein"])])
              for r in core_of[c["cluster"]] if (c["assembly"], r["protein"]) in seqs]
        hs = [h for h in hs if h.size]
        keep = np.unique(np.concatenate(hs))[:SKETCH] if hs else np.zeros(0, dtype=np.uint64)
        sk.append(keep); ids.append(c["cluster"])

    n = len(sk)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for i in range(n):
        if sk[i].size == 0:
            continue
        for j in range(i + 1, n):
            if sk[j].size == 0:
                continue
            inter = np.intersect1d(sk[i], sk[j], assume_unique=True).size
            if inter / min(sk[i].size, sk[j].size) >= CUT:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
        if (i + 1) % 200 == 0:
            print(f"  linked {i+1}/{n}", flush=True)

    groups = collections.defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)
    order = sorted(groups.values(), key=len, reverse=True)

    refs = {n_: s for n_, s in read_fasta(os.path.join(HERE, "data", "refs.faa"))}
    ref_h = {n_: set(hashes(s).tolist()) for n_, s in refs.items()}

    fam_of, rows = {}, []
    for k, idx in enumerate(order, 1):
        fid = f"F{k:03d}"
        words = collections.Counter()
        arch = collections.Counter()
        prot_h = set()
        for i in idx:
            c = cl[i]
            fam_of[c["cluster"]] = fid
            for r in by_cl[c["cluster"]]:
                for w in r["product"].lower().replace("/", " ").split():
                    if w not in STOP and len(w) > 3:
                        words[w] += 1
                if r["domains"]:
                    for d in r["domains"].split("-"):
                        arch[d] += 1
            for r in core_of[c["cluster"]]:
                s = seqs.get((c["assembly"], r["protein"]))
                if s is not None:
                    prot_h |= set(hashes(s).tolist())
        hits = []
        for rn, rh in ref_h.items():
            if rh and len(rh & prot_h) / len(rh) >= 0.5:
                hits.append(rn)
        rows.append(dict(
            family=fid, n_clusters=len(idx),
            n_assemblies=len({cl[i]["assembly"] for i in idx}),
            on_plasmid=sum(1 for i in idx if cl[i]["replicon_type"] == "plasmid"),
            kind=collections.Counter(cl[i]["kind"] for i in idx).most_common(1)[0][0],
            median_span=int(np.median([int(cl[i]["span"]) for i in idx])),
            median_cds=int(np.median([int(cl[i]["n_cds"]) for i in idx])),
            # every domain, not the eight commonest: a family with nine or more
            # otherwise loses its rarest, and the rare one is the interesting one
            top_domains=";".join(f"{d}x{c}" for d, c in arch.most_common()),
            ref_hits=",".join(sorted(hits)),
            top_products=";".join(w for w, _ in words.most_common(8))))
    with open(os.path.join(OUT, "families.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    json.dump(fam_of, open(os.path.join(OUT, "cluster_family.json"), "w"))
    print("families:", len(rows), "-> out/families.tsv", flush=True)


if __name__ == "__main__":
    main()
