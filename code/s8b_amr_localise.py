#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Put every called resistance determinant on its replicon.

For each assembly, the determinants called by s8 are matched back to the single
protein that carries them (the protein holding the most 10-residue words of the
best-matching reference allele), and that protein is looked up in the CDS table
to give its replicon. Output: out/amr_located.csv
"""
import collections, csv, gzip, os, sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(HERE, "data", "faa")
OUT = os.path.join(HERE, "out")
AMR = os.path.join(HERE, "data", "amr", "AMRProt.fa")
K, CALL = 10, 0.5


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


REFS = POS = None


def _init(refs, pos):
    global REFS, POS
    REFS, POS = refs, pos


def locate(acc):
    fams = POS.get(acc, [])
    if not fams:
        return acc, []
    where = collections.defaultdict(collections.Counter)
    prot_kmers = {}
    for h, s in read_fasta(os.path.join(FAA, acc + ".faa")):
        pid = h.split()[0]
        prot_kmers[pid] = s
    index = collections.defaultdict(list)
    for pid, s in prot_kmers.items():
        for j in range(len(s) - K + 1):
            index[s[j:j + K]].append(pid)
    out = []
    for fam in fams:
        best, best_prot, best_cov = None, None, 0.0
        for allele in REFS[fam]:
            tot = len(allele) - K + 1
            if tot <= 0:
                continue
            c = collections.Counter()
            for j in range(tot):
                for pid in index.get(allele[j:j + K], ()):
                    c[pid] += 1
            if not c:
                continue
            pid, hit = c.most_common(1)[0]
            if hit / tot > best_cov:
                best_cov, best_prot, best = hit / tot, pid, allele
        if best_prot and best_cov >= CALL:
            out.append((fam, best_prot, round(best_cov, 3)))
    return acc, out


def main():
    amr_rows = list(csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv"))))
    fams_all = [k for k in amr_rows[0] if k != "assembly"]
    pos = {r["assembly"]: [f for f in fams_all if float(r[f] or 0) >= CALL]
           for r in amr_rows}
    want = {f for v in pos.values() for f in v}
    refs = collections.defaultdict(list)
    for h, s in read_fasta(AMR):
        f = h.split("|")
        if len(f) >= 10 and f[4] in want:
            refs[f[4]].append(s)
    refs = {k: v[:40] for k, v in refs.items()}
    print(f"{len(want)} determinants to place in {sum(1 for v in pos.values() if v)} genomes",
          file=sys.stderr, flush=True)

    loc = {}
    with Pool(min(8, os.cpu_count() or 2), initializer=_init, initargs=(refs, pos)) as pool:
        for i, (acc, hits) in enumerate(pool.imap_unordered(locate, sorted(pos), chunksize=4), 1):
            loc[acc] = hits
            if i % 100 == 0:
                print(f"  {i}/{len(pos)}", file=sys.stderr, flush=True)

    rep = {}
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rep[(f[0], f[5])] = (f[1], f[2], f[4], f[6], f[7])
    with open(os.path.join(OUT, "amr_located.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["assembly", "determinant", "protein", "containment",
                    "replicon", "replicon_type", "plasmid", "start", "end"])
        n = 0
        for acc in sorted(loc):
            for fam, pid, cov in loc[acc]:
                r = rep.get((acc, pid), ("", "", "", "", ""))
                w.writerow([acc, fam, pid, cov, *r]); n += 1
    print("wrote out/amr_located.csv", n, "rows", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
