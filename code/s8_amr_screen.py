#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-genome acquired-resistance content, by k-mer containment of the
AMRFinderPlus reference proteins.

Containment = the fraction of the reference protein's 10-residue words that
occur in the genome's proteome. A true carrier (usually >95% identity) lands
near 1.0; a homologue that is not the determinant falls close to 0. The call
threshold is 0.5, with 0.3 and 0.7 kept for a sensitivity analysis.
"""
import csv, glob, os, re, sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(HERE, "data", "faa")
OUT = os.path.join(HERE, "out")
AMR = os.path.join(HERE, "data", "amr", "AMRProt.fa")
K = 10
# 葡萄球菌版：加上 β-內醯胺（mecA/mecC/blaZ）、mupirocin（mupA=ileS2）、
# fusidic acid（fusB/fusC）、trimethoprim（dfrG/dfrK）與四級銨消毒劑（qacA/qacC）
CLASSES = {"GLYCOPEPTIDE", "MACROLIDE", "LINCOSAMIDE", "STREPTOGRAMIN",
           "AMINOGLYCOSIDE", "TETRACYCLINE", "PHENICOL", "OXAZOLIDINONE",
           "BETA-LACTAM", "FUSIDIC_ACID", "TRIMETHOPRIM", "MUPIROCIN",
           "QUATERNARY_AMMONIUM", "PLEUROMUTILIN",
           "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN", "LINCOSAMIDE/MACROLIDE",
           "MACROLIDE/STREPTOGRAMIN", "LINCOSAMIDE/STREPTOGRAMIN",
           "LINCOSAMIDE/PLEUROMUTILIN/STREPTOGRAMIN",
           "LINCOSAMIDE/OXAZOLIDINONE/PHENICOL/PLEUROMUTILIN/STREPTOGRAMIN"}
FAM_RE = re.compile(r"^(mecA|mecB|mecC|mecD|blaZ|mupA|mupB|ileS2|fusB|fusC|"
                    r"fusD|dfrG|dfrK|dfrA|qacA|qacB|qacC|qacG|qacJ|smr|"
                    r"van[A-Z]|erm|msr|mph|mef|lsa|lnu|vga|vat|sal|aac|aph|ant|"
                    r"aad|str|tet|optrA|poxtA|cfr|fex|cat|sat|spc)", re.I)


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def load_refs():
    refs = {}
    for h, s in read_fasta(AMR):
        f = h.split("|")
        if len(f) < 10:
            continue
        family, cls = f[4], f[8].upper()
        if not FAM_RE.match(family):
            continue
        if not (set(cls.split("/")) & CLASSES):
            continue
        refs.setdefault(family, []).append(s)
    # cap the number of alleles per family so the scan stays linear
    return {k: v[:40] for k, v in refs.items()}


REFS = None


def _init(r):
    global REFS
    REFS = r


def scan(acc):
    kmers = set()
    for _, s in read_fasta(os.path.join(FAA, acc + ".faa")):
        for j in range(len(s) - K + 1):
            kmers.add(s[j:j + K])
    row = {"assembly": acc}
    for fam, seqs in REFS.items():
        best = 0.0
        for s in seqs:
            tot = len(s) - K + 1
            if tot <= 0:
                continue
            hit = sum(1 for j in range(tot) if s[j:j + K] in kmers)
            v = hit / tot
            if v > best:
                best = v
        row[fam] = round(best, 3)
    return row


def main():
    refs = load_refs()
    print(f"{len(refs)} families, {sum(len(v) for v in refs.values())} alleles",
          file=sys.stderr, flush=True)
    accs = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(FAA, "*.faa")))
    print(f"{len(accs)} proteomes", file=sys.stderr, flush=True)
    rows = []
    nproc = max(1, min(os.cpu_count() or 2, 8))
    with Pool(nproc, initializer=_init, initargs=(refs,)) as pool:
        for i, row in enumerate(pool.imap_unordered(scan, accs, chunksize=4), 1):
            rows.append(row)
            if i % 100 == 0:
                print(f"  {i}/{len(accs)}", file=sys.stderr, flush=True)
    rows.sort(key=lambda r: r["assembly"])
    fams = sorted(refs)
    with open(os.path.join(OUT, "amr_per_genome.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["assembly"] + fams)
        w.writeheader(); w.writerows(rows)
    print("wrote out/amr_per_genome.csv", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
