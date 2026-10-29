#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sequence types and clonal complexes for every genome in the panel.

The seven-locus schemes are pulled from PubMLST and each genome is typed by
exact match of an allele to the assembly, on either strand. Unauthenticated
PubMLST access is capped at records submitted up to 31 December 2024, which is
stated in the methods; a genome whose alleles are all newer than that comes
back untyped rather than mistyped.
"""
import collections, csv, glob, gzip, io, json, os, sys, time, urllib.request

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FNA = os.path.join(HERE, "data", "fna")
MLST = os.path.join(HERE, "data", "mlst")
# 五個有七基因座方案的物種。lugdunensis 的方案不在 PubMLST，在 Institut Pasteur
# 的 BIGSdb（同一套 REST API，不同主機），所以每個物種連 base URL 一起記。
PUB = "https://rest.pubmlst.org/db"
PAST = "https://bigsdb.pasteur.fr/api/db"
DB = {"epidermidis": (PUB, "pubmlst_sepidermidis_seqdef"),
      "capitis": (PUB, "pubmlst_scapitis_seqdef"),
      "haemolyticus": (PUB, "pubmlst_shaemolyticus_seqdef"),
      "hominis": (PUB, "pubmlst_shominis_seqdef"),
      "lugdunensis": (PAST, "pubmlst_staphlugdunensis_seqdef")}
K = 31
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def get(url, binary=False):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as fh:
                return fh.read() if binary else fh.read().decode()
        except Exception as e:
            if attempt == 4:
                raise
            time.sleep(4 * (attempt + 1))


def fetch_scheme(sp):
    d = os.path.join(MLST, sp)
    os.makedirs(d, exist_ok=True)
    host, dbname = DB[sp]
    base = f"{host}/{dbname}"
    meta = json.loads(get(f"{base}/schemes/1"))
    loci = [u.rsplit("/", 1)[1] for u in meta["loci"]]
    prof = os.path.join(d, "profiles.tsv")
    if not os.path.exists(prof):
        open(prof, "w").write(get(f"{base}/schemes/1/profiles_csv"))
    for loc in loci:
        p = os.path.join(d, loc + ".fas")
        if not os.path.exists(p):
            try:
                open(p, "w").write(get(f"{base}/loci/{loc}/alleles_fasta"))
            except Exception as e:
                # 未認證存取只到 2024-12-31；某些方案（S. capitis）所有等位基因都比
                # 那個日期新，整個 locus 回 404。該物種就不分型，不要硬猜。
                os.path.exists(p) and os.remove(p)
                print(f"  {sp}/{loc}: no public alleles ({e}) -> 略過此物種",
                      flush=True)
                return None
            time.sleep(0.4)
    print(f"{sp}: {len(loci)} loci, profiles {os.path.getsize(prof)} bytes", flush=True)
    return loci


def read_fasta_text(txt):
    n, b = None, []
    for line in txt.split("\n"):
        if line.startswith(">"):
            if n:
                yield n, "".join(b)
            n, b = line[1:].strip(), []
        else:
            b.append(line.strip().upper())
    if n:
        yield n, "".join(b)


def hashes(seq):
    a = np.frombuffer(seq.encode(), dtype=np.uint8).astype(np.uint64)
    if a.size < K:
        return np.zeros(0, dtype=np.uint64)
    from numpy.lib.stride_tricks import sliding_window_view
    w = np.arange(K, dtype=np.uint64) * np.uint64(7)
    h = (sliding_window_view(a, K) << (w % np.uint64(13))).sum(axis=1)
    h ^= h >> np.uint64(29)
    h *= np.uint64(0x94d049bb133111eb)
    return h


def seed_hash(s):
    return hashes(s[:K])[:1]


def main():
    alleles = {}
    usable = []
    for sp in list(DB):
        loci = fetch_scheme(sp)
        if not loci:
            continue
        usable.append(sp)
        for loc in loci:
            txt = open(os.path.join(MLST, sp, loc + ".fas")).read()
            for name, seq in read_fasta_text(txt):
                num = name.split("_")[-1]
                alleles.setdefault((sp, loc), []).append((num, seq))
        print(f"  {sp}: " + ", ".join(f"{loc}={len(alleles[(sp, loc)])}" for loc in loci),
              flush=True)

    profiles = {}
    for sp in usable:
        rows = list(csv.reader(open(os.path.join(MLST, sp, "profiles.tsv")), delimiter="\t"))
        head = rows[0]
        li = [i for i, h in enumerate(head) if h not in ("ST", "clonal_complex", "species")]
        cc_i = head.index("clonal_complex") if "clonal_complex" in head else None
        for r in rows[1:]:
            key = tuple(r[i] for i in li)
            profiles[(sp, key)] = (r[0], r[cc_i] if cc_i is not None else "")
        profiles[(sp, "loci")] = [head[i] for i in li]
        print(f"  {sp}: {sum(1 for k in profiles if k[0] == sp) - 1} profiles", flush=True)

    mat = {r["assembly"]: r["species"] for r in
           csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    rows = []
    for n, (acc, sp) in enumerate(sorted(mat.items()), 1):
        if sp not in usable:      # 沒有（可公開取用的）方案就不分型
            continue
        seqs, name, buf = {}, None, []
        with gzip.open(os.path.join(FNA, acc + ".fna.gz"), "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(buf)
                    name, buf = line[1:].split()[0], []
                else:
                    buf.append(line.strip().upper())
        if name:
            seqs[name] = "".join(buf)
        genome = "".join(seqs.values())
        hs = np.unique(np.concatenate([hashes(genome), hashes(rc(genome))]))
        call = {}
        for loc in profiles[(sp, "loci")]:
            got = ""
            for num, seq in alleles.get((sp, loc), []):
                sh = seed_hash(seq)
                if sh.size and np.searchsorted(hs, sh[0]) < hs.size and \
                        hs[min(np.searchsorted(hs, sh[0]), hs.size - 1)] == sh[0]:
                    if seq in genome or rc(seq) in genome:
                        got = num
                        break
            call[loc] = got
        key = tuple(call[loc] for loc in profiles[(sp, "loci")])
        st, cc = profiles.get((sp, key), ("", ""))
        rows.append(dict(assembly=acc, species=sp,
                         **{loc: call[loc] for loc in profiles[(sp, "loci")]},
                         ST=st, clonal_complex=cc,
                         complete=int(all(call.values()))))
        if n % 50 == 0:
            typed = sum(1 for r in rows if r["ST"])
            print(f"  {n}/{len(mat)} genomes, {typed} with an ST", flush=True)

    with open(os.path.join(OUT, "mlst.csv"), "w", newline="") as fh:
        keys = ["assembly", "species", "ST", "clonal_complex", "complete"]
        loci_all = sorted({k for r in rows for k in r if k not in keys})
        w = csv.DictWriter(fh, fieldnames=keys + loci_all, restval="")
        w.writeheader(); w.writerows(rows)
    typed = sum(1 for r in rows if r["ST"])
    print(f"typed {typed}/{len(rows)}; top STs:", flush=True)
    for sp in usable:
        c = collections.Counter(r["ST"] for r in rows if r["species"] == sp and r["ST"])
        print(f"  S. {sp}: " + ", ".join(f"ST{s} x{n}" for s, n in c.most_common(6)),
              flush=True)


if __name__ == "__main__":
    main()
