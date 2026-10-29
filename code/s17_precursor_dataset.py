#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Training set for the precursor detector.

Positives
  a. UniProt reviewed bacteriocins, antimicrobial peptides and lantibiotics of
     the Bacillota, 10-120 residues: curated, and mostly from outside this panel
  b. every short protein in the panel that a RiPP structural Pfam profile
     already recognises

Negatives
  short proteins of the panel with no detected domain that sit more than 20 kb
  from any maturation or export anchor, share less than SIM of their 10-residue
  words with any positive, and are not named as a bacteriocin by RefSeq. The
  distance rule alone was not enough: the first version of this set swallowed 17
  distinct Blp-family bacteriocins, which the detector then re-found and which
  were reported as discoveries.

Held-out annotation set
  short proteins that RefSeq names a bacteriocin, lantibiotic, lanthipeptide,
  microcin or similar but that no Pfam profile in the library recognises. These
  enter neither side of the training set. They are exactly the peptides the
  profile library misses, with a label that owes nothing to the model, so recall
  on them is an honest test of what the detector adds.

Sequences are deduplicated and then greedily clustered on 5-mer Jaccard at 0.5,
and the cluster is the unit of cross-validation: a near-duplicate of a training
peptide must not be scored as a held-out success.
"""
import collections, csv, gzip, json, os, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FAA = os.path.join(HERE, "data", "faa")
DATA = os.path.join(HERE, "data")
MAXLEN = 120
NEG_MAXLEN = 105
ANCHOR_DOMS = {"C39", "LanM_cyc", "LanC", "YcaO", "LanB_N", "LanB_C",
               "Peptidase_S8", "Radical_SAM"}
PRECURSOR_DOMS = {"Mersacidin", "Bacteriocin_II", "Bacteriocin_IIc", "Bacteriocin_IId",
                  "Bacteriocin_IIi", "BacteriocIIc_cy", "Antimicrobial17",
                  "L_biotic_typeA", "MccV",
                  # 葡萄球菌版新增
                  "Gallidermin", "Antimicrobial18", "Lanthi_prec",
                  "Lactococcin_972"}
FAR_BP = 20000
SIM = 0.30                 # 10-mer containment against any positive
NAME_RE = __import__("re").compile(
    r"bacteriocin|lantibiotic|lanthipeptide|microcin|enterocin|plantaricin|"
    r"lactococcin|sactipeptide|thiopeptide|linear azol|autoinducer|pheromone",
    __import__("re").I)
KW = 10
QUERIES = [
    "(taxonomy_id:1239) AND (reviewed:true) AND (keyword:KW-0078) AND (length:[10 TO 120])",
    "(taxonomy_id:1239) AND (reviewed:true) AND (keyword:KW-0929) AND (length:[10 TO 120])",
    "(reviewed:true) AND (keyword:KW-0425) AND (length:[10 TO 120])",
]


def uniprot(query):
    url = ("https://rest.uniprot.org/uniprotkb/search?format=tsv&size=500"
           "&fields=accession,protein_name,organism_name,length,sequence&query="
           + urllib.parse.quote(query))
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as fh:
                rows = fh.read().decode().strip().split("\n")
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(5)
    out = []
    for line in rows[1:]:
        f = line.split("\t")
        if len(f) >= 5 and f[4].isalpha():
            out.append((f[0], f[1], f[2], f[4]))
    return out


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


def cluster(seqs, k=5, cut=0.5):
    """Greedy 5-mer Jaccard clustering; centroid order is longest first."""
    def kmers(s):
        return {s[i:i + k] for i in range(len(s) - k + 1)} or {s}
    order = sorted(seqs, key=len, reverse=True)
    cents, out = [], {}
    for s in order:
        ks = kmers(s)
        hit = None
        for ci, (cs, cks) in enumerate(cents):
            j = len(ks & cks) / len(ks | cks)
            if j >= cut:
                hit = ci
                break
        if hit is None:
            cents.append((s, ks))
            hit = len(cents) - 1
        out[s] = hit
    return out


def main():
    # --- panel domains and coordinates
    dom = collections.defaultdict(dict)
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            dom[f[0]].setdefault(f[1], set()).add(f[3])
    pos_of, anchors, prod_of = {}, collections.defaultdict(list), {}
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, rep, prot, s, e = f[0], f[1], f[5], int(f[6]), int(f[7])
            pos_of[(a, prot)] = (rep, s, e)
            prod_of[(a, prot)] = f[10]
            if dom.get(a, {}).get(prot, set()) & ANCHOR_DOMS:
                anchors[(a, rep)].append((s, e))

    def far_from_anchor(a, prot):
        rep, s, e = pos_of[(a, prot)]
        for (as_, ae) in anchors.get((a, rep), ()):
            if s <= ae and e >= as_:
                return False
            if min(abs(s - ae), abs(as_ - e)) <= FAR_BP:
                return False
        return True

    rows, holdout = [], []
    seen = {}
    pos_words = set()
    # --- positives from UniProt
    up = {}
    for q in QUERIES:
        for acc, name, org, seq in uniprot(q):
            up[acc] = (name, org, seq)
        time.sleep(1)
    for acc, (name, org, seq) in up.items():
        if seq in seen:
            continue
        seen[seq] = 1
        rows.append(dict(seq=seq, label=1, source="uniprot", id=acc,
                         note=f"{name} | {org}"))
        pos_words |= {seq[i:i + KW] for i in range(len(seq) - KW + 1)}
    n_up = len(rows)

    # --- positives and negatives from the panel
    n_panel_pos = n_neg = n_drop = 0
    for f in sorted(os.listdir(FAA)):
        a = f[:-4]
        for prot, seq in read_fasta(os.path.join(FAA, f)):
            if not (10 <= len(seq) <= MAXLEN) or (a, prot) not in pos_of:
                continue
            ds = dom.get(a, {}).get(prot, set())
            if ds & PRECURSOR_DOMS:
                if seq not in seen:
                    seen[seq] = 1
                    rows.append(dict(seq=seq, label=1, source="panel", id=f"{a}:{prot}",
                                     note=",".join(sorted(ds))))
                    pos_words |= {seq[i:i + KW] for i in range(len(seq) - KW + 1)}
                    n_panel_pos += 1
            elif not ds and len(seq) <= NEG_MAXLEN and far_from_anchor(a, prot):
                if seq in seen:
                    continue
                w = [seq[i:i + KW] for i in range(len(seq) - KW + 1)]
                if w and sum(1 for x in w if x in pos_words) / len(w) >= SIM:
                    n_drop += 1
                    continue
                if NAME_RE.search(prod_of.get((a, prot), "")):
                    seen[seq] = 1
                    holdout.append(dict(seq=seq, source="panel", id=f"{a}:{prot}",
                                        product=prod_of[(a, prot)]))
                    continue
                seen[seq] = 1
                rows.append(dict(seq=seq, label=0, source="panel", id=f"{a}:{prot}",
                                 note=""))
                n_neg += 1

    # only the positives need clustering: they are the side where a near
    # duplicate would leak between folds, and clustering 50k negatives greedily
    # would take longer than the rest of the pipeline put together
    cl = cluster([r["seq"] for r in rows if r["label"] == 1])
    nid = 10 ** 6
    for r in rows:
        if r["label"] == 1:
            r["cluster"] = cl[r["seq"]]
        else:
            r["cluster"] = nid
            nid += 1
    p = os.path.join(OUT, "precursor_trainset.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["seq", "label", "source", "id", "note", "cluster"])
        w.writeheader(); w.writerows(rows)
    ph = os.path.join(OUT, "precursor_holdout.csv")
    with open(ph, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["seq", "source", "id", "product"])
        w.writeheader(); w.writerows(holdout)
    print(f"positives: {n_up} from UniProt, {n_panel_pos} from the panel; "
          f"negatives: {n_neg} ({n_drop} dropped for resembling a positive, "
          f"{len(holdout)} held out because RefSeq names them a bacteriocin); "
          f"positive clusters: {len(set(cl.values()))}")
    print("->", ph)
    print("->", p)


if __name__ == "__main__":
    main()
