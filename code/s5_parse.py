#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build two tables from the GFFs and the domtblout files.

out/cds.tsv.gz      every CDS: assembly, replicon, replicon type/length, protein, coords, product
out/domains.tsv.gz  every Pfam domain hit at the gathering threshold
"""
import glob, gzip, os, re, sys, urllib.parse

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GFF = os.path.join(HERE, "data", "gff")
DOM = os.path.join(HERE, "out", "dom")
OUT = os.path.join(HERE, "out")

# Pfam accession -> short domain label used throughout the paper
LABEL = {
    "PF00109": "KS", "PF02801": "KS_C", "PF16197": "KS_assoc",
    "PF00501": "A", "PF13193": "A_C", "PF00550": "PCP/ACP", "PF00668": "C",
    "PF00698": "AT", "PF08659": "KR", "PF00975": "TE",
    "PF07993": "ER/NAD", "PF08242": "MT", "PF00202": "Aminotran_3",
    "PF05147": "LanC", "PF04738": "LanB_N", "PF14028": "LanB_C",
    "PF13575": "LanM_cyc", "PF16934": "Mersacidin", "PF04604": "L_biotic_typeA",
    "PF03412": "C39", "PF02624": "YcaO", "PF05402": "PqqD",
    "PF01721": "Bacteriocin_II", "PF10439": "Bacteriocin_IIc",
    "PF09221": "Bacteriocin_IId", "PF12173": "BacteriocIIc_cy",
    "PF11758": "Bacteriocin_IIi", "PF08129": "Antimicrobial17",
    "PF17508": "MccV", "PF08951": "EntA_Immun", "PF00082": "Peptidase_S8",
    "PF04055": "Radical_SAM",
    "PF02052": "Gallidermin", "PF08130": "Antimicrobial18",
    "PF04369": "Lactococcin_972", "PF19402": "Lanthi_prec",
    "PF00589": "Integrase", "PF13495": "Integrase", "PF01609": "Transposase",
    "PF00665": "Transposase", "PF13683": "Transposase",
    "PF02914": "Recombinase", "PF07508": "Recombinase",
    "PF03050": "Transposase", "PF01548": "Transposase",
}


def attrs(s):
    d = {}
    for kv in s.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = urllib.parse.unquote(v)
    return d


def parse_gff(path):
    rep, rows = {}, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "region":
                a = attrs(f[8])
                rep[f[0]] = (a.get("genome", "unknown"), int(f[4]),
                             a.get("plasmid-name", ""))
            elif f[2] == "CDS":
                a = attrs(f[8])
                pid = a.get("protein_id")
                if not pid:
                    continue
                rows.append((f[0], pid, int(f[3]), int(f[4]), f[6],
                             a.get("product", ""), a.get("locus_tag", "")))
    return rep, rows


def parse_dom(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            pf = f[4].split(".")[0]
            out.append((f[0], pf, LABEL.get(pf, pf), int(f[17]), int(f[18]),
                        f[6], f[7], int(f[2])))
    return out


def main():
    accs = sorted(os.path.basename(p)[:-7] for p in glob.glob(os.path.join(DOM, "*.domtbl")))
    print(len(accs), "assemblies", flush=True)
    fc = gzip.open(os.path.join(OUT, "cds.tsv.gz"), "wt")
    fd = gzip.open(os.path.join(OUT, "domains.tsv.gz"), "wt")
    fc.write("assembly\treplicon\treplicon_type\treplicon_len\tplasmid\tprotein\tstart\tend\tstrand\tlocus_tag\tproduct\n")
    fd.write("assembly\tprotein\tpfam\tdomain\tali_from\tali_to\tievalue\tscore\tprot_len\n")
    for i, a in enumerate(accs):
        g = os.path.join(GFF, a + ".gff.gz")
        if not os.path.exists(g):
            print("no gff", a, flush=True); continue
        rep, rows = parse_gff(g)
        for r, pid, s, e, st, prod, lt in rows:
            ty, ln, pl = rep.get(r, ("unknown", 0, ""))
            fc.write(f"{a}\t{r}\t{ty}\t{ln}\t{pl}\t{pid}\t{s}\t{e}\t{st}\t{lt}\t{prod}\n")
        for pid, pf, lab, af, at, ev, sc, plen in parse_dom(os.path.join(DOM, a + ".domtbl")):
            fd.write(f"{a}\t{pid}\t{pf}\t{lab}\t{af}\t{at}\t{ev}\t{sc}\t{plen}\n")
        if (i + 1) % 100 == 0:
            print(i + 1, flush=True)
    fc.close(); fd.close()
    print("wrote cds.tsv.gz and domains.tsv.gz", flush=True)


if __name__ == "__main__":
    main()
