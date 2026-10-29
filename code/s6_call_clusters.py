#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Call biosynthetic gene clusters from the domain table and the CDS table.

A cluster is a run of domain-carrying CDS on one replicon, chained while the
gap between neighbours stays below MAXGAP, that contains at least one anchor
(a maturation enzyme) and at least two domain-carrying proteins in total.

Writes out/clusters.tsv (one row per cluster) and out/cluster_cds.tsv
(every CDS inside a cluster, with its domain architecture).
"""
import collections, csv, gzip, os, re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
MAXGAP = 10000          # nt between neighbouring domain-carrying CDS
FLANK = 5000            # 核心（anchor）基因往外各留 5 kb 收前驅肽／輸出／免疫基因

# 初級代謝的「假 anchor」：DltA 帶 AMP-binding、BioB 帶 Radical_SAM，
# 串連時會把 dlt 操縱子或 bio 操縱子一起吞進來，界定基因簇時必須排除。
# （2026-09-20 踩過：S. lugdunensis 的 F003 因此長到 26.7 kb、含整個 dlt 操縱子。）
PRIMARY_RE = re.compile(
    r"D-alanyl|poly\(phosphoribitol\)|lipoyl synthase|biotin synthase|"
    r"dihydroxyacetone|acyl-CoA|fatty acid|acetyl-CoA carboxylase|"
    r"3-oxoacyl-\[acyl-carrier-protein\]|beta-ketoacyl-ACP", re.I)

# maturation enzymes: their presence defines a cluster
ANCHOR = {"LanC", "LanM_cyc", "LanB_N", "LanB_C", "YcaO", "C39",
          "A", "C", "KS", "AT"}
# structural peptides, immunity and accessory markers
MARKER = {"Mersacidin", "L_biotic_typeA", "Bacteriocin_II", "Bacteriocin_IIc",
          "Bacteriocin_IId", "Bacteriocin_IIi", "BacteriocIIc_cy",
          "Antimicrobial17", "MccV", "EntA_Immun", "PqqD", "Peptidase_S8",
          "PCP/ACP", "KR", "TE", "MT", "KS_C", "A_C", "KS_assoc", "ER/NAD",
          "Aminotran_3", "Radical_SAM", "Gallidermin", "Antimicrobial18",
          "Lactococcin_972", "Lanthi_prec"}
MOBILE = {"Integrase", "Transposase", "Recombinase"}
NRPS = {"A", "C", "KS", "AT", "PCP/ACP", "KR", "TE", "KS_C", "A_C"}


def load_domains():
    dom = collections.defaultdict(lambda: collections.defaultdict(list))
    plen = {}
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            a, prot, pf, lab, af, at, ev, sc, pl = line.rstrip("\n").split("\t")
            dom[a][prot].append((int(af), lab))
            plen[(a, prot)] = int(pl)
    return dom, plen


def load_cds():
    cds = collections.defaultdict(lambda: collections.defaultdict(list))
    rep = {}
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, r, ty, ln, pl, prot, s, e, st, lt, prod = f
            cds[a][r].append((int(s), int(e), st, prot, lt, prod))
            rep[(a, r)] = (ty, int(ln), pl)
    return cds, rep


def main():
    dom, plen = load_domains()
    cds, rep = load_cds()
    cl_rows, cds_rows = [], []
    n = 0
    for a in sorted(cds):
        for r in cds[a]:
            genes = sorted(cds[a][r])
            hits = [(i, g) for i, g in enumerate(genes)
                    if g[3] in dom.get(a, {})]
            if not hits:
                continue
            run = []
            runs = []
            for i, g in hits:
                if run and g[0] - run[-1][1][1] > MAXGAP:
                    runs.append(run); run = []
                run.append((i, g))
            if run:
                runs.append(run)
            for run in runs:
                arch = collections.Counter()
                for _, g in run:
                    for _, lab in dom[a][g[3]]:
                        arch[lab] += 1
                anchors = {d for d in arch if d in ANCHOR}
                # 「至少兩條帶域蛋白」這一關會把**單獨一條巨型合成酶**刷掉：
                # 71 株 S. epidermidis 的 ausA 旁邊 10 kb 內剛好沒有第二條帶域蛋白，
                # 整個座就消失了（2026-09-20 由 draft 基因體的外部 PCR 驗證抓到：
                # 引子在這些「非帶菌者」身上照樣長出產物，因為基因其實在）。
                big_here = any(plen[(a, g[3])] >= 700
                               and len({lab for _, lab in dom[a][g[3]]} & NRPS) >= 2
                               for _, g in run)
                if not anchors or (len(run) < 2 and not big_here):
                    continue
                # a lone fatty-acid-synthesis operon is not a BGC: require a
                # true multi-domain module before calling an NRPS/PKS cluster
                kind = "RiPP"
                if anchors <= {"A", "C", "KS", "AT"}:
                    big = any(plen[(a, g[3])] >= 700 and
                              len({lab for _, lab in dom[a][g[3]]} & NRPS) >= 2
                              for _, g in run)
                    if not big:
                        continue
                    kind = "NRPS/PKS"
                elif anchors & {"A", "C", "KS"}:
                    kind = "hybrid"
                # 界線縮到最外側的核心基因，再各留 FLANK
                core = [g for _, g in run
                        if ({lab for _, lab in dom[a][g[3]]} & ANCHOR)
                        and not PRIMARY_RE.search(g[5] or "")]
                if not core:
                    continue
                clo, chi = min(g[0] for g in core), max(g[1] for g in core)
                inside0 = [g for g in genes
                           if g[1] >= clo - FLANK and g[0] <= chi + FLANK]
                # 報告的島長是核心基因的跨距；兩側 FLANK 只用來收鄰近基因，不計入 span
                lo, hi = clo, chi
                wlo = max(1, min(g[0] for g in inside0))
                whi = max(g[1] for g in inside0)
                arch = collections.Counter()
                for g in inside0:
                    for _, lab in dom[a].get(g[3], []):
                        arch[lab] += 1
                anchors = {d for d in arch if d in ANCHOR}
                ndom = sum(1 for g in inside0 if g[3] in dom.get(a, {}))
                # 單獨一條多模組巨型合成酶（例如 AusA 2,387 aa）本身就是一個座，
                # 不該被「至少兩條帶域蛋白」的規則刷掉（RiPP 操縱子才需要那條）。
                big_nrps = any(plen[(a, g[3])] >= 700 and
                               len({lab for _, lab in dom[a][g[3]]} & NRPS) >= 2
                               for g in inside0 if g[3] in dom.get(a, {}))
                if not anchors or (ndom < 2 and not big_nrps):
                    continue
                kind = "RiPP"
                if anchors <= {"A", "C", "KS", "AT"}:
                    if not any(plen[(a, g[3])] >= 700 and
                               len({lab for _, lab in dom[a][g[3]]} & NRPS) >= 2
                               for g in inside0 if g[3] in dom.get(a, {})):
                        continue
                    kind = "NRPS/PKS"
                elif anchors & {"A", "C", "KS"}:
                    kind = "hybrid"
                ty, rlen, pl = rep[(a, r)]
                n += 1
                cid = f"C{n:05d}"
                inside = inside0
                mob = sum(1 for g in inside if g[3] in dom.get(a, {})
                          and {lab for _, lab in dom[a][g[3]]} & MOBILE)
                cl_rows.append(dict(
                    cluster=cid, assembly=a, replicon=r, replicon_type=ty,
                    plasmid=pl, replicon_len=rlen, start=lo, end=hi,
                    span=hi - lo + 1, window_start=wlo, window_end=whi,
                    n_cds=len(inside), kind=kind,
                    anchors=",".join(sorted(anchors)),
                    arch=";".join(f"{d}x{c}" for d, c in sorted(arch.items())),
                    mobile_cds=mob))
                for s, e, st, prot, lt, prod in inside:
                    ds = [lab for _, lab in sorted(dom[a].get(prot, []))]
                    cds_rows.append(dict(
                        cluster=cid, assembly=a, replicon=r, start=s, end=e,
                        strand=st, protein=prot, locus_tag=lt, product=prod,
                        aa=plen.get((a, prot), ""), domains="-".join(ds)))
    for name, rows in (("clusters.tsv", cl_rows), ("cluster_cds.tsv", cds_rows)):
        with open(os.path.join(OUT, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(name, len(rows), flush=True)


if __name__ == "__main__":
    main()
