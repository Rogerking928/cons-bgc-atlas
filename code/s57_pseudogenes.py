#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pseudogenised machinery: what a proteome-based atlas cannot see.

這支腳本是 antiSMASH 對照（s54／s55）追出來的結果，不是另起的題目。
對照發現 antiSMASH 在圖譜判為「非帶菌」的少數基因體上照樣呼叫 NRPS 區段；
追到那些區段，核心基因是 RefSeq 標成 **pseudogene** 的 `ausA`。原因是結構性的：

  圖譜從 protein.faa 出發（§2.2），而 RefSeq 不為假基因產生蛋白質條目，
  所以假基因化的合成酶對圖譜而言不存在；antiSMASH 讀的是 GFF 與核酸序列，
  所以它看得到。

兩者都沒有錯，是**偵測基質不同**。這支腳本把這件事量化：掃全 panel 的 GFF，
找出落在各家族標記基因名下的假基因，並與圖譜的帶菌判定交叉比對，
算出每個家族「基因座存在」與「合成酶完整」兩個數字。

標記基因名不是手打的清單：先從圖譜每個家族的成員蛋白的 RefSeq `gene=` 名稱
歸納出該家族用到的基因名，再拿這組名字去找假基因，所以新家族不必改程式。

Output: out/pseudogenes.csv      每個命中的假基因一列
        out/pseudogene_summary.json
"""
import collections, csv, glob, gzip, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s6_call_clusters import ANCHOR, MOBILE  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
GFF = os.path.join(HERE, "data", "gff")

ATTR_NAME = re.compile(r"(?:^|;)Name=([^;]+)")
ATTR_GENE = re.compile(r"(?:^|;)gene=([^;]+)")
ATTR_NOTE = re.compile(r"(?:^|;)Note=([^;]+)")
ATTR_LOCUS = re.compile(r"(?:^|;)locus_tag=([^;]+)")


def attr(pat, s, default=""):
    m = pat.search(s)
    return m.group(1) if m else default


def family_gene_names():
    """每個家族用到的 RefSeq 基因名，從圖譜自己的成員蛋白歸納。

    out/cluster_cds.tsv 是每個簇的成員 CDS；它帶 product 但不一定帶 gene 名，
    所以基因名從 GFF 取（以 protein id 對上）。為了只讀一次 GFF，這裡先把
    需要的 protein id 收齊。
    """
    fam = json.load(open(os.path.join(OUT, "cluster_family.json")))
    want = collections.defaultdict(set)          # assembly -> {protein}
    prot_fam = {}                                # (assembly, protein) -> family
    with open(os.path.join(OUT, "cluster_cds.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            f = fam.get(r["cluster"])
            if not f:
                continue
            # 只取**帶 anchor 域**的成員。cluster_cds.tsv 也收簇旁 5 kb 內的
            # 鄰居（會混進 arsB、blaR1 這類名字），而可動性基因雖然帶域，
            # 卻不是生合成標記（否則 ccrA／ccrB 的假基因會被當成基因座存在）。
            labels = {d for part in (r.get("domains") or "").split(";")
                      for d in part.split("-") if d}
            if not (labels & ANCHOR) or labels <= MOBILE:
                continue
            want[r["assembly"]].add(r["protein"])
            prot_fam[(r["assembly"], r["protein"])] = f
    names = collections.defaultdict(collections.Counter)
    # 家族 -> 它出現過的物種。家族是物種專一的，所以基因名要連物種一起對，
    # 否則 ausA（四個家族都有這個名字）會被算進每一個家族。
    for acc, prots in want.items():
        path = os.path.join(GFF, acc + ".gff.gz")
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS":
                    continue
                pid = attr(ATTR_NAME, f[8])
                if pid in prots:
                    g = attr(ATTR_GENE, f[8])
                    if g:
                        names[prot_fam[(acc, pid)]][g] += 1
    return {k: dict(v) for k, v in names.items()}


def family_species(mat, md):
    """每個家族實際出現的物種集合，取自株 x 家族矩陣。"""
    sp = collections.defaultdict(set)
    fams = [c for c in (next(iter(mat.values())).keys() if mat else [])
            if c.startswith("F") and c[1:].isdigit()]
    for acc, row in mat.items():
        for f in fams:
            if row.get(f):
                sp[f].add(md[acc]["species"])
    return sp


def main():
    fam_names = family_gene_names()
    md = {r["assembly"]: r for r in
          csv.DictReader(open(os.path.join(OUT, "metadata.csv")))}
    mat = {r["assembly"]: r for r in
           csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    fam_sp = family_species(mat, md)

    gene_fam = {}
    for f, counts in fam_names.items():
        for g in counts:
            gene_fam.setdefault(g, set()).add(f)
    if not gene_fam:
        sys.exit("no gene names recovered from the atlas - run s11/s12 first")

    rows = []
    for acc in sorted(md):
        path = os.path.join(GFF, acc + ".gff.gz")
        if not os.path.exists(path):
            continue
        # 假基因化的原因寫在子 CDS 行的 Note，不在 pseudogene 行上，
        # 所以先按 locus_tag 把子行的 Note 收起來。
        with gzip.open(path, "rt") as fh:
            notes = {}
            pseudo = []
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[2] == "CDS":
                    n = attr(ATTR_NOTE, f[8])
                    if n:
                        notes[attr(ATTR_LOCUS, f[8])] = (
                            n.replace("%3B", ";").replace("%2C", ","))
                elif f[2] == "pseudogene":
                    pseudo.append(f)
        for f in pseudo:
                g = attr(ATTR_NAME, f[8]) or attr(ATTR_GENE, f[8])
                if g not in gene_fam:
                    continue
                note = notes.get(attr(ATTR_LOCUS, f[8]), "")
                for family in sorted(gene_fam[g]):
                    # 家族專一於某些物種，別的物種帶同名基因不算這個家族的座
                    if md[acc]["species"] not in fam_sp.get(family, set()):
                        continue
                    carrier = mat.get(acc, {}).get(family, "")
                    rows.append([acc, md[acc]["species"], md[acc]["strain"],
                                 family, g, f[0], int(f[3]), int(f[4]),
                                 int(f[4]) - int(f[3]) + 1, f[6], note,
                                 attr(ATTR_LOCUS, f[8]),
                                 int(bool(carrier))])

    with open(os.path.join(OUT, "pseudogenes.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["assembly", "species", "strain", "family", "gene",
                    "replicon", "start", "end", "span", "strand", "note",
                    "locus_tag", "atlas_calls_carrier"])
        w.writerows(rows)

    # ------------------------------------------------------------- 摘要
    # 每個家族：圖譜帶菌數、帶假基因但圖譜非帶菌的數（= 基因座在、合成酶不完整）
    per_fam = collections.defaultdict(lambda: {"pseudo": 0, "pseudo_noncarrier": 0,
                                               "genes": collections.Counter(),
                                               "species": collections.Counter()})
    for r in rows:
        d = per_fam[r[3]]
        d["pseudo"] += 1
        d["genes"][r[4]] += 1
        d["species"][r[1]] += 1
        if not r[12]:
            d["pseudo_noncarrier"] += 1
    summary = {}
    for family, d in sorted(per_fam.items()):
        carriers = sum(1 for a in mat if mat[a].get(family))
        # 同一基因體可能有同家族的多個假基因，非帶菌者要按基因體去重
        nc = len({r[0] for r in rows if r[3] == family and not r[12]})
        summary[family] = {
            "atlas_carriers": carriers,
            "pseudogenes": d["pseudo"],
            "genomes_with_pseudogene_but_not_called_carrier": nc,
            "locus_present_including_pseudogenised": carriers + nc,
            "genes": dict(d["genes"].most_common()),
            "species": dict(d["species"].most_common()),
        }
    notes = collections.Counter(
        ("frameshifted" if "frameshift" in r[10] else
         "internal stop" if "internal stop" in r[10] else
         "partial" if "incomplete" in r[10] or "partial" in r[10] else
         "other") for r in rows)
    out = {"families": summary,
           "family_gene_names": fam_names,
           "pseudogene_rows": len(rows),
           "genomes_affected": len({r[0] for r in rows}),
           "reasons": dict(notes.most_common())}
    with open(os.path.join(OUT, "pseudogene_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"{len(rows)} pseudogenised marker genes in "
          f"{len({r[0] for r in rows})} genomes")
    for family, s in summary.items():
        print(f"  {family}: atlas carriers {s['atlas_carriers']}, "
              f"pseudogene in a non-carrier "
              f"{s['genomes_with_pseudogene_but_not_called_carrier']}, "
              f"locus present {s['locus_present_including_pseudogenised']} "
              f"({', '.join(s['genes'])})")
    print(f"  reasons: {dict(notes)}")
    print("-> out/pseudogenes.csv, out/pseudogene_summary.json")


if __name__ == "__main__":
    main()
