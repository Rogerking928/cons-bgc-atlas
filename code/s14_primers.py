#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A PCR typing scheme for the main cluster families, designed and validated
entirely in silico on the 834 complete genomes.

For each family: take the marker gene carried by (almost) every member, keep
the 20-mers that are exact in at least MINPREV of the carriers' copies, pair
them into amplicons of 150-900 bp that satisfy simple primer rules, and then
run the pair against every genome in the panel (both strands, exact match) to
get sensitivity and specificity.
"""
import collections, csv, gzip, json, math, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FNA = os.path.join(HERE, "data", "fna")
PLEN_RANGE, MINPREV = (18, 25), 0.90
AMP = (120, 1200)
# 每個家族挑一個幾乎人人都有的標記基因（家族編號見 out/families.tsv）
TARGETS = {"F001-ausA": ("F001", "A"),        # epidermidis aureusimine
           "F002-ausA": ("F002", "A"),        # capitis aureusimine
           "F003-nrps": ("F003", "A"),        # lugdunensis 第二套 NRPS
           "F004-lanB": ("F004", "LanB_N"),   # capitis gallidermin 脫水酶
           "F005-ycaO": ("F005", "YcaO"),     # lugdunensis hybrid
           "F006-lugC": ("F006", "C"),        # lugdunin（A domain 那條在核苷酸層
           "F006-lugTE": ("F006", "TE"),      # 上對不齊，改用 C domain 那條）
           "F007-c39": ("F007", "C39"),
           "F009-lanM": ("F009", "LanM_cyc"), # lacticin 481 型（質體）
           "F011-lanC": ("F011", "LanC")}

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


# SantaLucia (1998) unified nearest-neighbour parameters, dH kcal/mol, dS cal/mol/K
NN = {"AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4),
      "TA": (-7.2, -21.3), "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
      "GT": (-8.4, -22.4), "AC": (-8.4, -22.4), "CT": (-7.8, -21.0),
      "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
      "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9),
      "CC": (-8.0, -19.9)}
NA, CT_PRIMER = 0.05, 0.25e-6          # 50 mM Na+, 0.25 uM primer


def tm(p):
    """Nearest-neighbour melting temperature, salt-corrected (SantaLucia 1998)."""
    dh = ds = 0.0
    for i in range(len(p) - 1):
        d = NN.get(p[i:i + 2])
        if d is None:
            return 0.0
        dh += d[0]; ds += d[1]
    # initiation
    dh += 0.1 if p[0] in "GC" else 2.3
    ds += -2.8 if p[0] in "GC" else 4.1
    dh += 0.1 if p[-1] in "GC" else 2.3
    ds += -2.8 if p[-1] in "GC" else 4.1
    ds += 0.368 * (len(p) - 1) * math.log(NA)
    return dh * 1000 / (ds + 1.987 * math.log(CT_PRIMER / 4)) - 273.15


def ok_primer(p):
    gc = sum(c in "GC" for c in p) / len(p)
    if not 0.35 <= gc <= 0.65:
        return False
    if not 55.0 <= tm(p) <= 65.0:
        return False
    if any(b * 4 in p for b in "ACGT"):
        return False
    if p[-1] not in "GC":
        return False
    if sum(c in "GC" for c in p[-5:]) > 3:
        return False
    return True


def load_seqs(acc):
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
    return seqs


def find_all(s, sub):
    out, i = [], s.find(sub)
    while i >= 0:
        out.append(i)
        i = s.find(sub, i + 1)
    return out


def amplify(seqs, fwd, rev):
    """Exact-match in-silico PCR on both strands; returns list of product sizes."""
    hits = []
    rrc = rc(rev)
    frc = rc(fwd)
    for s in seqs.values():
        for i in find_all(s, fwd):
            for j in find_all(s, rrc):
                n = j + len(rev) - i
                if AMP[0] <= n <= AMP[1]:
                    hits.append(n)
        for i in find_all(s, rev):
            for j in find_all(s, frc):
                n = j + len(fwd) - i
                if AMP[0] <= n <= AMP[1]:
                    hits.append(n)
    return hits


def main():
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cds = list(csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"))
    mat = list(csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv"))))
    carriers = {lab: {r["assembly"] for r in mat if r[fam]}
                for lab, (fam, _) in TARGETS.items()}
    have = {os.path.basename(p)[:-7] for p in os.listdir(FNA)} if os.path.isdir(FNA) else set()
    have = {a[:-7] if a.endswith(".fna.gz") else a for a in os.listdir(FNA)}
    print(f"{len(have)} genomes with nucleotide sequence", flush=True)

    rows, cand = [], []
    for lab, (fam, dom) in TARGETS.items():
        # marker CDS of every member cluster
        genes = collections.defaultdict(list)
        for r in cds:
            if fam_of.get(r["cluster"]) != fam:
                continue
            if dom in (r["domains"] or "").split("-"):
                genes[r["assembly"]].append(r)
        # 一個座裡常有好幾條同域蛋白（lug 座就有 LugA/LugC/LugD 三條 A domain）。
        # 取「位置上的第一條」會在不同基因體挑到不同基因，序列自然對不起來，
        # 保守 20-mer 就一個都找不到。改成先選家族內最常見的產物名，再按名字挑。
        pnames = collections.Counter(r["product"] for v in genes.values() for r in v)
        modal = pnames.most_common(1)[0][0] if pnames else ""
        # 只留「產物名就是眾數那個基因」的基因體：退而取位置上的第一條會在不同
        # 基因體挑到不同的旁系基因（lugdunin 座有 LugA/LugC/LugD 三條 C domain），
        # 保守 20-mer 就一個都找不到。敏感度照樣在全 panel 上實測，不受影響。
        targets = {}
        for a, v in genes.items():
            if a not in have:
                continue
            same = [r for r in v if r["product"] == modal]
            if same:
                targets[a] = same[0]
        if len(targets) < 5:
            targets = {a: v[0] for a, v in genes.items() if a in have}
        print(f"{lab}: marker {dom} in {len(genes)}/{len(carriers[lab])} carriers, "
              f"{len(targets)} with sequence", flush=True)
        if len(targets) < 5 or len(genes) < 5:
            continue
        seqs = {}
        for a, r in targets.items():
            s = load_seqs(a).get(r["replicon"], "")
            sub = s[int(r["start"]) - 1:int(r["end"])]
            seqs[a] = sub if r["strand"] == "+" else rc(sub)
        ref = max(seqs, key=lambda a: len(seqs[a]))
        refseq = seqs[ref]
        others = [v for k, v in seqs.items() if k != ref]

        def conserved(minprev):
            need = minprev * len(others)
            out = []
            for i in range(0, len(refseq) - PLEN_RANGE[0] + 1):
                for L in range(*PLEN_RANGE):
                    p = refseq[i:i + L]
                    if len(p) < L or "N" in p or not ok_primer(p):
                        continue
                    n = sum(1 for o in others if p in o)
                    if n >= need:
                        out.append((i, p, n / max(1, len(others))))
                        break
            return out

        # 保守度門檻降級：有些家族（lugdunin 的 lugA）在 90% 這一關找不到任何
        # 共通 20-mer，硬卡著就等於沒有引子。降到 80% 再試一次，實際保守度
        # 逐組記在 conserved_frac 欄，敏感度照樣在全 panel 上實測。
        used_prev = MINPREV
        cons = conserved(MINPREV)
        if len(cons) < 4:
            used_prev = 0.80
            cons = conserved(used_prev)
        print(f"  {len(cons)} conserved candidate primers "
              f"(prevalence threshold {used_prev:.2f})", flush=True)
        pairs = []
        for a in range(len(cons)):
            ia, pa, fa = cons[a]
            for b in range(len(cons) - 1, a, -1):
                ib, pb, fb = cons[b]
                size = ib + len(pb) - ia
                if not AMP[0] <= size <= AMP[1]:
                    continue
                rev = rc(pb)
                if not ok_primer(rev) or abs(tm(pa) - tm(rev)) > 3.0:
                    continue
                pairs.append((min(fa, fb), size, pa, rev, ia, ib))
                break
        pairs.sort(key=lambda t: (-t[0], abs(t[1] - 400)))
        chosen = []
        for score, size, fwd, rev, ia, ib in pairs:
            # keep alternatives genuinely apart, not the same primer shifted 1 nt
            if any(abs(ia - x) < 60 for x in chosen):
                continue
            chosen.append(ia)
            if len(chosen) > 2:
                break
            cand.append(dict(target=lab, family=fam, marker=dom, forward=fwd, reverse=rev,
                             product_bp=size, tm_f=round(tm(fwd), 1),
                             tm_r=round(tm(rev), 1), conserved_frac=round(score, 4)))

    # one pass over the panel, testing every candidate pair on every genome
    print(f"validating {len(cand)} candidate pairs on {len(have)} genomes", flush=True)
    tp = collections.Counter(); fp = collections.Counter()
    sizes = collections.defaultdict(collections.Counter)
    disc = collections.defaultdict(list)
    for n, a in enumerate(sorted(have), 1):
        seqs = load_seqs(a)
        for k, c in enumerate(cand):
            hits = amplify(seqs, c["forward"], c["reverse"])
            carrier = a in carriers[c["target"]]
            if hits:
                sizes[k][hits[0]] += 1
                if carrier:
                    tp[k] += 1
                else:
                    fp[k] += 1
                    disc[k].append(("false positive", a))
            elif carrier:
                disc[k].append(("false negative", a))
        if n % 100 == 0:
            print(f"  {n}/{len(have)}", flush=True)
    for k, c in enumerate(cand):
        npos = len(carriers[c["target"]] & have)
        nneg = len(have) - npos
        c.update(carriers_tested=npos, amplified_carriers=tp[k],
                 sensitivity=round(tp[k] / npos, 4) if npos else "",
                 non_carriers_tested=nneg, amplified_non_carriers=fp[k],
                 specificity=round(1 - fp[k] / nneg, 4) if nneg else "",
                 product_sizes=";".join(f"{x}:{y}" for x, y in sizes[k].most_common(3)),
                 discordant=";".join(f"{w}:{a}" for w, a in disc[k]))
        print("   ", c, flush=True)
    rows = cand
    with open(os.path.join(OUT, "primers.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("wrote out/primers.csv", flush=True)


if __name__ == "__main__":
    main()
