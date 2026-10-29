#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""每個 RiPP 家族的前驅肽：切割位、核心肽序列與理論質量（MALDI-TOF 可查的 m/z）。

臨床實驗室手上有 MALDI-TOF（菌種鑑定用的那台），所以除了 PCR 之外，
再給一組「該看哪個 m/z」的數字。做法全部寫死在這裡，可複核：

1. 前驅肽＝基因簇內帶前驅肽域（Gallidermin／L_biotic_typeA／Lanthi_prec／
   Bacteriocin_*）或長度 ≤ 90 aa 且在基因簇核心 ±5 kb 內的蛋白。
2. 切割位：class II（雙甘胺酸型，LanM／C39 系統）找最後一個 GG/GA/GS；
   class I（LanB/LanC 系統，gallidermin／epidermin 型）文獻上的成熟位點在
   前導肽末端的 -1 位 R/K 之後或 GA 之後，這裡兩種都試，取較長的核心肽。
3. 質量：單一同位素與平均質量各算一次；羊毛硫肽每形成一個硫醚環少一個水
   （-18.0106 / -18.0153），所以同時列出 0..n 次脫水的質量，n＝核心肽裡的
   Ser+Thr 數。**不宣稱哪一個是真的**，只說「落在這個梯子上」。
4. 對照：用已知的 gallidermin（2164 Da）與 epidermin（2164 Da）核對梯子上
   對得到哪一格；對不上就要在稿件裡講清楚。
"""
import collections, csv, json, os, re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FAA = os.path.join(HERE, "data", "faa")
MAXPREC = 90
PREC_DOM = {"Gallidermin", "L_biotic_typeA", "Lanthi_prec", "Mersacidin",
            "Bacteriocin_II", "Bacteriocin_IIc", "Bacteriocin_IId",
            "Bacteriocin_IIi", "BacteriocIIc_cy", "Antimicrobial17",
            "Antimicrobial18", "MccV", "Lactococcin_972"}

# 單一同位素 / 平均 殘基質量（Da）
MONO = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841,
            T=101.04768, C=103.00919, L=113.08406, I=113.08406, N=114.04293,
            D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049,
            H=137.05891, F=147.06841, R=156.10111, Y=163.06333, W=186.07931)
AVG = dict(G=57.0519, A=71.0788, S=87.0782, P=97.1167, V=99.1326, T=101.1051,
           C=103.1388, L=113.1594, I=113.1594, N=114.1038, D=115.0886,
           Q=128.1307, K=128.1741, E=129.1155, M=131.1926, H=137.1411,
           F=147.1766, R=156.1875, Y=163.1760, W=186.2132)
H2O_MONO, H2O_AVG = 18.010565, 18.01528
PROTON = 1.007276
# epidermin／gallidermin 型的 C 端 AviCys：EpiD 把末端 Cys 氧化脫羧，
# 質量少掉 CO2 與 2H。核心肽以 Cys 結尾時另外列一條「脫羧」梯子。
DECARB_AVG, DECARB_MONO = 46.026, 45.9929
# 已發表的成熟質量，用來標出梯子上對得到的那一格。
# **質量型態很重要**：gallidermin／epidermin 文獻報的 2164 Da 是**單一同位素質量**
# （nano-ESI-MS 的 [M+H]+ 落在 2165），不是平均質量。2.2 kDa 的胜肽 average 比
# monoisotopic 大約大 1.7 Da，拿 average 去對會憑空多出 1.6 Da 的差，看起來像「差不多」
# 其實是比錯量。兩種質量都算、都報，並在 published_type 欄寫明對的是哪一種。
PUBLISHED = {"IASKFLCTPGCAKTGSFNSYCC": ("gallidermin", 2164.0, "monoisotopic"),
             "IASKFICTPGCAKTGSFNSYCC": ("epidermin", 2164.0, "monoisotopic")}


def mass(seq, table, water):
    return sum(table[c] for c in seq if c in table) + water


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


# 已驗證的切割規則：取「最後一個」G[GAS] 或 P[QRK] 之後開始，核心肽長度 12–45 aa。
# 這條規則是照 UniProt 的 PEPTIDE 註記反推的，下面 CHECKS 會在每次執行時重跑，
# 四個已知肽只要有一個對不上就整支停掉（規則改壞了要馬上知道）。
# [AVE]R 是 type-AI 前導肽的另一種結尾（capitis 的 44 aa 前驅肽 …DVKVNAR↓）；
# 加進來之後上面五個已知肽的切點都沒有跑掉（self_check 每次會驗）。
CLEAVE_RE = re.compile(r"G[GAS]|P[QRK]|[AVE]R")
CHECKS = [
    # (UniProt accession, 前驅肽全長序列, 文獻／UniProt 註記的成熟核心肽)
    ("P21838", "MEAVKEKNELFDLDVKVNAKESNDSGAEPRIASKFLCTPGCAKTGSFNSYCC",
     "IASKFLCTPGCAKTGSFNSYCC"),                      # gallidermin
    ("P08136", "MEAVKEKNDLFNLDVKVNAKESNDSGAEPRIASKFICTPGCAKTGSFNSYCC",
     "IASKFICTPGCAKTGSFNSYCC"),                      # epidermin
    ("Q9KWM4", "MENSKVMKDIEVANLLEEVQEDELNEVLGAKKKSGVIPTVSHDCHMNSFQFVFTCCS",
     "KKKSGVIPTVSHDCHMNSFQFVFTCCS"),                 # nukacin
    ("P19578", "MKNNKNLFDLEIKKETSQNTDELEPQTAGPAIRASVKQCQKTLKATRLFTVSCKGKNGCK",
     "TAGPAIRASVKQCQKTLKATRLFTVSCKGKNGCK"),          # Pep5
    ("P86047", "MKKELFDLNLNKDIEAQKSDLNPQSASIVKTTIKASKKLCRGFTLTCGCHFTGKK",
     "SASIVKTTIKASKKLCRGFTLTCGCHFTGKK"),             # epilancin 15X
]


def cleave(seq):
    """回傳 (核心肽, 切割位描述)；找不到合格切點就回 ('', 'none')。"""
    best = ("", "none")
    for m in CLEAVE_RE.finditer(seq):
        core = seq[m.end():]
        if 12 <= len(core) <= 45:
            best = (core, f"{m.group(0)} at {m.end()}")
    return best


def self_check():
    bad = []
    for acc, prec, core in CHECKS:
        got, how = cleave(prec)
        if got != core:
            bad.append(f"{acc}: 預期 {core}，實得 {got} ({how})")
    if bad:
        raise SystemExit("切割規則自我檢查失敗：\n  " + "\n  ".join(bad))
    print(f"切割規則自我檢查通過（{len(CHECKS)} 個已知肽）", flush=True)


def main():
    self_check()
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cds = list(csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"))
    seqs = {}
    need = collections.defaultdict(list)
    for r in cds:
        fam = fam_of.get(r["cluster"], "")
        doms = set((r["domains"] or "").split("-"))
        aa = int(r["aa"] or 0)
        if (doms & PREC_DOM) or (0 < aa <= MAXPREC):
            need[r["assembly"]].append(r)
    for a, rs in need.items():
        want = {r["protein"] for r in rs}
        for n, s in read_fasta(os.path.join(FAA, a + ".faa")):
            if n in want:
                seqs[(a, n)] = s

    rows = []
    for a, rs in sorted(need.items()):
        for r in rs:
            s = seqs.get((a, r["protein"]), "")
            if not s or len(s) > MAXPREC:
                continue
            core, how = cleave(s)
            if not core:
                continue
            st = sum(core.count(c) for c in "ST")
            mono = mass(core, MONO, H2O_MONO)
            avg = mass(core, AVG, H2O_AVG)
            rows.append(dict(
                family=fam_of.get(r["cluster"], ""), cluster=r["cluster"],
                assembly=a, protein=r["protein"], product=r["product"],
                precursor_aa=len(s), core_aa=len(core), cleavage=how,
                core_seq=core, n_ser_thr=st,
                mono_core=round(mono, 3), avg_core=round(avg, 2),
                mono_mh_plus=round(mono + PROTON, 3),
                ladder_avg=";".join(
                    f"{k}:{avg - k * H2O_AVG:.1f}" for k in range(0, st + 1)),
                ladder_mono=";".join(
                    f"{k}:{mono - k * H2O_MONO:.2f}" for k in range(0, st + 1)),
                ladder_avg_decarb="" if not core.endswith("C") else ";".join(
                    f"{k}:{avg - k * H2O_AVG - DECARB_AVG:.1f}"
                    for k in range(0, st + 1)),
                ladder_mono_decarb="" if not core.endswith("C") else ";".join(
                    f"{k}:{mono - k * H2O_MONO - DECARB_MONO:.2f}"
                    for k in range(0, st + 1)),
                published=PUBLISHED.get(core, ("", "", ""))[0],
                published_mass=PUBLISHED.get(core, ("", "", ""))[1],
                published_type=PUBLISHED.get(core, ("", "", ""))[2],
                # 比對一律同型：已發表的是單一同位素，就拿單一同位素的梯子去比
                best_ladder_gap=("" if core not in PUBLISHED else round(min(
                    abs(mono - k * H2O_MONO - (DECARB_MONO if core.endswith("C") else 0)
                        - PUBLISHED[core][1]) for k in range(0, st + 1)), 2)),
                best_ladder_gap_if_average=("" if core not in PUBLISHED else round(min(
                    abs(avg - k * H2O_AVG - (DECARB_AVG if core.endswith("C") else 0)
                        - PUBLISHED[core][1]) for k in range(0, st + 1)), 2))))
    p = os.path.join(OUT, "precursor_masses.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("wrote", p, len(rows), "precursors", flush=True)
    byfam = collections.defaultdict(list)
    for r in rows:
        byfam[r["family"]].append(r)
    for f in sorted(byfam):
        rs = byfam[f]
        uniq = collections.Counter(r["core_seq"] for r in rs)
        seq, n = uniq.most_common(1)[0]
        ex = next(r for r in rs if r["core_seq"] == seq)
        print(f"  {f}: {len(rs)} precursors, {len(uniq)} distinct cores; "
              f"modal {seq} ({n}x) avg {ex['avg_core']:.1f} Da, "
              f"ladder {ex['ladder_avg']}", flush=True)


if __name__ == "__main__":
    main()
