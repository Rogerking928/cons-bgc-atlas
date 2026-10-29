#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一條呼叫規則可以造出一個「顯著關聯」（2026-09-20，由外部驗證抓到）。

原本的基因簇呼叫要求「視窗內至少兩條帶域蛋白」。S. epidermidis 的 ausA 是單獨一條
2,387 aa 的巨型合成酶，旁邊 10 kb 內常常沒有第二條帶域蛋白，於是 71 株帶著基因卻被
判成不帶菌——而那 71 株剛好偏在某些譜系（ST23 全部落在裡面），結果做出
「F001 與 erm(A) 強烈負相關（OR 0.045, P=8.8e-13）」這種看起來很漂亮的結論。

draft 基因體的外部 PCR 先發現了破綻：引子在這些「非帶菌者」身上照樣長出產物。

這支把三種帶菌定義並排算出來，給 Figure 3 用：
  strict  舊規則（視窗內至少兩條帶域蛋白）
  fixed   現行規則（單獨一條多模組巨型合成酶也算）
  gene    只要 RefSeq 註記裡有那個合成酶蛋白就算
輸出 out/calling_rule.json。
"""
import collections, csv, gzip, json, os

from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
CALL = 0.5
MAXGAP = 10000
NRPS = {"A", "C", "KS", "AT", "PCP/ACP", "KR", "TE", "KS_C", "A_C"}
GENE = "AusA"
FAM = "F001"
SPECIES = "epidermidis"
DETS = ["erm(A)", "ant(9)-Ia", "mecA", "blaZ"]


def main():
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    amr = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv")))}
    mlst = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "mlst.csv")))}
    pool = [a for a, r in mat.items() if r["species"] == SPECIES]

    dom = collections.defaultdict(dict)
    plen = {}
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            a, prot, pf, lab, af, at, ev, sc, pl = line.rstrip("\n").split("\t")
            dom[a].setdefault(prot, set()).add(lab)
            plen[(a, prot)] = int(pl)

    genes = collections.defaultdict(lambda: collections.defaultdict(list))
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genes[f[0]][f[1]].append((int(f[6]), int(f[7]), f[5], f[10]))

    strict, fixed, gene = set(), set(), set()
    for a in pool:
        for r, gg in genes[a].items():
            gg.sort()
            hits = [g for g in gg if g[2] in dom.get(a, {})]
            for g in gg:
                if GENE in (g[3] or ""):
                    gene.add(a)
            # 以 ausA 為核心的 run
            for i, g in enumerate(hits):
                if GENE not in (g[3] or ""):
                    continue
                run = [g]
                j = i - 1
                while j >= 0 and run[0][0] - hits[j][1] <= MAXGAP:
                    run.insert(0, hits[j]); j -= 1
                j = i + 1
                while j < len(hits) and hits[j][0] - run[-1][1] <= MAXGAP:
                    run.append(hits[j]); j += 1
                big = any(plen[(a, x[2])] >= 700 and len(dom[a][x[2]] & NRPS) >= 2
                          for x in run)
                if len(run) >= 2:
                    strict.add(a)
                if len(run) >= 2 or big:
                    fixed.add(a)

    defs = {"strict (>=2 domain-carrying CDS)": strict,
            "fixed (a lone megasynthetase counts)": fixed,
            "gene present (RefSeq annotation)": gene}

    res = {"n_genomes": len(pool), "carriage": {}, "association": {}, "by_ST": {}}
    for name, s in defs.items():
        res["carriage"][name] = len(s & set(pool))
        row = {}
        for det in DETS:
            def has(x):
                return float(amr[x].get(det, 0) or 0) >= CALL
            a_ = sum(1 for x in pool if x in s and has(x))
            b_ = sum(1 for x in pool if x in s and not has(x))
            c_ = sum(1 for x in pool if x not in s and has(x))
            d_ = sum(1 for x in pool if x not in s and not has(x))
            orr, p = stats.fisher_exact([[a_, b_], [c_, d_]])
            ha = ((a_ + .5) * (d_ + .5)) / ((b_ + .5) * (c_ + .5))
            # Woolf 區間（四格各加 0.5）：讀者要比較 0.05 和 0.17，沒有區間比不了
            import math
            A, B_, C_, D_ = a_ + .5, b_ + .5, c_ + .5, d_ + .5
            se = math.sqrt(1 / A + 1 / B_ + 1 / C_ + 1 / D_)
            lo = math.exp(math.log((A * D_) / (B_ * C_)) - 1.96 * se)
            hi = math.exp(math.log((A * D_) / (B_ * C_)) + 1.96 * se)
            row[det] = {"a": a_, "b": b_, "c": c_, "d": d_,
                        "odds": float(orr), "odds_haldane": float(ha),
                        "ci": [float(lo), float(hi)], "p": float(p)}
        res["association"][name] = row
        sts = collections.Counter(mlst.get(x, {}).get("ST", "") for x in pool)
        res["by_ST"][name] = {
            f"ST{st}": {"n": n, "carriers": sum(1 for x in pool
                                                if mlst.get(x, {}).get("ST") == st and x in s)}
            for st, n in sts.most_common() if st and n >= 7}
        print(f"  {name:40s} carriers {len(s & set(pool)):3d}/{len(pool)}  "
              f"erm(A) OR {row['erm(A)']['odds']:.3f} (P {row['erm(A)']['p']:.2g})",
              flush=True)

    json.dump(res, open(os.path.join(OUT, "calling_rule.json"), "w"), indent=1)
    print("-> out/calling_rule.json", flush=True)


if __name__ == "__main__":
    main()
