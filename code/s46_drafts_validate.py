#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""外部驗證：把 PCR 方案拿到沒有進過設計集的基因體上重跑（2026-09-20 審閱意見）。

正文的敏感度／特異度是在設計引子的那 554 株完整基因體上量的，必然樂觀。這支在
兩批外部資料上重測：

  * 同 12 種 CoNS 的 draft 組裝（scaffold／contig），完全沒進過設計集；
  * 非目標物種（S. aureus、Cutibacterium、Enterococcus、E. coli、Streptococcus、
    Micrococcus），測的是「不該長東西的地方有沒有長」。

外部基因體上誰是帶菌者，用蛋白質層的 10-mer containment 對該家族的標記蛋白認定
（draft 的基因簇會被 contig 切斷，用基因簇呼叫會低估）。
輸出 out/external_pcr.csv、out/external_pcr.json。
"""
import collections, csv, json, math, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
EXT = os.path.join(HERE, "data", "ext")
FAA = os.path.join(HERE, "data", "faa")
K = 10
CALL = 0.5
AMP = (120, 1200)
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


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


def wilson(k, n):
    if n == 0:
        return (float("nan"),) * 2
    p, z = k / n, 1.96
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def best_primers():
    best = {}
    for r in csv.DictReader(open(os.path.join(OUT, "primers.csv"))):
        k = r["target"]
        sc = (float(r["sensitivity"]), float(r["specificity"]), -int(r["product_bp"]))
        if k not in best or sc > best[k][0]:
            best[k] = (sc, r)
    return [r for _, r in sorted(best.values(), key=lambda v: v[1]["target"])]


def marker_kmers():
    """每個目標家族的標記蛋白（panel 裡的每一份複本）的 10-mer 集合。"""
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cds = list(csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"))
    want = collections.defaultdict(lambda: collections.defaultdict(set))
    for r in cds:
        fam = fam_of.get(r["cluster"])
        if fam:
            want[fam][r["assembly"]].add((r["protein"], r["domains"] or ""))
    return want


def main():
    meta = json.load(open(os.path.join(EXT, "meta.json")))
    pri = best_primers()
    targets = {r["target"]: (r["family"], r["marker"], r) for r in pri}
    want = marker_kmers()

    # panel 裡每個家族的標記蛋白序列
    refs = collections.defaultdict(list)
    for tgt, (fam, marker, _) in targets.items():
        seen = set()
        for acc, prots in want[fam].items():
            names_ = {p for p, doms in prots if marker in doms.split("-")}
            if not names_:
                continue
            for n, s in read_fasta(os.path.join(FAA, acc + ".faa")):
                if n in names_ and s not in seen:
                    seen.add(s)
                    refs[tgt].append(s)
            if len(refs[tgt]) >= 25:
                break
        print(f"{tgt:12s} family {fam:5s} marker {marker:10s} "
              f"{len(refs[tgt])} reference proteins", flush=True)

    rows = []
    accs = sorted(a[:-4] for a in os.listdir(os.path.join(EXT, "faa")))
    for i, acc in enumerate(accs, 1):
        fna = os.path.join(EXT, "fna", acc + ".fna")
        if not os.path.exists(fna):
            continue
        prot_kmers = set()
        for _, s in read_fasta(os.path.join(EXT, "faa", acc + ".faa")):
            for j in range(len(s) - K + 1):
                prot_kmers.add(s[j:j + K])
        contigs = [s.upper() for _, s in read_fasta(fna)]
        row = {"assembly": acc, **meta.get(acc, {})}
        for tgt, (fam, marker, pr) in targets.items():
            best = 0.0
            for s in refs[tgt]:
                tot = len(s) - K + 1
                if tot <= 0:
                    continue
                hit = sum(1 for j in range(tot) if s[j:j + K] in prot_kmers)
                best = max(best, hit / tot)
            f, r_ = pr["forward"], pr["reverse"]
            amp = 0
            for c in contigs:
                for a, b in ((f, rc(r_)), (r_, rc(f))):
                    x = c.find(a)
                    while x >= 0:
                        y = c.find(b, x)
                        if y >= 0 and AMP[0] <= y + len(b) - x <= AMP[1]:
                            amp = 1
                            break
                        x = c.find(a, x + 1)
                    if amp:
                        break
                if amp:
                    break
            row[f"{tgt}_carrier"] = int(best >= CALL)
            row[f"{tgt}_containment"] = round(best, 3)
            row[f"{tgt}_amplified"] = amp
        rows.append(row)
        if i % 50 == 0:
            print(f"  {i}/{len(accs)} assemblies", flush=True)

    with open(os.path.join(OUT, "external_pcr.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    res = {}
    drafts = [r for r in rows if r.get("group") == "draft CoNS"]
    nont = [r for r in rows if r.get("group") == "non-target"]
    for tgt in targets:
        car = [r for r in drafts if r[f"{tgt}_carrier"]]
        non = [r for r in drafts if not r[f"{tgt}_carrier"]]
        tp = sum(r[f"{tgt}_amplified"] for r in car)
        fp = sum(r[f"{tgt}_amplified"] for r in non)
        nt = sum(r[f"{tgt}_amplified"] for r in nont)
        res[tgt] = {
            "draft_carriers": len(car), "draft_non_carriers": len(non),
            "sensitivity": tp / len(car) if car else None,
            "sensitivity_ci": wilson(tp, len(car)) if car else None,
            "specificity": (len(non) - fp) / len(non) if non else None,
            "specificity_ci": wilson(len(non) - fp, len(non)) if non else None,
            "non_target_tested": len(nont), "non_target_amplified": nt}
        r_ = res[tgt]
        print(f"  {tgt:12s} drafts: {tp}/{len(car)} carriers amplified, "
              f"{fp}/{len(non)} non-carriers; non-target {nt}/{len(nont)}", flush=True)
    json.dump({"per_target": res,
               "n_draft": len(drafts), "n_non_target": len(nont),
               "species": collections.Counter(r.get("species", "") for r in drafts)},
              open(os.path.join(OUT, "external_pcr.json"), "w"), indent=1)
    print("-> out/external_pcr.csv, out/external_pcr.json", flush=True)


if __name__ == "__main__":
    main()
