#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把「拿掉一個定序型」當敏感度分析（2026-09-20 審閱建議）。

F001 與 erm(A) 的負相關看起來很強，但 ST23 這一支 22 株全部帶 erm(A)＋mecA 而且
完全不帶 F001。把 ST23 拿掉之後關聯還在不在，是這件事最直接的檢定。
輸出 out/st_sensitivity.json。
"""
import collections, csv, json, os

from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
CALL = 0.5
PAIRS = [("epidermidis", "F001", "erm(A)"), ("epidermidis", "F001", "ant(9)-Ia"),
         ("epidermidis", "F001", "blaZ"), ("epidermidis", "F001", "qacC")]


def main():
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    amr = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv")))}
    mlst = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "mlst.csv")))}

    def has(a, d):
        return float(amr.get(a, {}).get(d, 0) or 0) >= CALL

    def table(sub, fam, det):
        a = sum(1 for x in sub if mat[x].get(fam) and has(x, det))
        b = sum(1 for x in sub if mat[x].get(fam) and not has(x, det))
        c = sum(1 for x in sub if not mat[x].get(fam) and has(x, det))
        d = sum(1 for x in sub if not mat[x].get(fam) and not has(x, det))
        orr, p = stats.fisher_exact([[a, b], [c, d]])
        return dict(n=len(sub), a=a, b=b, c=c, d=d, odds=float(orr), p=float(p))

    out = {}
    for sp, fam, det in PAIRS:
        pool = [a for a, r in mat.items() if r["species"] == sp]
        sts = collections.Counter(mlst.get(a, {}).get("ST", "") for a in pool)
        res = {"all": table(pool, fam, det), "leave_one_ST_out": {}}
        for st, n in sts.most_common():
            if not st or n < 5:
                continue
            sub = [a for a in pool if mlst.get(a, {}).get("ST") != st]
            res["leave_one_ST_out"][f"ST{st}"] = table(sub, fam, det) | {"removed": n}
        worst = max(res["leave_one_ST_out"].items(), key=lambda kv: kv[1]["p"])
        res["worst_ST"] = worst[0]
        res["worst_p"] = worst[1]["p"]
        res["worst_odds"] = worst[1]["odds"]
        out[f"{sp}:{fam}:{det}"] = res
        print(f"  {sp} {fam} x {det}: crude OR {res['all']['odds']:.3f} "
              f"(P {res['all']['p']:.2g}); worst when {worst[0]} removed "
              f"(n={worst[1]['removed']}): OR {worst[1]['odds']:.3f}, P {worst[1]['p']:.3g}",
              flush=True)
    json.dump(out, open(os.path.join(OUT, "st_sensitivity.json"), "w"), indent=1)
    print("-> out/st_sensitivity.json", flush=True)


if __name__ == "__main__":
    main()
