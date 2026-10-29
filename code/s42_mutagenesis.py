#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""因果測試：把雙甘胺酸位點改掉，看偵測器的分數掉多少（2026-09-20 審閱意見）。

Integrated gradients 只說「模型看了哪裡」，是相關性。這支做的是擾動：

  1. 候選與 profile 抓到的前驅肽：把最後一個 G[GAS] 改成 A[AAS]（兩個位置都動），
     另外做一組「等長度、等數量、隨機位置」的對照突變當基準。
  2. 陰性：在同樣的相對位置插入一個 GG，看分數升多少。

分數用的是正文那個凍結模型（ESM-2 150M ＋ 邏輯迴歸，全訓練集重訓一次）。
輸出 out/mutagenesis.json、out/mutagenesis.csv。
"""
import csv, json, os, random, re

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

import esmutil

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SEED = 20260920
GG = re.compile(r"G[GAS]")
rng = random.Random(SEED)


def last_gg(s):
    m = None
    for mm in GG.finditer(s):
        m = mm
    return m


def mutate_gg(s):
    m = last_gg(s)
    if m is None:
        return None, None
    i = m.start()
    return s[:i] + "AA" + s[i + 2:], i


def mutate_random(s, i_skip):
    """對照：改同樣兩個殘基，但位置隨機（避開 GG 本身）。"""
    pos = [i for i in range(len(s) - 1) if abs(i - i_skip) > 2]
    if not pos:
        return None
    i = rng.choice(pos)
    return s[:i] + "AA" + s[i + 2:]


def insert_gg(s):
    """陰性：在序列 40% 處插入 GG（前導肽末端的典型位置）。"""
    i = max(2, int(len(s) * 0.4))
    return s[:i] + "GG" + s[i:]


def main():
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    thr20 = rep["threshold_q20"]
    scan = list(csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv"))))

    # 凍結模型，全訓練集重訓（和正文的偵測器一致）
    X = esmutil.embed([r["seq"] for r in tr], verbose=False)
    Xm = np.stack([X[r["seq"]] for r in tr])
    y = np.array([int(r["label"]) for r in tr])
    sc = StandardScaler().fit(Xm)
    clf = LogisticRegression(max_iter=4000, class_weight="balanced").fit(sc.transform(Xm), y)

    def score(seqs):
        """回 log-odds，不回機率：候選的機率已經貼在 0 或 1，機率差會被天花板吃掉，
        看不出擾動有沒有效果。"""
        e = esmutil.embed(seqs)
        M = np.stack([e[s] for s in seqs])
        return clf.decision_function(sc.transform(M))

    seen, groups = set(), {"candidate": [], "profile precursor": []}
    for r in scan:
        s = r["seq"]
        if s in seen:
            continue
        seen.add(s)
        if r["known_domain"]:
            groups["profile precursor"].append(s)
        elif float(r["score"]) >= thr20:
            groups["candidate"].append(s)
    negs = [r["seq"] for r in tr if int(r["label"]) == 0]
    rng.shuffle(negs)
    groups["negative"] = negs[:120]

    rows, res = [], {}
    for grp, seqs in groups.items():
        seqs = [s for s in seqs if last_gg(s) or grp == "negative"][:120]
        if not seqs:
            continue
        base = score(seqs)
        if grp == "negative":
            mut = [insert_gg(s) for s in seqs]
            delta = score(mut) - base
            kind = "GG inserted"
            # 對照：插入同樣兩個殘基，但插 AA（不是 motif），看單純「變長兩個殘基」
            # 本身會不會抬分數
            ctrl = [s[:max(2, int(len(s) * 0.4))] + "AA" + s[max(2, int(len(s) * 0.4)):]
                    for s in seqs]
            ctrl_delta = score(ctrl) - base
        else:
            pairs = [mutate_gg(s) for s in seqs]
            mut = [p[0] for p in pairs]
            idx = [p[1] for p in pairs]
            delta = score(mut) - base
            ctrl = [mutate_random(s, i) for s, i in zip(seqs, idx)]
            ok = [j for j, c in enumerate(ctrl) if c]
            ctrl_delta = score([ctrl[j] for j in ok]) - base[ok]
            kind = "GG -> AA"
        res[grp] = {"n": len(seqs), "kind": kind,
                    "baseline_median_logodds": float(np.median(base)),
                    "delta_median": float(np.median(delta)),
                    "delta_mean": float(np.mean(delta)),
                    "fraction_dropping": float(np.mean(delta < 0)),
                    "control_delta_median": (None if ctrl_delta is None
                                             else float(np.median(ctrl_delta)))}
        for s, b, d in zip(seqs, base, delta):
            rows.append(dict(group=grp, kind=kind, seq=s, base=round(float(b), 4),
                             mutated=round(float(b + d), 4), delta=round(float(d), 4)))
        print(f"  {grp:18s} n={len(seqs):4d} {kind:12s} "
              f"median log-odds {np.median(base):+.2f} -> {np.median(base + delta):+.2f} "
              f"(median change {np.median(delta):+.3f}; "
              f"control {('n/a' if ctrl_delta is None else f'{np.median(ctrl_delta):+.3f}')})",
              flush=True)

    json.dump(res, open(os.path.join(OUT, "mutagenesis.json"), "w"), indent=1)
    with open(os.path.join(OUT, "mutagenesis.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("-> out/mutagenesis.json, out/mutagenesis.csv", flush=True)


if __name__ == "__main__":
    main()
