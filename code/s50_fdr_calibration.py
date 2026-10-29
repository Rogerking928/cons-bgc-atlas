#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""宣稱的 q 值站不站得住，以及為什麼誘餌比真的還容易過關（2026-09-20 審閱意見）。

第一版只問「誘餌越過門檻的比例是多少」，結果是真實序列的八倍，但沒有診斷。
這一版把誘餌分層，並把兩邊的長度窗口對齊：

  decoy-from-positive   把 149 條正例逐條打散
  decoy-from-background 把面板裡的一般短 CDS 逐條打散（長度分布對齊真實集）
  real                  面板裡的真實短 CDS（同一個長度窗口）

如果「順序被打亂反而更像前驅肽」，兩組誘餌都會高；如果是抽樣或長度造成的偏差，
對齊之後就會消失。輸出 out/fdr_calibration.json。
"""
import collections, csv, json, os, random

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

import esmutil

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
N_DECOY = 1500
LEN_LO, LEN_HI = 20, 120          # 真實集與誘餌集用同一個長度窗口
SEED = 20260920
rng = random.Random(SEED)


def shuffled(s):
    lst = list(s)
    rng.shuffle(lst)
    return "".join(lst)


def main():
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    scan = list(csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv"))))

    X = esmutil.embed([r["seq"] for r in tr], verbose=False)
    Xm = np.stack([X[r["seq"]] for r in tr])
    y = np.array([int(r["label"]) for r in tr])
    sc = StandardScaler().fit(Xm)
    clf = LogisticRegression(max_iter=4000, class_weight="balanced").fit(sc.transform(Xm), y)

    pos_seqs = [r["seq"] for r in tr if int(r["label"]) == 1]
    real = {}
    for r in scan:
        real.setdefault(r["seq"], float(r["score"]))
    real = {s: v for s, v in real.items() if LEN_LO <= len(s) <= LEN_HI}
    rlen = collections.Counter(len(s) for s in real)

    # 背景誘餌：依真實集的長度分布抽樣，再逐條打散
    by_len = collections.defaultdict(list)
    for s in real:
        if s not in set(pos_seqs):
            by_len[len(s)].append(s)
    lengths = list(rlen.elements())
    rng.shuffle(lengths)
    bg = []
    for L in lengths:
        if len(bg) >= N_DECOY:
            break
        if by_len.get(L):
            bg.append(shuffled(rng.choice(by_len[L])))
    dpos = [shuffled(s) for s in pos_seqs]

    emb = esmutil.embed(bg + dpos)
    def score(seqs):
        M = np.stack([emb[s] for s in seqs])
        return clf.predict_proba(sc.transform(M))[:, 1]

    s_bg, s_pos = score(bg), score(dpos)
    s_real = np.array(list(real.values()))

    out = {"length_window": [LEN_LO, LEN_HI],
           "n_real_distinct": len(s_real), "n_decoy_background": len(s_bg),
           "n_decoy_from_positives": len(s_pos), "thresholds": {}}
    for name, thr in (("q<0.05", rep["threshold_999th_percentile_of_negatives"]),
                      ("q<0.20", rep["threshold_q20"])):
        r_rate = float((s_real >= thr).mean())
        b_rate = float((s_bg >= thr).mean())
        p_rate = float((s_pos >= thr).mean())
        out["thresholds"][name] = {
            "score": float(thr),
            "real_above": int((s_real >= thr).sum()), "real_rate": r_rate,
            "decoy_background_above": int((s_bg >= thr).sum()), "decoy_background_rate": b_rate,
            "decoy_from_positives_above": int((s_pos >= thr).sum()),
            "decoy_from_positives_rate": p_rate,
            "empirical_fdp_background": (b_rate / r_rate) if r_rate else float("nan")}
        print(f"  {name}: real {100 * r_rate:.2f}%  background decoys {100 * b_rate:.2f}%  "
              f"decoys made from positives {100 * p_rate:.1f}%  "
              f"-> FDP estimate {b_rate / r_rate if r_rate else float('nan'):.2f}",
              flush=True)

    grid = np.linspace(0.3, 1.0, 71)
    out["curve"] = [{"threshold": float(t),
                     "real_rate": float((s_real >= t).mean()),
                     "decoy_background_rate": float((s_bg >= t).mean()),
                     "decoy_from_positives_rate": float((s_pos >= t).mean())}
                    for t in grid]
    json.dump(out, open(os.path.join(OUT, "fdr_calibration.json"), "w"), indent=1)
    print("-> out/fdr_calibration.json", flush=True)


if __name__ == "__main__":
    main()
