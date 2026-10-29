#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""胜肽層的領域基準：規則式前驅肽啟發法（2026-09-20 審閱意見）。

已發表的胜肽層 RiPP 預測器沒有一個能在這台機器離線跑（NeuRiPP 的權重是
TensorFlow 2.0-alpha、RiPPMiner 與 DeepRiPP 只有網頁服務），所以這裡實作一組
**規則式的啟發法**，形式與 BAGEL／antiSMASH 在前驅肽偵測上使用的規則同類：

  1. 長度 20–120 aa
  2. 前導肽末端有雙甘胺酸（G[GAS]）或 P[QRK]／[AVE]R 切割基序
  3. 切下來的核心 12–45 aa
  4. 核心富含 Cys/Ser/Thr（羊毛硫肽與硫肽的修飾位點）
  5. 前導肽淨電荷為負、核心淨電荷為正或中性

分數＝滿足的條件數（0–5），可以畫 PR 曲線。**這不是那些工具本身**，只是同類規則，
用來回答「語言模型贏的是不是只有樸素的統計基準」。
輸出 out/rule_baseline.json。
"""
import csv, json, os, re

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
CLEAVE = re.compile(r"G[GAS]|P[QRK]|[AVE]R")
POS = set("KR")
NEG = set("DE")


def core_of(s):
    best = ""
    for m in CLEAVE.finditer(s):
        c = s[m.end():]
        if 12 <= len(c) <= 45:
            best = c
    return best


def score(s):
    n = 0
    if 20 <= len(s) <= 120:
        n += 1
    core = core_of(s)
    if core:
        n += 1
        if 12 <= len(core) <= 45:
            n += 1
        frac = sum(core.count(c) for c in "CST") / len(core)
        if frac >= 0.15:
            n += 1
        leader = s[:len(s) - len(core)]
        lc = sum(1 for c in leader if c in POS) - sum(1 for c in leader if c in NEG)
        cc = sum(1 for c in core if c in POS) - sum(1 for c in core if c in NEG)
        if lc < 0 <= cc:
            n += 1
    return n


def main():
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    y = np.array([int(r["label"]) for r in tr])
    s = np.array([score(r["seq"]) for r in tr], dtype=float)
    # 同分的先後不該影響 AP：加一點點 jitter 之前先確認每一級的正例比例
    res = {"n": len(y), "n_positives": int(y.sum()),
           "ap": float(average_precision_score(y, s)),
           "auroc": float(roc_auc_score(y, s)),
           "by_score": {}}
    for k in range(6):
        m = s == k
        if m.sum():
            res["by_score"][k] = {"n": int(m.sum()), "positives": int(y[m].sum()),
                                  "precision": float(y[m].mean())}
    hold = list(csv.DictReader(open(os.path.join(OUT, "precursor_holdout.csv"))))
    hs = np.array([score(h["seq"]) for h in hold], dtype=float)
    neg = s[y == 0]
    thr = 4
    res["holdout"] = {"n": len(hs), "recall_at_score>=4": float((hs >= thr).mean()),
                      "negatives_at_score>=4": float((neg >= thr).mean())}
    json.dump(res, open(os.path.join(OUT, "rule_baseline.json"), "w"), indent=1)
    print(f"  rule-based heuristic: AP {res['ap']:.3f}  AUROC {res['auroc']:.3f}", flush=True)
    for k, v in res["by_score"].items():
        print(f"    score {k}: n={v['n']:6d} positives={v['positives']:4d} "
              f"precision={v['precision']:.3f}", flush=True)
    print(f"  held-out annotated set: {100 * res['holdout']['recall_at_score>=4']:.0f}% "
          f"reach score ≥4, against "
          f"{100 * res['holdout']['negatives_at_score>=4']:.1f}% of negatives", flush=True)
    print("-> out/rule_baseline.json", flush=True)


if __name__ == "__main__":
    main()
