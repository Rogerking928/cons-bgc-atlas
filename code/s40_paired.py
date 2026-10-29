#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""兩個配對比較（審閱意見 2026-09-20）。

1. **語言模型 vs 最強的古典基線**：之前只報了對 composition（最弱那個）的差值，
   看起來像挑了有利的比較。這裡改成對 3-mer counts（最強）做配對 bootstrap，
   同一批樣本、同一組 fold 的 out-of-fold 分數。

2. **凍結 vs LoRA 微調**：兩根長條不夠，要 across-fold 的配對差值與區間。
   重現 s33 的抽樣（同一個 SEED、同樣 3,000 個陰性）與同一組 GroupKFold，
   用同一份 35M 嵌入跑凍結版，逐折相減。

輸出 out/paired_tests.json。
"""
import csv, json, os, random

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SEED, N_NEG = 0, 3000           # 與 s33_finetune.py 相同
B = 2000                        # bootstrap 次數
rng = np.random.default_rng(20260920)


def paired_bootstrap(y, s1, s2, groups):
    """以序列簇為單位重抽（同一簇的樣本一起進出），回傳差值與百分位區間。"""
    obs = average_precision_score(y, s1) - average_precision_score(y, s2)
    uniq = np.unique(groups)
    idx_by_g = {gg: np.where(groups == gg)[0] for gg in uniq}
    diffs = []
    for _ in range(B):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([idx_by_g[gg] for gg in pick])
        if y[idx].sum() < 2 or y[idx].sum() == len(idx):
            continue
        diffs.append(average_precision_score(y[idx], s1[idx])
                     - average_precision_score(y[idx], s2[idx]))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    p = 2 * min(np.mean(np.asarray(diffs) <= 0), np.mean(np.asarray(diffs) >= 0))
    return dict(difference=float(obs), ci=[float(lo), float(hi)],
                p_bootstrap=float(min(1.0, p)), n_resamples=len(diffs))


def frozen_folds():
    """s33 的抽樣與 fold，改用凍結的 35M 嵌入 ＋ 邏輯迴歸。"""
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    X = np.load(os.path.join(OUT, "emb_train_35M.npy")).astype(np.float32)
    pos = [i for i, r in enumerate(tr) if int(r["label"]) == 1]
    neg = [i for i, r in enumerate(tr) if int(r["label"]) == 0]
    random.Random(SEED).shuffle(neg)
    sub = pos + neg[:N_NEG]
    y = np.array([int(tr[i]["label"]) for i in sub])
    g = np.array([int(tr[i]["cluster"]) for i in sub])
    Xs = X[sub]
    aps = []
    for a, b in GroupKFold(n_splits=5).split(Xs, y, groups=g):
        sc = StandardScaler().fit(Xs[a])
        m = LogisticRegression(max_iter=4000, C=1.0, class_weight="balanced")
        m.fit(sc.transform(Xs[a]), y[a])
        p = m.predict_proba(sc.transform(Xs[b]))[:, 1]
        aps.append(float(average_precision_score(y[b], p)))
    return aps, int(y.sum()), len(y)


def main():
    out = {}

    # 1. 語言模型 vs 最強古典基線
    z = np.load(os.path.join(OUT, "baseline_oof.npz"))
    y = z["y"]
    groups = z["groups"]
    esm = z["ESM-2__logistic"] if "ESM-2__logistic" in z else z[[k for k in z if "ESM" in k][0]]
    k3 = z[[k for k in z if "3-mer" in k][0]]
    comp = z[[k for k in z if "composition" in k][0]]
    out["language_model_vs_3mer"] = paired_bootstrap(y, esm, k3, groups)
    out["language_model_vs_composition"] = paired_bootstrap(y, esm, comp, groups)
    out["ap"] = {"language_model": float(average_precision_score(y, esm)),
                 "3mer": float(average_precision_score(y, k3)),
                 "composition": float(average_precision_score(y, comp))}

    # 2. 凍結 vs 微調（逐折配對）
    ft = json.load(open(os.path.join(OUT, "finetune_ap.json")))
    fr, npos, ntot = frozen_folds()
    d = np.array(ft["folds"]) - np.array(fr)
    n = len(d)
    sd = float(d.std(ddof=1))
    se = sd / np.sqrt(n)
    t975 = 2.776                       # t_{0.975, df=4}
    out["finetuned_vs_frozen_same_folds"] = {
        "n_sequences": ntot, "n_positives": npos,
        "frozen_folds": fr, "finetuned_folds": ft["folds"],
        "frozen_mean": float(np.mean(fr)), "finetuned_mean": float(np.mean(ft["folds"])),
        "paired_difference": float(d.mean()),
        "ci_t": [float(d.mean() - t975 * se), float(d.mean() + t975 * se)],
        "note": "五折的配對差值，t 區間（df=4）。n 很小，區間跨 0 也不奇怪。"}

    json.dump(out, open(os.path.join(OUT, "paired_tests.json"), "w"), indent=1)
    for k, v in out.items():
        print(k, json.dumps(v)[:200], flush=True)


if __name__ == "__main__":
    main()
