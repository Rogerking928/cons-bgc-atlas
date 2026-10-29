#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""任務到底有多難：難負例、相似度分層、學習曲線（2026-09-20 審閱意見）。

原本的負例是「離機器 20 kb 以外、和正例不像的短蛋白」——那是 easy negatives，
模型可能只是在學「這是不是一條小的分泌型胜肽」。這支腳本把三件事補起來：

1. **難負例**：(a) 落在基因簇 10 kb 內但沒有前驅肽域的短 CDS、(b) 保留胺基酸組成的
   洗牌序列、(c) 長度分布對齊正例的隨機短 CDS。easy 與 hard 分開報 AP。
2. **長度混淆**：只用「長度」當唯一特徵的邏輯迴歸當基準，看 AP 有多少是長度給的。
3. **相似度分層**：依「與最近的訓練正例的 10-mer 相似度」分 bin 報 AP。
4. **學習曲線**：AP 對正例數量，看 149 條到底夠不夠。

輸出 out/hard_negatives.json。
"""
import collections, csv, json, os, random, re

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

import esmutil

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
SEED = 20260920
NEAR = 10000
K = 10
NAME_RE = re.compile(r"bacteriocin|lantibiotic|lanthipeptide|microcin|"
                     r"sactipeptide|thiopeptide|linear azol", re.I)
rng = random.Random(SEED)


def fit_oof(X, y, g):
    oof = np.zeros(len(y))
    for a, b in GroupKFold(n_splits=5).split(X, y, groups=g):
        sc = StandardScaler().fit(X[a])
        m = LogisticRegression(max_iter=4000, class_weight="balanced")
        m.fit(sc.transform(X[a]), y[a])
        oof[b] = m.predict_proba(sc.transform(X[b]))[:, 1]
    return oof


def kmers(s, k=K):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def main():
    idx = json.load(open(os.path.join(OUT, "esm_index.json")))
    emb = np.load(os.path.join(OUT, "esm_emb.npy"), mmap_mode="r")
    row = {s: i for i, s in enumerate(idx["seqs"])}
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    pos = [r for r in tr if int(r["label"]) == 1]
    neg = [r for r in tr if int(r["label"]) == 0]
    pos_seqs = {r["seq"] for r in pos}
    plen = [len(r["seq"]) for r in pos]

    # ---------------------------------------------------------------- 難負例
    scan = list(csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv"))))
    seen, near_neg = set(), []
    for r in scan:
        s = r["seq"]
        if s in seen or s in pos_seqs or r["known_domain"]:
            continue
        if NAME_RE.search(r["product"] or ""):
            continue                      # 這些是保留集，不能當負例
        d = r["nearest_anchor_bp"]
        if d in ("", "None") or int(float(d)) > NEAR:
            continue
        seen.add(s)
        near_neg.append(s)
    # 組成洗牌：保留胺基酸組成，破壞順序
    shuffled = []
    for r in pos:
        lst = list(r["seq"])
        rng.shuffle(lst)
        shuffled.append("".join(lst))
    # 長度對齊的隨機短 CDS
    by_len = collections.defaultdict(list)
    for r in neg:
        by_len[len(r["seq"])].append(r["seq"])
    lmatched = []
    for L in plen:
        cand = [s for dl in range(0, 6)
                for s in by_len.get(L + dl, []) + by_len.get(L - dl, [])]
        if cand:
            lmatched.append(rng.choice(cand))

    # 洗牌序列不在面板裡，要另外算嵌入
    extra = esmutil.embed(shuffled)
    vec = {}

    def get(s):
        if s in row:
            return np.asarray(emb[row[s]], dtype=np.float32)
        return extra.get(s)

    def block(seqs):
        return [s for s in seqs if get(s) is not None]

    sets = {"easy (far from machinery)": [r["seq"] for r in neg],
            "hard: inside a cluster window": block(near_neg),
            "hard: composition-shuffled positives": block(shuffled),
            "hard: length-matched short CDS": block(lmatched)}
    res = {"n_positives": len(pos), "negative_sets": {}}

    P = np.stack([np.asarray(emb[row[r["seq"]]], dtype=np.float32) for r in pos])
    gpos = np.array([int(r["cluster"]) for r in pos])
    for name, ss in sets.items():
        ss = [s for s in dict.fromkeys(ss) if s not in pos_seqs]
        if len(ss) < 50:
            continue
        N = np.stack([get(s) for s in ss])
        X = np.vstack([P, N])
        y = np.r_[np.ones(len(P)), np.zeros(len(N))]
        # 負例自己一組一群，避免同群跨折
        g = np.r_[gpos, np.arange(len(N)) + gpos.max() + 1]
        oof_set = fit_oof(X, y, g)
        ap = float(average_precision_score(y, oof_set))
        auroc = float(roc_auc_score(y, oof_set))
        # 只用長度的基準
        Xl = np.r_[[len(r["seq"]) for r in pos], [len(s) for s in ss]].reshape(-1, 1).astype(float)
        ap_len = float(average_precision_score(y, fit_oof(Xl, y, g)))
        res["negative_sets"][name] = {"n_negatives": len(ss), "ap": ap,
                                      "auroc": auroc, "ap_length_only": ap_len,
                                      "prevalence": float(y.mean())}
        # AP 會隨正例比例上下跑，所以每一組都要把 prevalence（＝隨機猜的 AP）
        # 一起報，另外報與 prevalence 無關的 AUROC
        print(f"  {name:38s} n={len(ss):6d} prev={y.mean():.3f}  AP {ap:.3f} "
              f"(chance {y.mean():.3f})  AUROC {auroc:.3f}  length-only AP {ap_len:.3f}",
              flush=True)

    # ------------------------------------------------- 相似度分層（easy 負例）
    oof = np.load(os.path.join(OUT, "precursor_oof_language.npy"))
    y_all = np.load(os.path.join(OUT, "precursor_y.npy"))
    pos_km = [kmers(r["seq"]) for r in pos]
    bins = {"<30%": [], "30-50%": [], "50-70%": [], ">=70%": []}
    for i, r in enumerate(tr):
        if int(r["label"]) != 1:
            continue
        ks = kmers(r["seq"])
        best = 0.0
        for j, other in enumerate(pos_km):
            if j == i or not ks:
                continue
            best = max(best, len(ks & other) / len(ks))
        b = ("<30%" if best < 0.3 else "30-50%" if best < 0.5
             else "50-70%" if best < 0.7 else ">=70%")
        bins[b].append(i)
    strat = {}
    negs_idx = [i for i, r in enumerate(tr) if int(r["label"]) == 0]
    for b, ii in bins.items():
        if not ii:
            continue
        keep = ii + negs_idx
        # 每個 bin 的 AP 要有區間，n=9 那一格沒有誤差棒會被問
        yk, sk = y_all[keep], oof[keep]
        boots = []
        rs = np.random.default_rng(20260920)
        for _ in range(500):
            idx = rs.integers(0, len(yk), len(yk))
            if yk[idx].sum() < 2:
                continue
            boots.append(average_precision_score(yk[idx], sk[idx]))
        lo, hi = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))) \
            if boots else (None, None)
        strat[b] = {"n_positives": len(ii),
                    "ap": float(average_precision_score(yk, sk)),
                    "ap_ci": [lo, hi]}
        print(f"  identity {b:8s} n_pos={len(ii):4d}  AP {strat[b]['ap']:.3f}", flush=True)
    res["identity_strata"] = strat

    # ------------------------------------------------------------ 學習曲線
    curve = []
    Xall = np.stack([np.asarray(emb[row[r["seq"]]], dtype=np.float32)
                     for r in tr if r["seq"] in row])
    keep = [i for i, r in enumerate(tr) if r["seq"] in row]
    yall = np.array([int(tr[i]["label"]) for i in keep])
    gall = np.array([int(tr[i]["cluster"]) for i in keep])
    pos_groups = sorted({gall[i] for i in range(len(yall)) if yall[i] == 1})
    for frac in (0.1, 0.25, 0.5, 0.75, 1.0):
        take = set(rng.sample(pos_groups, max(2, int(len(pos_groups) * frac))))
        sel = [i for i in range(len(yall)) if yall[i] == 0 or gall[i] in take]
        ap = float(average_precision_score(
            yall[sel], fit_oof(Xall[sel], yall[sel], gall[sel])))
        n = int(yall[sel].sum())
        curve.append({"fraction": frac, "n_positives": n, "ap": ap})
        print(f"  learning curve {n:4d} positives  AP {ap:.3f}", flush=True)
    res["learning_curve"] = curve

    json.dump(res, open(os.path.join(OUT, "hard_negatives.json"), "w"), indent=1)
    print("-> out/hard_negatives.json", flush=True)


if __name__ == "__main__":
    main()
