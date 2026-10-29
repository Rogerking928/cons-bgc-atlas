#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Confidence intervals for the model ladder.

The ladder claims that size stops paying by 150M and that fine-tuning pays
more. Neither claim survives without an interval, so the out-of-fold
predictions are recomputed from the stored embeddings and bootstrapped over
sequence clusters, the same resampling used for the main comparison. The
fine-tuned model is the exception: its folds would have to be retrained to
bootstrap, so its five per-fold values are carried instead.
"""
import csv, json, os, re, sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
B = 1000


def oof_for(X, y, g):
    out = np.zeros(len(y))
    for a, b in GroupKFold(n_splits=5).split(X, y, groups=g):
        sc = StandardScaler().fit(X[a])
        m = LogisticRegression(max_iter=4000, C=0.05, class_weight="balanced")
        m.fit(sc.transform(X[a]), y[a])
        out[b] = m.predict_proba(sc.transform(X[b]))[:, 1]
    return out


def main():
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    y = np.array([int(r["label"]) for r in tr])
    g = np.array([int(r["cluster"]) for r in tr])
    seqs = [r["seq"] for r in tr]
    idx = json.load(open(os.path.join(OUT, "esm_index.json")))
    row = {s: i for i, s in enumerate(idx["seqs"])}
    emb150 = np.load(os.path.join(OUT, "esm_emb.npy")).astype(np.float32)

    AA = "ACDEFGHIKLMNPQRSTVWY"
    KD = dict(zip("AVLIPFMWGSTCYNQDEKRH",
                  [1.8, 4.2, 3.8, 4.5, -1.6, 2.8, 1.9, -0.9, -0.4, -0.8,
                   -0.7, 2.5, -1.3, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5, -3.2]))

    def comp(s):
        n = len(s)
        return [s.count(a) / n for a in AA] + [
            n / 120, (s.count("K") + s.count("R") - s.count("D") - s.count("E")) / n,
            sum(KD.get(c, 0) for c in s) / n, s.count("C") / n]

    mats = {"composition": np.array([comp(s) for s in seqs], dtype=np.float32),
            "150M": np.stack([emb150[row[s]] for s in seqs])}
    for tag in ("8M", "35M"):
        p = os.path.join(OUT, f"emb_train_{tag}.npy")
        if os.path.exists(p):
            mats[tag] = np.load(p).astype(np.float32)

    rng = np.random.default_rng(20260919)
    uniq = np.unique(g)
    idx_of = {u: np.where(g == u)[0] for u in uniq}
    res = {}
    for tag, X in mats.items():
        oof = oof_for(X, y, g)
        aps = []
        for _ in range(B):
            pick = rng.choice(uniq, size=len(uniq), replace=True)
            ii = np.concatenate([idx_of[u] for u in pick])
            if y[ii].sum() < 5:
                continue
            aps.append(average_precision_score(y[ii], oof[ii]))
        res[tag] = {"ap": float(average_precision_score(y, oof)),
                    "ci": [float(np.percentile(aps, 2.5)),
                           float(np.percentile(aps, 97.5))]}
        print(f"  {tag:12s} AP {res[tag]['ap']:.3f} "
              f"({res[tag]['ci'][0]:.3f}–{res[tag]['ci'][1]:.3f})", flush=True)

    folds = [float(m.group(1)) for m in
             re.finditer(r"fold \d: AP (\d\.\d+)",
                         open(os.path.join(OUT, "finetune.log")).read())]
    ft = json.load(open(os.path.join(OUT, "finetune_summary.json")))
    res["LoRA 35M"] = {"ap": ft["finetuned_cv_ap"], "ci": None, "folds": folds}
    print(f"  LoRA 35M     AP {ft['finetuned_cv_ap']:.3f}; folds "
          + ", ".join(f"{v:.3f}" for v in folds), flush=True)
    json.dump(res, open(os.path.join(OUT, "ladder_ci.json"), "w"), indent=1)
    print("-> out/ladder_ci.json", flush=True)


if __name__ == "__main__":
    main()
