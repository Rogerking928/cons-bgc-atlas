#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""What the language model is worth against simpler detectors.

Methods compared, all on the same folds and the same held-out annotation set:

  profile library   the incumbent: the 41 Pfam profiles used to build the atlas
  composition       amino-acid frequencies, length, charge, hydropathy
  3-mer counts      a bag-of-words classical text model
  nearest neighbour 10-mer containment against the training positives, which is
                    what a BLAST-style search of a peptide database gives
  ESM-2 + logistic  this paper

antiSMASH 7, DeepBGC, RODEO and NeuRiPP were not run: this machine has no
antiSMASH dependency stack, its Python is externally managed (PEP 668), and the
published NeuRiPP weights target TensorFlow 2.0-alpha. The comparison therefore
runs against the incumbent profile method and three classical baselines, and
the absence of the published deep tools is stated as a limitation rather than
papered over.
"""
import collections, csv, gzip, itertools, json, os, sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, precision_recall_curve
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
AA = "ACDEFGHIKLMNPQRSTVWY"
KD = dict(zip("AVLIPFMWGSTCYNQDEKRH",
              [1.8, 4.2, 3.8, 4.5, -1.6, 2.8, 1.9, -0.9, -0.4, -0.8,
               -0.7, 2.5, -1.3, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5, -3.2]))
KW = 10


def comp(s):
    n = len(s)
    v = [s.count(a) / n for a in AA]
    return v + [n / 120,
                (s.count("K") + s.count("R") - s.count("D") - s.count("E")) / n,
                sum(KD.get(c, 0) for c in s) / n, s.count("C") / n]


def kmer3(seqs):
    vocab = {"".join(t): i for i, t in enumerate(itertools.product(AA, repeat=2))}
    X = np.zeros((len(seqs), len(vocab)), dtype=np.float32)
    for i, s in enumerate(seqs):
        for j in range(len(s) - 1):
            k = vocab.get(s[j:j + 2])
            if k is not None:
                X[i, k] += 1
        X[i] /= max(1, len(s) - 1)
    return X


def fit_predict(Xtr, ytr, Xte):
    sc = StandardScaler().fit(Xtr)
    m = LogisticRegression(max_iter=4000, C=0.05, class_weight="balanced")
    m.fit(sc.transform(Xtr), ytr)
    return m.predict_proba(sc.transform(Xte))[:, 1], (sc, m)


def recall_at_fpr(p_pos, p_neg, fpr=0.001):
    """Recall at a fixed false-positive rate, with ties handled honestly.

    Taking the 99.9th percentile of the negatives as a cut fails when a method
    piles most negatives on one value: the cut lands on that value and the
    method appears to recover everything. Here the cut is the smallest score
    whose false-positive rate is at or below the target, and a method that
    cannot reach that rate at any cut returns nothing.
    """
    vals = np.unique(np.r_[p_pos, p_neg])[::-1]
    for v in vals:
        if (p_neg >= v).mean() <= fpr:
            return float((p_pos >= v).mean()), float(v)
    return float("nan"), float("nan")


def nn_score(seqs, pos_words):
    out = np.zeros(len(seqs))
    for i, s in enumerate(seqs):
        w = [s[j:j + KW] for j in range(len(s) - KW + 1)]
        out[i] = sum(1 for x in w if x in pos_words) / max(1, len(w))
    return out


def main():
    idx = json.load(open(os.path.join(OUT, "esm_index.json")))
    emb = np.load(os.path.join(OUT, "esm_emb.npy")).astype(np.float32)
    row = {s: i for i, s in enumerate(idx["seqs"])}
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    hold = list(csv.DictReader(open(os.path.join(OUT, "precursor_holdout.csv"))))
    y = np.array([int(r["label"]) for r in tr])
    g = np.array([int(r["cluster"]) for r in tr])
    seqs = [r["seq"] for r in tr]
    Xe = np.stack([emb[row[s]] for s in seqs])
    Xc = np.array([comp(s) for s in seqs], dtype=np.float32)
    Xk = kmer3(seqs)

    hold_seqs = [h["seq"] for h in hold if h["seq"] in row]
    He = np.stack([emb[row[s]] for s in hold_seqs])
    Hc = np.array([comp(s) for s in hold_seqs], dtype=np.float32)
    Hk = kmer3(hold_seqs)

    methods = {"ESM-2 + logistic": Xe, "3-mer counts": Xk, "composition": Xc}
    res, curves, oofs = {}, {}, {}
    for name, X in methods.items():
        oof = np.zeros(len(tr))
        for a, b in GroupKFold(n_splits=5).split(X, y, groups=g):
            oof[b] = fit_predict(X[a], y[a], X[b])[0]
        # 存下每個方法的 out-of-fold 分數：要拿最強的古典基線跟語言模型做
        # 配對比較（同一批樣本、同一組 fold），只有 AP 的點估計不夠
        oofs[name] = oof.copy()
        ap = average_precision_score(y, oof)
        pr, rc, _ = precision_recall_curve(y, oof)
        curves[name] = [pr.tolist()[::7], rc.tolist()[::7]]
        # hold-out: train on everything, score the held-out annotated peptides
        H = {"ESM-2 + logistic": He, "3-mer counts": Hk, "composition": Hc}[name]
        p_hold, (sc, m) = fit_predict(X, y, H)
        # the false-positive rate is taken from the out-of-fold negative scores,
        # not from the in-sample ones: a model scores its own training negatives
        # far lower than it scores unseen ones, which would set the cut too high
        p_neg = oof[y == 0]
        rec, thr = recall_at_fpr(p_hold, p_neg)
        res[name] = {"cv_ap": float(ap),
                     "holdout_recall": rec,
                     "operating_threshold": thr,
                     "holdout_ap": float(average_precision_score(
                         np.r_[np.ones(len(p_hold)), np.zeros(len(p_neg))],
                         np.r_[p_hold, p_neg]))}

    # nearest neighbour by word containment, cross-validated the same way
    oof = np.zeros(len(tr))
    for a, b in GroupKFold(n_splits=5).split(Xe, y, groups=g):
        pw = set()
        for i in a[y[a] == 1]:
            s = seqs[i]
            pw |= {s[j:j + KW] for j in range(len(s) - KW + 1)}
        oof[b] = nn_score([seqs[i] for i in b], pw)
    ap = average_precision_score(y, oof)
    pr, rc, _ = precision_recall_curve(y, oof)
    curves["nearest neighbour (10-mer)"] = [pr.tolist()[::7], rc.tolist()[::7]]
    pw = set()
    for i in np.where(y == 1)[0]:
        s = seqs[i]
        pw |= {s[j:j + KW] for j in range(len(s) - KW + 1)}
    p_hold = nn_score(hold_seqs, pw)
    p_neg = oof[y == 0]
    rec, thr = recall_at_fpr(p_hold, p_neg)
    res["nearest neighbour (10-mer)"] = {
        "cv_ap": float(ap), "holdout_recall": rec, "operating_threshold": thr,
        "holdout_ap": float(average_precision_score(
            np.r_[np.ones(len(p_hold)), np.zeros(len(p_neg))], np.r_[p_hold, p_neg]))}

    # the incumbent: by construction the profile library finds none of the
    # held-out peptides, because that is how they were selected
    res["Pfam profile library"] = {"cv_ap": None, "holdout_recall": 0.0,
                                   "holdout_ap": None}

    np.savez(os.path.join(OUT, "baseline_oof.npz"), y=y,
             groups=np.asarray(g),
             **{k.replace(" ", "_").replace("+", "").replace("(", "").replace(")", ""): v
                for k, v in oofs.items()})
    json.dump({"methods": res, "curves": curves,
               "n_holdout": len(hold_seqs)},
              open(os.path.join(OUT, "baselines.json"), "w"), indent=1)
    print(f"held-out annotated peptides: {len(hold_seqs)}")
    for k, v in sorted(res.items(), key=lambda kv: -(kv[1]["cv_ap"] or 0)):
        ap_s = f"{v['cv_ap']:.3f}" if v["cv_ap"] is not None else "   -  "
        print(f"  {k:28s} CV AP {ap_s}   recall at 0.1% FPR "
              f"{100*v['holdout_recall']:5.1f}%   hold-out AP "
              f"{v['holdout_ap'] if v['holdout_ap'] is None else round(v['holdout_ap'],3)}")


if __name__ == "__main__":
    main()
