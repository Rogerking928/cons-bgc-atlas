#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Train the precursor detector, then scan every short protein in the panel.

Evaluation is deliberately hard on the model:
  * five-fold cross-validation grouped by sequence cluster, so a near-duplicate
    of a training peptide is never scored as a held-out success
  * a cross-source test: train on the curated UniProt peptides alone and test on
    the panel's own precursors, and the reverse
  * a baseline of composition features, to show what the language model adds
  * an external check that uses no label at all: do the high-scoring proteins
    that no Pfam profile recognises sit next to maturation and export genes more
    often than short proteins in general?
"""
import collections, csv, gzip, json, os, sys

import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FAA = os.path.join(HERE, "data", "faa")
AA = "ACDEFGHIKLMNPQRSTVWY"
ANCHOR_DOMS = {"C39", "LanM_cyc", "LanC", "YcaO", "LanB_N", "LanB_C", "Peptidase_S8"}
# 要和 s17 的那份一致，否則「Pfam 本來就抓到的前驅肽」會被當成模型新找到的
# （2026-09-20 踩過：漏掉 Gallidermin，F004 的 39 個座全被算成模型的功勞）
PRECURSOR_DOMS = {"Mersacidin", "Bacteriocin_II", "Bacteriocin_IIc", "Bacteriocin_IId",
                  "Bacteriocin_IIi", "BacteriocIIc_cy", "Antimicrobial17",
                  "L_biotic_typeA", "MccV",
                  "Gallidermin", "Antimicrobial18", "Lanthi_prec",
                  "Lactococcin_972"}
NEAR = 10000
KD = dict(zip("AVLIPFMWGSTCYNQDEKRH",
              [1.8, 4.2, 3.8, 4.5, -1.6, 2.8, 1.9, -0.9, -0.4, -0.8,
               -0.7, 2.5, -1.3, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5, -3.2]))


def comp_features(s):
    n = len(s)
    v = [s.count(a) / n for a in AA]
    charge = (s.count("K") + s.count("R") - s.count("D") - s.count("E")) / n
    gravy = sum(KD.get(c, 0) for c in s) / n
    return v + [n / 120, charge, gravy, s.count("C") / n,
                float("GG" in s[:40]), float("GA" in s[:40]), float(s.startswith("M"))]


def read_fasta(path):
    n, b = None, []
    for line in open(path):
        if line.startswith(">"):
            if n:
                yield n, "".join(b)
            n, b = line[1:].split()[0], []
        else:
            b.append(line.strip())
    if n:
        yield n, "".join(b)


def fit(X, y):
    sc = StandardScaler().fit(X)
    m = LogisticRegression(max_iter=4000, C=0.05, class_weight="balanced")
    m.fit(sc.transform(X), y)
    return sc, m


def score(sc, m, X):
    return m.predict_proba(sc.transform(X))[:, 1]


def main():
    idx = json.load(open(os.path.join(OUT, "esm_index.json")))
    emb = np.load(os.path.join(OUT, "esm_emb.npy")).astype(np.float32)
    row_of = {s: i for i, s in enumerate(idx["seqs"])}
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    y = np.array([int(r["label"]) for r in tr])
    g = np.array([int(r["cluster"]) for r in tr])
    src = np.array([r["source"] for r in tr])
    Xe = np.stack([emb[row_of[r["seq"]]] for r in tr])
    Xc = np.array([comp_features(r["seq"]) for r in tr], dtype=np.float32)
    print(f"{len(tr)} sequences, {y.sum()} positives, {len(set(g[y==1]))} positive clusters",
          flush=True)

    report = {"model": idx["model"], "n_train": len(tr), "n_pos": int(y.sum()),
              "n_pos_clusters": int(len(set(g[y == 1])))}
    for name, X in (("language model", Xe), ("composition baseline", Xc)):
        aps, aucs = [], []
        oof = np.zeros(len(tr))
        for tr_i, te_i in GroupKFold(n_splits=5).split(X, y, groups=g):
            sc, m = fit(X[tr_i], y[tr_i])
            p = score(sc, m, X[te_i])
            oof[te_i] = p
            aps.append(average_precision_score(y[te_i], p))
            aucs.append(roc_auc_score(y[te_i], p))
        report[name] = {"auprc": float(np.mean(aps)), "auprc_sd": float(np.std(aps)),
                        "auroc": float(np.mean(aucs)), "auroc_sd": float(np.std(aucs))}
        print(f"  {name:22s} AUPRC {np.mean(aps):.3f} +- {np.std(aps):.3f}   "
              f"AUROC {np.mean(aucs):.3f}", flush=True)
        np.save(os.path.join(OUT, "precursor_oof_%s.npy"
                             % name.split()[0]), oof)
        np.save(os.path.join(OUT, "precursor_y.npy"), y)

    # bootstrap confidence intervals on the out-of-fold predictions, resampled
    # by sequence cluster so a large cluster cannot be counted many times over
    oof_lm = np.load(os.path.join(OUT, "precursor_oof_language.npy"))
    oof_bl = np.load(os.path.join(OUT, "precursor_oof_composition.npy"))
    rng = np.random.default_rng(20260919)
    groups = np.array(g)
    uniq_g = np.unique(groups)
    idx_of = {u: np.where(groups == u)[0] for u in uniq_g}
    aps_lm, aps_bl, diffs = [], [], []
    for _ in range(1000):
        pick = rng.choice(uniq_g, size=len(uniq_g), replace=True)
        ii = np.concatenate([idx_of[u] for u in pick])
        if y[ii].sum() < 5 or y[ii].sum() == len(ii):
            continue
        a1 = average_precision_score(y[ii], oof_lm[ii])
        a2 = average_precision_score(y[ii], oof_bl[ii])
        aps_lm.append(a1); aps_bl.append(a2); diffs.append(a1 - a2)
    def ci(v):
        return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    report["bootstrap"] = {
        "ap_language_model": float(average_precision_score(y, oof_lm)),
        "ap_language_model_ci": ci(aps_lm),
        "ap_composition": float(average_precision_score(y, oof_bl)),
        "ap_composition_ci": ci(aps_bl),
        "ap_difference": float(np.mean(diffs)), "ap_difference_ci": ci(diffs),
        "p_two_sided": float(2 * min((np.array(diffs) <= 0).mean(),
                                     (np.array(diffs) >= 0).mean()))}
    print(f"  AP language model {report['bootstrap']['ap_language_model']:.3f} "
          f"95% CI {ci(aps_lm)[0]:.3f}-{ci(aps_lm)[1]:.3f}; composition "
          f"{report['bootstrap']['ap_composition']:.3f} "
          f"({ci(aps_bl)[0]:.3f}-{ci(aps_bl)[1]:.3f}); difference "
          f"{np.mean(diffs):.3f} ({ci(diffs)[0]:.3f}-{ci(diffs)[1]:.3f}), "
          f"P = {report['bootstrap']['p_two_sided']:.4f}", flush=True)

    # calibration of the out-of-fold scores
    iso = IsotonicRegression(out_of_bounds="clip").fit(oof_lm, y)
    cal = iso.predict(oof_lm)
    bins = np.linspace(0, 1, 11)
    which = np.clip(np.digitize(oof_lm, bins) - 1, 0, 9)
    rel = [[float(np.mean(oof_lm[which == b])) if (which == b).any() else None,
            float(np.mean(y[which == b])) if (which == b).any() else None,
            int((which == b).sum())] for b in range(10)]
    report["calibration"] = {
        "brier_raw": float(brier_score_loss(y, oof_lm)),
        "brier_isotonic": float(brier_score_loss(y, cal)),
        "reliability": rel}
    print(f"  Brier {report['calibration']['brier_raw']:.4f} raw, "
          f"{report['calibration']['brier_isotonic']:.4f} after isotonic", flush=True)

    # cross-source generalisation. The negatives are split in half as well, so
    # the test half is unseen on both sides of the label.
    rng0 = np.random.default_rng(20260918)
    neg = np.where(y == 0)[0]
    rng0.shuffle(neg)
    neg_tr, neg_te = neg[:len(neg) // 2], neg[len(neg) // 2:]
    for a, b in (("uniprot", "panel"), ("panel", "uniprot")):
        tr_i = np.concatenate([np.where((src == a) & (y == 1))[0], neg_tr])
        te_i = np.concatenate([np.where((src == b) & (y == 1))[0], neg_te])
        sc, m = fit(Xe[tr_i], y[tr_i])
        p = score(sc, m, Xe[te_i])
        ap = average_precision_score(y[te_i], p)
        report[f"train_{a}_test_{b}"] = {"auprc": float(ap),
                                         "n_test_pos": int(y[te_i].sum())}
        print(f"  train on {a:8s} test on {b:8s} AUPRC {ap:.3f} "
              f"({int(y[te_i].sum())} positives)", flush=True)

    # the annotation hold-out: peptides RefSeq calls bacteriocins that no
    # profile recognises and that never entered training
    hold = []
    ph = os.path.join(OUT, "precursor_holdout.csv")
    if os.path.exists(ph):
        hold = list(csv.DictReader(open(ph)))

    # final model on everything, then score the whole panel
    sc, m = fit(Xe, y)
    import joblib
    joblib.dump({"scaler": sc, "model": m, "esm": idx["model"]},
                os.path.join(OUT, "precursor_model.joblib"))
    allp = score(sc, m, emb)
    seq_score = {s: float(p) for s, p in zip(idx["seqs"], allp)}

    # The calling threshold is a false-discovery rate, not a percentile of the
    # training negatives: a model scores its own training negatives far lower
    # than unseen ones, so an in-sample percentile is too permissive. The null
    # is the out-of-fold score of the negatives, which is what an unseen
    # non-precursor looks like.
    null = oof_lm[y == 0]
    scored = sorted({(s_, seq_score[s_]) for s_ in idx["seqs"]}, key=lambda t: -t[1])
    scores_all = np.array([v for _, v in scored])
    q = np.ones(len(scores_all))
    for i, v in enumerate(scores_all):
        q[i] = min(1.0, (null >= v).mean() * len(scores_all) / (i + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    n_q05 = int((q < 0.05).sum())
    thr = float(scores_all[q < 0.05].min()) if n_q05 else 1.0
    in_sample = float(np.quantile([seq_score[r["seq"]] for r in tr
                                   if int(r["label"]) == 0], 0.999))
    thr20 = float(scores_all[q < 0.20].min()) if (q < 0.20).any() else 1.0
    report["threshold"] = {"rule": "empirical-null q value",
                           "score": thr, "n_distinct_called": n_q05,
                           "score_q20": thr20,
                           "n_distinct_called_q20": int((q < 0.20).sum()),
                           "in_sample_999th_percentile_for_comparison": in_sample,
                           "n_distinct_scored": len(scores_all)}
    report["threshold_q20"] = thr20
    report["threshold_999th_percentile_of_negatives"] = thr   # kept for downstream
    with open(os.path.join(OUT, "precursor_qvalues.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["seq", "score", "q"])
        for (s_, v), qq in zip(scored, q):
            if qq < 0.5:
                w.writerow([s_, round(v, 4), round(float(qq), 5)])
    print(f"  calling threshold: score >= {thr:.4f} (q < 0.05), {n_q05} distinct "
          f"sequences of {len(scores_all)}; the in-sample 99.9th percentile would "
          f"have been {in_sample:.3f}", flush=True)

    if hold:
        hs = np.array([seq_score[h["seq"]] for h in hold if h["seq"] in seq_score])
        neg_scores = np.array([seq_score[r["seq"]] for r in tr if int(r["label"]) == 0])
        rec = float((hs >= thr).mean())
        yy = np.r_[np.ones(len(hs)), np.zeros(len(neg_scores))]
        ss = np.r_[hs, neg_scores]
        report["holdout_annotation"] = {
            "n": int(len(hs)), "recall_at_threshold": rec,
            "median_score": float(np.median(hs)),
            "auprc_against_negatives": float(average_precision_score(yy, ss))}
        print(f"  annotation hold-out: {int((hs >= thr).sum())}/{len(hs)} recovered at "
              f"the threshold ({100*rec:.0f}%), AP against the negatives "
              f"{report['holdout_annotation']['auprc_against_negatives']:.3f}", flush=True)

    # map back to the genomes
    dom = collections.defaultdict(dict)
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            dom[f[0]].setdefault(f[1], set()).add(f[3])
    pos_of, anchors = {}, collections.defaultdict(list)
    prod = {}
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a, rep, prot, s, e = f[0], f[1], f[5], int(f[6]), int(f[7])
            pos_of[(a, prot)] = (rep, s, e, f[8], f[2])
            prod[(a, prot)] = f[10]
            if dom.get(a, {}).get(prot, set()) & ANCHOR_DOMS:
                anchors[(a, rep)].append((s, e, prot))

    def nearest_anchor(a, rep, s, e):
        best, who = None, ""
        for (as_, ae, prot) in anchors.get((a, rep), ()):
            d = 0 if (s <= ae and e >= as_) else min(abs(s - ae), abs(as_ - e))
            if best is None or d < best:
                best, who = d, prot
        return best, who

    rows = []
    for f in sorted(os.listdir(FAA)):
        a = f[:-4]
        for prot, s_ in read_fasta(os.path.join(FAA, f)):
            if not (10 <= len(s_) <= 120) or (a, prot) not in pos_of:
                continue
            rep, s, e, strand, reptype = pos_of[(a, prot)]
            d, who = nearest_anchor(a, rep, s, e)
            ds = dom.get(a, {}).get(prot, set())
            rows.append(dict(assembly=a, replicon=rep, replicon_type=reptype,
                             protein=prot, start=s, end=e, strand=strand,
                             aa=len(s_), score=round(seq_score[s_], 4),
                             known_domain=",".join(sorted(ds & PRECURSOR_DOMS)),
                             any_domain=",".join(sorted(ds)),
                             nearest_anchor_bp="" if d is None else d,
                             anchor_protein=who, product=prod[(a, prot)], seq=s_))
    with open(os.path.join(OUT, "precursor_scan.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"  scored {len(rows)} short coding sequences", flush=True)

    # external check: do the high scorers sit beside the machinery?
    # Sequences that were in the training set are dropped from both groups: the
    # negatives were deliberately sampled far from the machinery, so leaving
    # them in would manufacture the very enrichment being tested.
    trained = {r["seq"] for r in tr}
    pool = [r for r in rows if not r["known_domain"] and r["seq"] not in trained]
    def near(r):
        return r["nearest_anchor_bp"] != "" and int(r["nearest_anchor_bp"]) <= NEAR
    hi = [r for r in pool if r["score"] >= thr]
    bg = [r for r in pool if r["score"] < thr]
    a_ = sum(1 for r in hi if near(r)); b_ = len(hi) - a_
    c_ = sum(1 for r in bg if near(r)); d_ = len(bg) - c_
    from scipy import stats
    orr, p = stats.fisher_exact([[a_, b_], [c_, d_]])

    # and a positional null: keep every short protein where it is, move the
    # machinery to random positions on its own replicon, 200 times
    rng = np.random.default_rng(20260918)
    rep_len = {}
    with gzip.open(os.path.join(OUT, "cds.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rep_len[(f[0], f[1])] = int(f[3])
    hi_by_rep = collections.defaultdict(list)
    for r in hi:
        hi_by_rep[(r["assembly"], r["replicon"])].append((r["start"], r["end"]))
    null = []
    for _ in range(200):
        k = 0
        for key, spans in hi_by_rep.items():
            L = rep_len.get(key, 0)
            anc = anchors.get(key, [])
            if not anc or L <= 0:
                continue
            fake = rng.integers(0, L, size=len(anc))
            for (s, e) in spans:
                if any(min(abs(s - (f0 + 1000)), abs(f0 - e)) <= NEAR for f0 in fake):
                    k += 1
        null.append(k)
    null = np.array(null)
    p_perm = float((np.sum(null >= a_) + 1) / (len(null) + 1))
    report["enrichment_near_machinery"] = {
        "high_near": a_, "high_far": b_, "rest_near": c_, "rest_far": d_,
        "odds_ratio": float(orr), "p": float(p),
        "permutation_null_mean": float(null.mean()),
        "permutation_p": p_perm}
    print(f"  high-scoring, no Pfam hit: {len(hi)}; {a_} within {NEAR//1000} kb of "
          f"machinery vs {c_}/{len(bg)} of the rest, OR {orr:.1f}, P = {p:.1e}; "
          f"permutation null {null.mean():.1f}, P = {p_perm:.3f}", flush=True)
    json.dump(report, open(os.path.join(OUT, "precursor_report.json"), "w"), indent=1)
    print("wrote out/precursor_scan.csv and out/precursor_report.json", flush=True)


if __name__ == "__main__":
    main()
