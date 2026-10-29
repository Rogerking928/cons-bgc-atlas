#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Does a bigger language model help?

The same training set, the same folds and the same classifier, over the ESM-2
size ladder. The 150M embeddings already exist from the panel-wide scan; the
others are computed here for the training sequences only, which is what the
comparison needs.
"""
import csv, json, os, sys, time

import numpy as np
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from transformers import AutoModel, AutoTokenizer

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
MODELS = [("esm2_t6_8M_UR50D", "8M"), ("esm2_t12_35M_UR50D", "35M"),
          ("esm2_t30_150M_UR50D", "150M")]
BIG = ("esm2_t33_650M_UR50D", "650M")
BATCH_TOKENS = 3000


def embed(name, seqs):
    torch.set_num_threads(4)
    tok = AutoTokenizer.from_pretrained("facebook/" + name)
    mod = AutoModel.from_pretrained("facebook/" + name).eval()
    out = np.zeros((len(seqs), mod.config.hidden_size), dtype=np.float32)
    i, t0 = 0, time.time()
    order = np.argsort([len(s) for s in seqs])
    srt = [seqs[j] for j in order]
    buf = np.zeros_like(out)
    while i < len(srt):
        n = max(1, BATCH_TOKENS // (len(srt[i]) + 2))
        batch = srt[i:i + n]
        with torch.no_grad():
            b = tok(batch, return_tensors="pt", padding=True)
            h = mod(**b).last_hidden_state
            m = b["attention_mask"].clone()
            m[:, 0] = 0
            idx = m.sum(1) - 1
            m[torch.arange(len(batch)), idx + 1] = 0
            m = m.unsqueeze(-1).float()
            buf[i:i + len(batch)] = ((h * m).sum(1) / m.sum(1).clamp(min=1)).numpy()
        i += len(batch)
        if (i // 4000) != ((i - len(batch)) // 4000):
            print(f"    {name}: {i}/{len(srt)}  {time.time() - t0:.0f}s", flush=True)
    out[order] = buf
    return out


def cv_ap(X, y, g):
    oof = np.zeros(len(y))
    for a, b in GroupKFold(n_splits=5).split(X, y, groups=g):
        sc = StandardScaler().fit(X[a])
        m = LogisticRegression(max_iter=4000, C=0.05, class_weight="balanced")
        m.fit(sc.transform(X[a]), y[a])
        oof[b] = m.predict_proba(sc.transform(X[b]))[:, 1]
    return float(average_precision_score(y, oof)), oof


def main():
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    y = np.array([int(r["label"]) for r in tr])
    g = np.array([int(r["cluster"]) for r in tr])
    seqs = [r["seq"] for r in tr]
    res = {}
    p = os.path.join(OUT, "ladder.json")
    if os.path.exists(p):
        res = json.load(open(p))

    idx = json.load(open(os.path.join(OUT, "esm_index.json")))
    emb150 = np.load(os.path.join(OUT, "esm_emb.npy")).astype(np.float32)
    row = {s: i for i, s in enumerate(idx["seqs"])}

    todo = list(MODELS) + ([BIG] if os.environ.get("LADDER_BIG") else [])
    for name, tag in todo:
        if tag in res:
            continue
        t0 = time.time()
        if tag == "150M":
            X = np.stack([emb150[row[s]] for s in seqs])
        else:
            X = embed(name, seqs)
            np.save(os.path.join(OUT, f"emb_train_{tag}.npy"), X.astype(np.float16))
        ap, _ = cv_ap(X, y, g)
        res[tag] = {"model": name, "dim": int(X.shape[1]), "cv_ap": ap,
                    "seconds": round(time.time() - t0, 1)}
        json.dump(res, open(p, "w"), indent=1)
        print(f"  {tag:5s} dim {X.shape[1]:4d}  CV AP {ap:.3f}  "
              f"({time.time() - t0:.0f}s)", flush=True)
    print("-> out/ladder.json")


if __name__ == "__main__":
    main()
