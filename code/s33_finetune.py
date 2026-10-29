#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fine-tune the language model instead of freezing it, then ask what it looks at.

LoRA adapters (rank 8) are attached to the query and value projections of the
last four transformer blocks of ESM-2 35M; everything else stays frozen. The
comparison is against the same model used frozen, on the same folds, so the
number that comes out is what fine-tuning is worth on this task.

Attribution is integrated gradients from a padding baseline, summed over the
embedding dimension, which gives one number per residue. The question it
answers is whether the model rests on the double-glycine leader motif that the
C39 peptidase cuts - a feature nobody told it about.
"""
import collections, csv, json, os, random, re, sys, time

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold
from transformers import AutoModel, AutoTokenizer

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
NAME = "facebook/esm2_t12_35M_UR50D"
N_NEG = 3000
EPOCHS = 2
LR = 1e-3
RANK = 8
LAST = 4
BATCH = 8
SEED = 20260919


class LoRALinear(nn.Module):
    def __init__(self, base, r=RANK):
        super().__init__()
        self.base = base
        for p in self.base.parameters():
            p.requires_grad = False
        self.A = nn.Parameter(torch.randn(r, base.in_features) * 0.01)
        self.B = nn.Parameter(torch.zeros(base.out_features, r))

    def forward(self, x):
        return self.base(x) + (x @ self.A.T) @ self.B.T


class Detector(nn.Module):
    def __init__(self):
        super().__init__()
        self.esm = AutoModel.from_pretrained(NAME)
        for p in self.esm.parameters():
            p.requires_grad = False
        layers = self.esm.encoder.layer
        for blk in layers[-LAST:]:
            att = blk.attention.self
            att.query = LoRALinear(att.query)
            att.value = LoRALinear(att.value)
        self.head = nn.Linear(self.esm.config.hidden_size, 1)

    def pooled(self, ids, mask, inputs_embeds=None):
        out = self.esm(input_ids=None if inputs_embeds is not None else ids,
                       attention_mask=mask, inputs_embeds=inputs_embeds)
        h = out.last_hidden_state
        m = mask.clone()
        m[:, 0] = 0
        idx = m.sum(1) - 1
        m[torch.arange(len(ids if ids is not None else inputs_embeds)), idx + 1] = 0
        m = m.unsqueeze(-1).float()
        return (h * m).sum(1) / m.sum(1).clamp(min=1)

    def forward(self, ids, mask, inputs_embeds=None):
        return self.head(self.pooled(ids, mask, inputs_embeds)).squeeze(-1)


def batches(seqs, labels, tok, shuffle=True, rng=None):
    idx = list(range(len(seqs)))
    if shuffle:
        rng.shuffle(idx)
    for i in range(0, len(idx), BATCH):
        j = idx[i:i + BATCH]
        b = tok([seqs[k] for k in j], return_tensors="pt", padding=True)
        yield b["input_ids"], b["attention_mask"], torch.tensor(
            [float(labels[k]) for k in j]), j


def train(seqs, labels, tok, log=""):
    torch.manual_seed(SEED)
    model = Detector()
    params = [p for p in model.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(params, lr=LR)
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(
        max(1.0, (len(labels) - sum(labels)) / max(1, sum(labels))) ** 0.5))
    rng = random.Random(SEED)
    model.train()
    for ep in range(EPOCHS):
        t0, n = time.time(), 0
        for ids, mask, y, _ in batches(seqs, labels, tok, True, rng):
            opt.zero_grad()
            loss = lossf(model(ids, mask), y)
            loss.backward()
            opt.step()
            n += len(y)
            if n % 800 < BATCH:
                print(f"    {log} epoch {ep + 1}: {n}/{len(seqs)} "
                      f"{time.time() - t0:.0f}s loss {loss.item():.4f}", flush=True)
    model.eval()
    return model


@torch.no_grad()
def predict(model, seqs, tok):
    out = np.zeros(len(seqs))
    for ids, mask, _, j in batches(seqs, [0] * len(seqs), tok, False):
        out[j] = torch.sigmoid(model(ids, mask)).numpy()
    return out


def integrated_gradients(model, seq, tok, steps=24):
    b = tok([seq], return_tensors="pt")
    ids, mask = b["input_ids"], b["attention_mask"]
    emb_layer = model.esm.embeddings.word_embeddings
    x = emb_layer(ids)
    base = torch.zeros_like(x)
    total = torch.zeros_like(x)
    for a in np.linspace(1.0 / steps, 1.0, steps):
        xi = (base + a * (x - base)).detach().requires_grad_(True)
        score = model(None, mask, inputs_embeds=xi)
        score.backward()
        total += xi.grad
    attr = ((x - base) * total / steps).sum(-1)[0].detach().numpy()
    return attr[1:1 + len(seq)]


def main():
    tok = AutoTokenizer.from_pretrained(NAME)
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    pos = [r for r in tr if int(r["label"]) == 1]
    neg = [r for r in tr if int(r["label"]) == 0]
    random.Random(SEED).shuffle(neg)
    sub = pos + neg[:N_NEG]
    seqs = [r["seq"] for r in sub]
    y = np.array([int(r["label"]) for r in sub])
    g = np.array([int(r["cluster"]) for r in sub])
    print(f"{len(sub)} sequences ({y.sum()} positives) for fine-tuning", flush=True)

    oof = np.zeros(len(sub))
    for k, (a, b) in enumerate(GroupKFold(n_splits=5).split(seqs, y, groups=g), 1):
        m = train([seqs[i] for i in a], [y[i] for i in a], tok, log=f"fold {k}")
        oof[b] = predict(m, [seqs[i] for i in b], tok)
        print(f"  fold {k}: AP {average_precision_score(y[b], oof[b]):.3f}", flush=True)
    ap = float(average_precision_score(y, oof))
    print(f"fine-tuned CV AP {ap:.3f}", flush=True)

    model = train(seqs, list(y), tok, log="full")
    torch.save({"state": model.state_dict()}, os.path.join(OUT, "finetuned_35M.pt"))

    # attribution: positives, candidates and negatives
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    t20 = rep["threshold"].get("score_q20", rep["threshold"]["score"])
    scan = list(csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv"))))
    seen, sets = set(), collections.defaultdict(list)
    for r in scan:
        if r["seq"] in seen:
            continue
        seen.add(r["seq"])
        if r["known_domain"]:
            sets["profile-detected precursor"].append(r["seq"])
        elif float(r["score"]) >= t20:
            sets["candidate"].append(r["seq"])
    sets["negative"] = [r["seq"] for r in neg[:40]]
    rows = []
    for grp, ss in sets.items():
        for s in ss[:60]:
            attr = integrated_gradients(model, s, tok)
            if attr.size < 12:
                continue
            gg = [m.start() for m in re.finditer(r"G[GA]", s)]
            top = set(np.argsort(attr)[-5:].tolist())
            near_gg = any(any(abs(t - p) <= 2 for p in gg) for t in top)
            rows.append(dict(group=grp, seq=s, length=len(s),
                             top_positions=";".join(str(t) for t in sorted(top)),
                             gg_positions=";".join(str(p) for p in gg),
                             top_near_gg=int(near_gg),
                             attribution=";".join(f"{v:.4f}" for v in attr)))
        print(f"  attributed {grp}: {sum(1 for r in rows if r['group'] == grp)}",
              flush=True)
    with open(os.path.join(OUT, "saliency.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    summary = {"finetuned_cv_ap": ap,
               "near_gg": {grp: [sum(r["top_near_gg"] for r in rows if r["group"] == grp),
                                 sum(1 for r in rows if r["group"] == grp)]
                           for grp in sets}}
    json.dump(summary, open(os.path.join(OUT, "finetune_summary.json"), "w"), indent=1)
    for grp, (k, n) in summary["near_gg"].items():
        print(f"  {grp:28s} top-attributed residue next to a GG motif in {k}/{n}",
              flush=True)
    print("-> out/saliency.csv, out/finetune_summary.json", flush=True)


if __name__ == "__main__":
    main()
