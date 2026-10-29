#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Embed every short protein in the panel with a protein language model.

ESM-2 (150M parameters) on CPU, mean-pooled over the real residues. The panel
holds 382k short coding sequences but only about 23k distinct ones, so the
embedding is computed per distinct sequence and reused.
"""
import csv, json, os, sys, time

import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FAA = os.path.join(HERE, "data", "faa")
MODEL = os.environ.get("ESM_MODEL", "facebook/esm2_t30_150M_UR50D")
MINLEN, MAXLEN = 10, 120
BATCH_TOKENS = 3000


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


def main():
    seqs = set()
    for f in sorted(os.listdir(FAA)):
        for _, s in read_fasta(os.path.join(FAA, f)):
            if MINLEN <= len(s) <= MAXLEN:
                seqs.add(s)
    for r in csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))):
        seqs.add(r["seq"])
    seqs = sorted(seqs, key=len)
    print(f"{len(seqs)} distinct sequences to embed with {MODEL}", flush=True)

    torch.set_num_threads(4)
    tok = AutoTokenizer.from_pretrained(MODEL)
    mod = AutoModel.from_pretrained(MODEL).eval()
    dim = mod.config.hidden_size
    emb = np.zeros((len(seqs), dim), dtype=np.float16)

    i, t0 = 0, time.time()
    while i < len(seqs):
        n = max(1, BATCH_TOKENS // (len(seqs[i]) + 2))
        batch = seqs[i:i + n]
        with torch.no_grad():
            b = tok(batch, return_tensors="pt", padding=True)
            h = mod(**b).last_hidden_state
            m = b["attention_mask"].clone()
            m[:, 0] = 0                                   # drop the CLS token
            idx = m.sum(1) - 1
            m[torch.arange(len(batch)), idx + 1] = 0      # and the EOS token
            m = m.unsqueeze(-1).float()
            pooled = (h * m).sum(1) / m.sum(1).clamp(min=1)
        emb[i:i + len(batch)] = pooled.numpy().astype(np.float16)
        i += len(batch)
        if (i // 2000) != ((i - len(batch)) // 2000):
            el = time.time() - t0
            print(f"  {i}/{len(seqs)}  {el:.0f}s elapsed, "
                  f"{el/max(i,1)*(len(seqs)-i):.0f}s left", flush=True)
    np.save(os.path.join(OUT, "esm_emb.npy"), emb)
    json.dump({"model": MODEL, "dim": dim, "seqs": seqs},
              open(os.path.join(OUT, "esm_index.json"), "w"))
    print(f"wrote out/esm_emb.npy {emb.shape} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
