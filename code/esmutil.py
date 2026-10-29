#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""替任意序列取 ESM-2 嵌入（和 s18_embed.py 同一個模型、同樣的 mean-pool）。

面板裡的序列在 out/esm_emb.npy 就有；洗牌序列、突變序列這些「不存在於基因體」的
序列要另外算，算完存進 out/esm_extra.npz 當快取。
"""
import json, os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
MODEL = os.environ.get("ESM_MODEL", "facebook/esm2_t30_150M_UR50D")
CACHE = os.path.join(OUT, "esm_extra.npz")
BATCH_TOKENS = 1200   # WSL 記憶體有限，批次放小一點


def _load_cache():
    if os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=True)
        return list(z["seqs"]), z["emb"]
    return [], np.zeros((0, 640), dtype=np.float16)


def embed(seqs, verbose=True):
    """回傳 {seq: vector}；已經算過的（面板內或快取內）不重算。"""
    idx = json.load(open(os.path.join(OUT, "esm_index.json")))
    panel = {s: i for i, s in enumerate(idx["seqs"])}
    # mmap：26,708 × 640 的矩陣不必整份進記憶體
    panel_emb = np.load(os.path.join(OUT, "esm_emb.npy"), mmap_mode="r")
    cseqs, cemb = _load_cache()
    cache = {s: cemb[i] for i, s in enumerate(cseqs)}

    out, todo = {}, []
    for s in dict.fromkeys(seqs):
        if s in panel:
            out[s] = panel_emb[panel[s]].astype(np.float32)
        elif s in cache:
            out[s] = cache[s].astype(np.float32)
        else:
            todo.append(s)
    if todo:
        import torch
        from transformers import AutoModel, AutoTokenizer
        if verbose:
            print(f"  embedding {len(todo)} new sequences with {MODEL}", flush=True)
        torch.set_num_threads(4)
        tok = AutoTokenizer.from_pretrained(MODEL)
        mod = AutoModel.from_pretrained(MODEL).eval()
        todo = sorted(todo, key=len)
        i = 0
        while i < len(todo):
            n = max(1, BATCH_TOKENS // (len(todo[i]) + 2))
            batch = todo[i:i + n]
            with torch.no_grad():
                b = tok(batch, return_tensors="pt", padding=True)
                h = mod(**b).last_hidden_state
                m = b["attention_mask"].clone()
                m[:, 0] = 0
                j = m.sum(1) - 1
                m[torch.arange(len(batch)), j + 1] = 0
                m = m.unsqueeze(-1).float()
                pooled = (h * m).sum(1) / m.sum(1).clamp(min=1)
            for s, v in zip(batch, pooled.numpy().astype(np.float16)):
                cache[s] = v
                out[s] = v.astype(np.float32)
            i += len(batch)
        ss = list(cache)
        np.savez(CACHE, seqs=np.array(ss, dtype=object),
                 emb=np.stack([cache[s] for s in ss]))
    return out
