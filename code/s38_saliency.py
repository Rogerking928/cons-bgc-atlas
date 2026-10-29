#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Attribution for the complete candidate set, and the aggregate profile.

The first pass capped each group at sixty sequences, which left the n in the
figure undefined. This one takes every q < 0.20 candidate and every
profile-detected precursor, and a fixed random sample of negatives, and adds
the aggregate that a single example cannot give: every sequence aligned on its
double-glycine motif and the median attribution taken across the group.
"""
import collections, csv, json, os, random, re, sys

import numpy as np
import torch
from transformers import AutoTokenizer

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s33_finetune import Detector, NAME, integrated_gradients  # noqa: E402

WIN = 30
N_NEG = 60


def main():
    tok = AutoTokenizer.from_pretrained(NAME)
    model = Detector()
    model.load_state_dict(torch.load(os.path.join(OUT, "finetuned_35M.pt"),
                                     map_location="cpu")["state"])
    model.eval()

    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    t05 = rep["threshold"]["score"]
    t20 = rep["threshold"].get("score_q20", t05)
    scan = list(csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv"))))
    seen, sets = set(), collections.defaultdict(list)
    for r in scan:
        if r["seq"] in seen:
            continue
        seen.add(r["seq"])
        if r["known_domain"]:
            sets["profile-detected precursor"].append(r["seq"])
        elif float(r["score"]) >= t20:
            sets["candidate, q < 0.20"].append(r["seq"])
            if float(r["score"]) >= t05:
                sets["candidate, q < 0.05"].append(r["seq"])
    tr = list(csv.DictReader(open(os.path.join(OUT, "precursor_trainset.csv"))))
    neg = [r["seq"] for r in tr if int(r["label"]) == 0]
    random.Random(20260919).shuffle(neg)
    sets["negative"] = neg[:N_NEG]
    print({k: len(v) for k, v in sets.items()}, flush=True)

    rows, prof = [], collections.defaultdict(list)
    for grp, ss in sets.items():
        for n, s in enumerate(ss, 1):
            attr = integrated_gradients(model, s, tok)
            if attr.size < 12:
                continue
            gg = [m.start() for m in re.finditer(r"G[GA]", s)]
            top = set(np.argsort(attr)[-5:].tolist())
            near = any(any(abs(t - p) <= 2 for p in gg) for t in top)
            rows.append(dict(group=grp, seq=s, length=len(s),
                             top_positions=";".join(str(t) for t in sorted(top)),
                             gg_positions=";".join(str(p) for p in gg),
                             top_near_gg=int(near),
                             attribution=";".join(f"{v:.4f}" for v in attr)))
            if gg:
                # align on the motif the model attributes most, falling back to
                # the first motif in the leader
                best = max(gg, key=lambda p: attr[max(0, p - 1):p + 3].max())
                z = (attr - attr.mean()) / (attr.std() + 1e-9)
                for d in range(-WIN, WIN + 1):
                    j = best + d
                    if 0 <= j < len(z):
                        prof[grp].append((d, float(z[j])))
            if n % 20 == 0:
                print(f"  {grp}: {n}/{len(ss)}", flush=True)

    with open(os.path.join(OUT, "saliency.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    meta = {}
    for grp, pts in prof.items():
        by = collections.defaultdict(list)
        for d, v in pts:
            by[d].append(v)
        meta[grp] = {str(d): [float(np.median(by[d])), len(by[d])]
                     for d in sorted(by)}
    summary = {"sets": {k: len(v) for k, v in sets.items()},
               "near_gg": {g: [sum(r["top_near_gg"] for r in rows if r["group"] == g),
                               sum(1 for r in rows if r["group"] == g)] for g in sets},
               "meta_profile": meta,
               "finetuned_cv_ap": json.load(open(os.path.join(
                   OUT, "finetune_summary.json")))["finetuned_cv_ap"]}
    json.dump(summary, open(os.path.join(OUT, "finetune_summary.json"), "w"), indent=1)
    for g, (k, n) in summary["near_gg"].items():
        print(f"  {g:30s} top-attributed residue next to a GG motif in {k}/{n}",
              flush=True)
    print("-> out/saliency.csv, out/finetune_summary.json", flush=True)


if __name__ == "__main__":
    main()
