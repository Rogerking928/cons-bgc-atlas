#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""跟真正的同類工具比：DeepBGC（2026-09-20 審閱意見）。

（英文原註保留如下）
An external tool on the same loci: DeepBGC.

The detector in this paper is compared against three baselines of our own,
which tells a reader nothing about how it stands against the field. DeepBGC is
the one published deep-learning BGC tool that installs and runs offline here
(0.1.10, models 0.1.0, Pfam-A 31.0, through micromamba because the system
Python is externally managed). It is run on a 10 kb window either side of one
representative cluster per family, plus control windows drawn at random from
the same genomes, and asked two things: does it call a BGC there, and does it
call it a RiPP.
"""
import collections, csv, gzip, json, os, random, subprocess, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FNA = os.path.join(HERE, "data", "fna")
WORK = os.path.join(OUT, "deepbgc")
MM = os.path.expanduser("~/.local/micromamba/bin/micromamba")
ENV = "bgc38"
PAD = 10000
FAMS = ["F001", "F003", "F004", "F005", "F006", "F007", "F009", "F011"]
N_CONTROL = 8


def load(acc):
    seqs, n, b = {}, None, []
    with gzip.open(os.path.join(FNA, acc + ".fna.gz"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    seqs[n] = "".join(b)
                n, b = line[1:].split()[0], []
            else:
                b.append(line.strip().upper())
    if n:
        seqs[n] = "".join(b)
    return seqs


def main():
    os.makedirs(WORK, exist_ok=True)
    fam = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    by_fam = collections.defaultdict(list)
    for c in cl:
        by_fam[fam[c["cluster"]]].append(c)

    jobs = []
    rng = random.Random(20260919)
    for f in FAMS:
        cand = sorted(by_fam[f], key=lambda c: int(c["span"]))
        if not cand:
            continue
        c = cand[len(cand) // 2]
        seqs = load(c["assembly"])
        s = seqs.get(c["replicon"], "")
        lo, hi = max(0, int(c["start"]) - PAD), min(len(s), int(c["end"]) + PAD)
        p = os.path.join(WORK, f"{f}_{c['assembly']}.fasta")
        with open(p, "w") as fh:
            fh.write(f">{f}_{c['assembly']}_{c['replicon']}_{lo}_{hi}\n{s[lo:hi]}\n")
        jobs.append((f, c["assembly"], p, int(c["start"]) - lo, int(c["end"]) - lo))
    # controls: windows of the same size with no called cluster on them
    used = {(j[1]) for j in jobs}
    for k, acc in enumerate(sorted(used)[:N_CONTROL], 1):
        seqs = load(acc)
        rep, s = max(seqs.items(), key=lambda kv: len(kv[1]))
        starts = [int(c["start"]) for c in cl if c["assembly"] == acc]
        for _ in range(50):
            lo = rng.randrange(0, max(1, len(s) - 40000))
            if all(abs(lo - x) > 60000 for x in starts):
                break
        p = os.path.join(WORK, f"control{k}_{acc}.fasta")
        with open(p, "w") as fh:
            fh.write(f">control{k}_{acc}_{rep}_{lo}\n{s[lo:lo + 40000]}\n")
        jobs.append((f"control{k}", acc, p, None, None))

    rows = []
    for tag, acc, path, lo, hi in jobs:
        outdir = os.path.join(WORK, tag)
        if not os.path.exists(os.path.join(outdir, os.path.basename(path)[:-6] + ".bgc.tsv")):
            r = subprocess.run([MM, "run", "-n", ENV, "deepbgc", "pipeline",
                                "--output", outdir, path],
                               capture_output=True, text=True,
                               env={**os.environ,
                                    "MAMBA_ROOT_PREFIX": os.path.expanduser("~/.local/micromamba")})
            if r.returncode != 0:
                print(f"  {tag}: failed\n{r.stderr[-600:]}", flush=True)
                rows.append(dict(region=tag, assembly=acc, status="failed",
                                 bgcs=0, ripp=0, best_score="", classes=""))
                continue
        hits = []
        for f in os.listdir(outdir):
            if f.endswith(".bgc.tsv"):
                hits = list(csv.DictReader(open(os.path.join(outdir, f)), delimiter="\t"))
        ripp = [h for h in hits
                if float(h.get("RiPP", 0) or 0) >= 0.5
                or "ripp" in (h.get("product_class", "") or "").lower()]
        rows.append(dict(region=tag, assembly=acc, status="ok", bgcs=len(hits),
                         ripp=len(ripp),
                         best_score=max((h.get("deepbgc_score", "") for h in hits),
                                        default=""),
                         classes=";".join(sorted({h.get("product_class", "") for h in hits}))))
        print(f"  {tag}: {len(hits)} BGC calls, {len(ripp)} RiPP", flush=True)

    with open(os.path.join(OUT, "deepbgc_results.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print("-> out/deepbgc_results.csv", flush=True)


if __name__ == "__main__":
    main()
