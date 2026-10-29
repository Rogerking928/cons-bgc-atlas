#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""hmmsearch every proteome against the Pfam profile set, at the curated
gathering threshold. One domtblout per assembly in out/dom/.

Usage: python3 s4_hmmscan.py [n_parallel]
"""
import glob, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(HERE, "data", "faa")
DOM = os.path.join(HERE, "out", "dom")
HMM = os.path.join(HERE, "data", "hmm", "library.hmm")
os.makedirs(DOM, exist_ok=True)

N = int(sys.argv[1]) if len(sys.argv) > 1 else 8


def scan(p):
    a = os.path.basename(p)[:-4]
    out = os.path.join(DOM, a + ".domtbl")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return a, "skip"
    tmp = out + ".tmp"
    r = subprocess.run(["hmmsearch", "--cut_ga", "--cpu", "1", "--noali",
                        "--domtblout", tmp, HMM, p],
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if r.returncode != 0:
        return a, "FAIL " + r.stderr.decode()[:200]
    os.replace(tmp, out)
    return a, "ok"


if __name__ == "__main__":
    faas = sorted(glob.glob(os.path.join(FAA, "*.faa")))
    print(f"{len(faas)} proteomes, {N} workers", flush=True)
    done = 0
    with ThreadPoolExecutor(N) as ex:
        for a, st in ex.map(scan, faas):
            done += 1
            if st != "skip":
                if done % 25 == 0 or st.startswith("FAIL"):
                    print(done, a, st, flush=True)
    print("done", len(os.listdir(DOM)), flush=True)
