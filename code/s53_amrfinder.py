#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run AMRFinderPlus itself over the panel, in place of k-mer containment.

審稿人（AoM #42 Reviewer 2，2026-10-28）指出用 10-residue word containment 當
抗藥基因呼叫、而不用策劃好的每基因門檻，是技術缺口。這支腳本改跑程式本身。

參考資料庫版本與 s8 當初用的 AMRProt.fa 同一版（2026-08-07.1），所以兩者的差異
是「containment vs 策劃門檻」，不是資料庫改版。

Protein mode with the RefSeq GFF, so every call carries a replicon and
coordinates and can go straight into the distance analysis.
Output: out/amrfinder/<assembly>.tsv (one per genome) + out/amrfinder.log
"""
import glob, gzip, os, shutil, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(HERE, "data", "faa")
GFF = os.path.join(HERE, "data", "gff")
OUT = os.path.join(HERE, "out")
DEST = os.path.join(OUT, "amrfinder")

ENV = os.path.expanduser("~/micromamba/envs/amrf")
AMRFINDER = os.path.join(ENV, "bin", "amrfinder")
DB = os.path.join(ENV, "share", "amrfinderplus", "data", "2026-08-07.1")

# 每株給 2 個執行緒、同時 6 株；RAM 只有 11 GB，blastp 同時太多會換頁。
PER_JOB_THREADS = 2
JOBS = 6


def run_one(acc):
    faa = os.path.join(FAA, acc + ".faa")
    gffgz = os.path.join(GFF, acc + ".gff.gz")
    dest = os.path.join(DEST, acc + ".tsv")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return acc, "cached", ""
    if not (os.path.exists(faa) and os.path.exists(gffgz)):
        return acc, "missing-input", ""
    tmpd = tempfile.mkdtemp(prefix="amrf_" + acc + "_")
    try:
        gff = os.path.join(tmpd, acc + ".gff")
        with gzip.open(gffgz, "rb") as fi, open(gff, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        cmd = [AMRFINDER, "-p", faa, "-g", gff, "-a", "pgap",
               "--plus", "-d", DB, "--threads", str(PER_JOB_THREADS),
               "-o", dest + ".part"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            return acc, "fail", (p.stderr or "")[-600:]
        os.replace(dest + ".part", dest)
        return acc, "ok", ""
    except Exception as exc:                                  # noqa: BLE001
        return acc, "error", repr(exc)[:300]
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def main():
    os.makedirs(DEST, exist_ok=True)
    accs = sorted(os.path.basename(p)[:-4]
                  for p in glob.glob(os.path.join(FAA, "*.faa")))
    log = open(os.path.join(OUT, "amrfinder.log"), "a", buffering=1)
    print(f"{len(accs)} proteomes, db={os.path.basename(DB)}", file=log)
    done = bad = 0
    with ProcessPoolExecutor(JOBS) as ex:
        futs = {ex.submit(run_one, a): a for a in accs}
        for fut in as_completed(futs):
            acc, status, msg = fut.result()
            done += 1
            if status not in ("ok", "cached"):
                bad += 1
                print(f"  {acc}\t{status}\t{msg}", file=log)
            if done % 25 == 0:
                print(f"  {done}/{len(accs)} ({bad} bad)", file=log)
    print(f"done {done}/{len(accs)}, {bad} bad", file=log)
    print(f"done {done}, {bad} bad", file=sys.stderr)


if __name__ == "__main__":
    main()
