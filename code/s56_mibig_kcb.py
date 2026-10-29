#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Known-cluster assignment against MIBiG, via antiSMASH knownclusterblast.

審稿人（AoM #42 Reviewer 1，2026-10-28）：「It is unclear how known BGCs are
determined, is this via MIBiG - if so cite it!」

原稿的家族命名來自 RefSeq product name 與 UniProt（§2.3），**不是** MIBiG。
這支腳本補上一個獨立的、以 MIBiG 為參照的指認：對每株再跑一次 antiSMASH
（同 s54 的設定：RefSeq GFF、--minimal），只多開 --cb-knownclusters，
把每個 region 對 MIBiG 的最佳命中記下來。s55 的 region 座標與這裡相同
（同一版 antiSMASH、同一份輸入），所以可以用 (assembly, replicon, region_number)
直接接回去。

similarity 是 antiSMASH 的定義：MIBiG 參考簇的基因裡，有幾成在 region 內
找得到 BLAST 命中（不是序列一致度）。

Output: out/antismash_kcb/<assembly>.kcb.tsv（每個 region 一列）
        out/antismash_kcb.log
"""
import csv, gzip, json, os, shutil, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s54_antismash import (ANTISMASH, FNA, GFF, OUT, RUNENV, panel,  # noqa: E402
                           write_filtered_gff)

DEST = os.path.join(OUT, "antismash_kcb")
PER_JOB_CPUS = 3
JOBS = 5          # 單株約 300 MB RAM、20 秒；11 GB RAM 撐得起 5 個
TIMEOUT = 3600
COLS = ["assembly", "replicon", "region_number", "mibig_accession",
        "mibig_description", "mibig_type", "similarity", "core_gene_hits",
        "hits", "mibig_version"]


def parse_kcb(acc, path):
    with open(path) as fh:
        doc = json.load(fh)
    rows = []
    for rec in doc.get("records", []):
        cb = rec.get("modules", {}).get("antismash.modules.clusterblast") or {}
        kc = cb.get("knowncluster") or {}
        ver = kc.get("data_version", "")
        for res in kc.get("results", []):
            rank = res.get("ranking") or []
            if rank:
                ref, score = rank[0]
                rows.append([acc, rec.get("id", ""), res.get("region_number"),
                             ref.get("accession", ""), ref.get("description", ""),
                             ref.get("cluster_type", ""), score.get("similarity", ""),
                             score.get("core_gene_hits", ""), score.get("hits", ""),
                             ver])
            else:
                rows.append([acc, rec.get("id", ""), res.get("region_number"),
                             "", "", "", "", "", "", ver])
    return rows


def run_one(acc):
    fnagz = os.path.join(FNA, acc + ".fna.gz")
    gffgz = os.path.join(GFF, acc + ".gff.gz")
    dest = os.path.join(DEST, acc + ".kcb.tsv")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return acc, "cached", ""
    tmpd = tempfile.mkdtemp(prefix="kcb_" + acc + "_")
    try:
        fna = os.path.join(tmpd, acc + ".fna")
        gff = os.path.join(tmpd, acc + ".gff")
        with gzip.open(fnagz, "rb") as fi, open(fna, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        write_filtered_gff(gffgz, gff)
        work = os.path.join(tmpd, "work")
        cmd = [ANTISMASH, "--minimal", "--cb-knownclusters",
               "--genefinding-tool", "none", "--genefinding-gff3", gff,
               "--cpus", str(PER_JOB_CPUS), "--output-dir", work,
               "--output-basename", acc,
               "--logfile", os.path.join(tmpd, "as.log"), fna]
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=TIMEOUT, env=RUNENV)
        res = os.path.join(work, acc + ".json")
        if p.returncode != 0 or not os.path.exists(res):
            return acc, "fail", (p.stderr or "")[-700:]
        rows = parse_kcb(acc, res)
        with open(dest + ".part", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(COLS)
            w.writerows(rows)
        os.replace(dest + ".part", dest)
        return acc, "ok", f"{len(rows)} regions"
    except subprocess.TimeoutExpired:
        return acc, "timeout", ""
    except Exception as exc:                                  # noqa: BLE001
        return acc, "error", repr(exc)[:300]
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def main():
    os.makedirs(DEST, exist_ok=True)
    accs = panel()
    log = open(os.path.join(OUT, "antismash_kcb.log"), "a", buffering=1)
    print(f"{len(accs)} genomes", file=log)
    done = bad = 0
    with ProcessPoolExecutor(JOBS) as ex:
        futs = {ex.submit(run_one, a): a for a in accs}
        for fut in as_completed(futs):
            acc, status, msg = fut.result()
            done += 1
            if status not in ("ok", "cached"):
                bad += 1
                print(f"  {acc}\t{status}\t{msg}", file=log)
            if done % 10 == 0:
                print(f"  {done}/{len(accs)} ({bad} bad)", file=log)
    print(f"done {done}/{len(accs)}, {bad} bad", file=log)


if __name__ == "__main__":
    main()
