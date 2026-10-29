#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run antiSMASH over the panel, as the detector comparison.

審稿人（AoM #42 Reviewer 1，2026-10-28）：稿件只比了 DeepBGC，沒有比 antiSMASH，
而 antiSMASH 是這個領域現在的標準。原稿在 §2.11 只把 antiSMASH 當「規則的來源」
引用，沒有當偵測器跑過，也沒有說明為什麼不跑。這支腳本把它跑起來。

兩個刻意的選擇，為了讓比較只差在「偵測器」這一項：

1. **餵 RefSeq 的 GFF**（--genefinding-gff3），不讓 antiSMASH 自己用 prodigal
   找基因。否則差異裡會混進基因預測的不同，而本圖譜用的是 RefSeq 的 CDS。
2. **--minimal**：只做基因簇偵測，不跑 clusterblast／知名簇比對等分析模組。
   要比的是「哪一段算一個基因簇」，不是註解的豐富程度；也省掉大半執行時間。

Output: out/antismash/<assembly>/ (regions.js + .gbk) + out/antismash.log
"""
import collections, csv, glob, gzip, json, os, re, shutil, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FNA = os.path.join(HERE, "data", "fna")
GFF = os.path.join(HERE, "data", "gff")
OUT = os.path.join(HERE, "out")
DEST = os.path.join(OUT, "antismash")

ENV = os.path.expanduser("~/micromamba/envs/as")
ANTISMASH = os.path.join(ENV, "bin", "antismash")
# antiSMASH 會從 PATH 找自己要用的 hmmscan／hmmpress／blastp／diamond。
# 用絕對路徑呼叫而不把環境的 bin 放進 PATH，它會回報
# "Failed to locate executable for 'hmmscan'" 而不是找不到資料庫 —— 踩過一次。
RUNENV = dict(os.environ, PATH=os.path.join(ENV, "bin") + os.pathsep
              + os.environ.get("PATH", ""))

# RAM 11 GB 是瓶頸，不是核心數：hmmscan 同時跑太多會換頁，反而更慢。
PER_JOB_CPUS = 4
JOBS = 3
TIMEOUT = 3600


def panel():
    """只跑收案的 554 株；data/fna 另外還有 draft 與 outgroup，不在這一步。

    收案名單以 out/amr_per_genome.csv 為準，它是每株一列的全panel 表。
    """
    with open(os.path.join(OUT, "amr_per_genome.csv")) as fh:
        return sorted({r["assembly"] for r in csv.DictReader(fh)})


REGION_COLS = ["assembly", "replicon", "region_number", "start", "end",
               "strand", "product", "contig_edge", "rules"]
LOC_RE = re.compile(r"\[(<?\d+):(>?\d+)\]\(([+-])\)")
# antiSMASH 的 secmet.Feature 要求 feature type 長度 < 16（16 字在 GenBank 格式裡
# 會和座標黏在一起），而 RefSeq 的 GFF 有 `pseudogenic_tRNA`（剛好 16 字），
# 整個基因體會因此被拒。這裡在送進去之前按長度濾掉，不是只針對那一個名字；
# 被濾掉的都是 tRNA／riboswitch 這類非 CDS 註解，不影響以 CDS 為基礎的簇偵測。
MAX_FEATURE_TYPE = 15


def write_filtered_gff(src_gz, dest):
    """解壓 GFF，順手濾掉 type 名稱過長的 feature 行。回傳被丟掉的 type 計數。"""
    dropped = collections.Counter()
    with gzip.open(src_gz, "rt") as fi, open(dest, "w") as fo:
        for line in fi:
            if line.startswith("#"):
                fo.write(line)
                continue
            f = line.split("\t")
            if len(f) >= 9 and len(f[2]) > MAX_FEATURE_TYPE:
                dropped[f[2]] += 1
                continue
            fo.write(line)
    return dropped


def parse_regions(acc, path):
    """把 antiSMASH 的 <acc>.json 讀成每個 region 一列。

    region feature 的 location 長這樣：'[181833:203934](+)'，是 0-based 半開區間；
    這裡轉成 1-based 起點，和 out/clusters.tsv（來自 GFF）同一套座標。
    """
    with open(path) as fh:
        doc = json.load(fh)
    rows = []
    for rec in doc.get("records", []):
        rid = rec.get("id", "")
        for feat in rec.get("features", []):
            if feat.get("type") != "region":
                continue
            m = LOC_RE.search(str(feat.get("location", "")))
            if not m:
                continue
            start = int(m.group(1).lstrip("<")) + 1
            end = int(m.group(2).lstrip(">"))
            q = feat.get("qualifiers", {})
            rows.append([acc, rid, (q.get("region_number") or [""])[0],
                         start, end, m.group(3),
                         ";".join(q.get("product") or []),
                         (q.get("contig_edge") or [""])[0],
                         ";".join(q.get("rules") or [])])
    return rows


def run_one(acc):
    fnagz = os.path.join(FNA, acc + ".fna.gz")
    gffgz = os.path.join(GFF, acc + ".gff.gz")
    dest = os.path.join(DEST, acc + ".regions.tsv")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return acc, "cached", ""
    if not (os.path.exists(fnagz) and os.path.exists(gffgz)):
        return acc, "missing-input", ""
    tmpd = tempfile.mkdtemp(prefix="as_" + acc + "_")
    try:
        fna = os.path.join(tmpd, acc + ".fna")
        gff = os.path.join(tmpd, acc + ".gff")
        with gzip.open(fnagz, "rb") as fi, open(fna, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        dropped = write_filtered_gff(gffgz, gff)
        work = os.path.join(tmpd, "work")
        cmd = [ANTISMASH, "--minimal", "--genefinding-tool", "none",
               "--genefinding-gff3", gff, "--cpus", str(PER_JOB_CPUS),
               "--output-dir", work, "--output-basename", acc,
               "--logfile", os.path.join(tmpd, "as.log"), fna]
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=TIMEOUT, env=RUNENV)
        res = os.path.join(work, acc + ".json")
        if p.returncode != 0 or not os.path.exists(res):
            err = (p.stderr or "")[-700:]
            try:
                err += "\n--log--\n" + open(os.path.join(tmpd, "as.log")).read()[-700:]
            except OSError:
                pass
            return acc, "fail", err
        # --minimal 不產生 regions.js（那是 HTML 輸出的一部分），結果在 <acc>.json。
        # 那個 json 含每條 replicon 的全長序列，554 份約 1.4 GB，所以在這裡就解析
        # 成精簡表，原始 json 不留。
        rows = parse_regions(acc, res)
        os.makedirs(DEST, exist_ok=True)
        with open(dest + ".part", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(REGION_COLS)
            w.writerows(rows)
        os.replace(dest + ".part", dest)
        note = f"{len(rows)} regions"
        if dropped:
            note += " | dropped " + ",".join(f"{k}x{v}" for k, v in dropped.items())
        return acc, "ok", note
    except subprocess.TimeoutExpired:
        return acc, "timeout", ""
    except Exception as exc:                                  # noqa: BLE001
        return acc, "error", repr(exc)[:300]
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def main():
    os.makedirs(DEST, exist_ok=True)
    accs = panel()
    log = open(os.path.join(OUT, "antismash.log"), "a", buffering=1)
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
    print(f"done {done}, {bad} bad", file=sys.stderr)


if __name__ == "__main__":
    main()
