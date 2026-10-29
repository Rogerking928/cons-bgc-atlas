#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""外部集：抓 draft（非完整）基因體與非目標物種，供獨立驗證用。

正文的一切都建立在 554 株完整基因體上，引子也是在同一批上設計的。這支抓兩批
完全沒有進過設計集的基因體：

  1. 同樣 12 種 CoNS 的 scaffold／contig 等級組裝，每種最多 CAP 株；
  2. 非目標物種（S. aureus、Cutibacterium、Enterococcus、E. coli），用來測特異度。

下載 protein.faa ＋ genomic fna。輸出在 data/ext/。
"""
import io, json, os, time, urllib.parse, urllib.request, zipfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data", "ext")
FAA = os.path.join(DATA, "faa")
FNA = os.path.join(DATA, "fna")
for d in (FAA, FNA):
    os.makedirs(d, exist_ok=True)
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{}/dataset_report"
DL = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/download"
CAP = 40
PANEL = ["Staphylococcus epidermidis", "Staphylococcus haemolyticus",
         "Staphylococcus capitis", "Staphylococcus lugdunensis",
         "Staphylococcus hominis", "Staphylococcus saprophyticus",
         "Staphylococcus simulans", "Staphylococcus warneri",
         "Staphylococcus pasteuri", "Staphylococcus cohnii",
         "Staphylococcus caprae", "Staphylococcus pettenkoferi"]
NONTARGET = {"Staphylococcus aureus": 30, "Cutibacterium acnes": 12,
             "Enterococcus faecalis": 10, "Escherichia coli": 8,
             "Streptococcus mitis": 6, "Micrococcus luteus": 6}


def reports(taxon, level, page_size=200):
    q = ["filters.assembly_source=refseq", "filters.exclude_atypical=true",
         f"filters.assembly_level={level}", f"page_size={page_size}"]
    url = API.format(urllib.parse.quote(taxon)) + "?" + "&".join(q)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as fh:
                return json.load(fh).get("reports", [])
        except Exception:
            if attempt == 3:
                return []
            time.sleep(5)


def fetch(accs):
    todo = [a for a in accs
            if not (os.path.exists(os.path.join(FAA, a + ".faa"))
                    and os.path.exists(os.path.join(FNA, a + ".fna")))]
    for i in range(0, len(todo), 10):
        chunk = todo[i:i + 10]
        body = json.dumps({"accessions": chunk,
                           "include_annotation_type": ["PROT_FASTA", "GENOME_FASTA"]}).encode()
        req = urllib.request.Request(DL, data=body,
                                     headers={"Content-Type": "application/json"})
        data = None
        for attempt in range(3):
            try:
                data = urllib.request.urlopen(req, timeout=1200).read()
                break
            except Exception as e:
                print("  retry", e, flush=True)
                time.sleep(8)
        if data is None:
            continue
        try:
            z = zipfile.ZipFile(io.BytesIO(data))
        except Exception:
            continue
        for n in z.namelist():
            p = n.split("/")
            if len(p) < 3:
                continue
            a = p[2]
            if n.endswith("protein.faa"):
                open(os.path.join(FAA, a + ".faa"), "wb").write(z.read(n))
            elif n.endswith(".fna") and "genomic" in n and "cds" not in n and "rna" not in n:
                open(os.path.join(FNA, a + ".fna"), "wb").write(z.read(n))
        print(f"  {i + len(chunk)}/{len(todo)}", flush=True)


def main():
    meta = {}
    accs = []
    for t in PANEL:
        rs = []
        for lvl in ("scaffold", "contig"):
            rs += reports(t, lvl)
            if len(rs) >= CAP:
                break
        rs = rs[:CAP]
        for r in rs:
            meta[r["accession"]] = {"species": t.split()[1], "group": "draft CoNS",
                                    "level": r["assembly_info"]["assembly_level"],
                                    "bioproject": r["assembly_info"].get("bioproject_accession", "")}
        accs += [r["accession"] for r in rs]
        print(f"{t:34s} {len(rs)} drafts", flush=True)
        time.sleep(0.5)
    for t, cap in NONTARGET.items():
        rs = reports(t, "complete_genome")[:cap]
        for r in rs:
            meta[r["accession"]] = {"species": t, "group": "non-target",
                                    "level": "complete",
                                    "bioproject": r["assembly_info"].get("bioproject_accession", "")}
        accs += [r["accession"] for r in rs]
        print(f"{t:34s} {len(rs)} non-target", flush=True)
        time.sleep(0.5)
    json.dump(meta, open(os.path.join(DATA, "meta.json"), "w"), indent=1)
    print(f"total {len(accs)} assemblies", flush=True)
    fetch(accs)
    print("faa", len(os.listdir(FAA)), "fna", len(os.listdir(FNA)), flush=True)


if __name__ == "__main__":
    main()
