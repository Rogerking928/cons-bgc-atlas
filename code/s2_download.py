#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Download RefSeq protein FASTA + genomic GFF for every assembly in the reports.

Proteins go to data/faa/<acc>.faa ; annotation to data/gff/<acc>.gff.gz
"""
import glob, gzip, io, json, os, time, zipfile, urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(HERE, "data", "faa")
GFF = os.path.join(HERE, "data", "gff")
os.makedirs(FAA, exist_ok=True); os.makedirs(GFF, exist_ok=True)
URL = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/download"

accs = []
for p in sorted(glob.glob(os.path.join(HERE, "data", "reports", "*.json"))):
    accs += [r["accession"] for r in json.load(open(p))]
accs = sorted(set(accs))
todo = [a for a in accs
        if not (os.path.exists(os.path.join(FAA, a + ".faa"))
                and os.path.exists(os.path.join(GFF, a + ".gff.gz")))]
print(f"total {len(accs)}, todo {len(todo)}", flush=True)

CHUNK = 20
for i in range(0, len(todo), CHUNK):
    chunk = todo[i:i + CHUNK]
    body = json.dumps({"accessions": chunk,
                       "include_annotation_type": ["PROT_FASTA", "GENOME_GFF"]}).encode()
    req = urllib.request.Request(URL, data=body,
                                 headers={"Content-Type": "application/json"})
    data = None
    for attempt in range(4):
        try:
            data = urllib.request.urlopen(req, timeout=1800).read()
            break
        except Exception as e:
            print("retry", i, e, flush=True)
            time.sleep(8 * (attempt + 1))
    if data is None:
        print("chunk FAIL", i, flush=True); continue
    try:
        z = zipfile.ZipFile(io.BytesIO(data))
    except Exception as e:
        print("zip FAIL", i, e, flush=True); continue
    for n in z.namelist():
        parts = n.split("/")
        if len(parts) < 3:
            continue
        a = parts[2]
        if n.endswith("protein.faa"):
            open(os.path.join(FAA, a + ".faa"), "wb").write(z.read(n))
        elif n.endswith("genomic.gff"):
            with gzip.open(os.path.join(GFF, a + ".gff.gz"), "wb") as fh:
                fh.write(z.read(n))
    print(f"chunk {i}: faa {len(os.listdir(FAA))} gff {len(os.listdir(GFF))}", flush=True)
    time.sleep(0.5)
print("done", len(os.listdir(FAA)), len(os.listdir(GFF)), flush=True)
