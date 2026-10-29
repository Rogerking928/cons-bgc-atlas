#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RefSeq 完整基因體的 assembly report（含 BioSample 屬性），每個物種一個 JSON。

收案名單＝人身上會碰到的凝固酶陰性葡萄球菌 12 種。
排除 sciuri（已改列 Mammaliicoccus，動物來源）、xylosus／succinus（食品／動物）。
"""
import json, os, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "data", "reports")
os.makedirs(OUT, exist_ok=True)
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{}/dataset_report"

TAXA = ["Staphylococcus epidermidis", "Staphylococcus haemolyticus",
        "Staphylococcus capitis", "Staphylococcus lugdunensis",
        "Staphylococcus hominis", "Staphylococcus saprophyticus",
        "Staphylococcus simulans", "Staphylococcus warneri",
        "Staphylococcus pasteuri", "Staphylococcus cohnii",
        "Staphylococcus caprae", "Staphylococcus pettenkoferi"]
LEVEL = "complete_genome"


def fetch(taxon, level=LEVEL, page_size=500):
    reports, token = [], None
    while True:
        q = ["filters.assembly_source=refseq", "filters.exclude_atypical=true",
             f"filters.assembly_level={level}", f"page_size={page_size}"]
        if token:
            q.append("page_token=" + token)
        url = API.format(urllib.parse.quote(taxon)) + "?" + "&".join(q)
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=180) as fh:
                    d = json.load(fh)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(5 * (attempt + 1))
        reports += d.get("reports", [])
        token = d.get("next_page_token")
        if not token:
            return reports
        time.sleep(1.0)


if __name__ == "__main__":
    tot = 0
    for t in TAXA:
        rs = fetch(t)
        p = os.path.join(OUT, t.replace(" ", "_") + f".{LEVEL}.json")
        json.dump(rs, open(p, "w"))
        tot += len(rs)
        print(f"{t:34s} {len(rs):4d} -> {os.path.basename(p)}", flush=True)
        time.sleep(1)
    print("total", tot, flush=True)
