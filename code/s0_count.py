#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""候選 CoNS 物種的完整基因體數與 isolation_source 覆蓋率（決定收案名單用）。"""
import json, time, urllib.parse, urllib.request

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{}/dataset_report"
CAND = ["Staphylococcus epidermidis", "Staphylococcus haemolyticus",
        "Staphylococcus capitis", "Staphylococcus lugdunensis",
        "Staphylococcus hominis", "Staphylococcus saprophyticus",
        "Staphylococcus warneri", "Staphylococcus pettenkoferi",
        "Staphylococcus caprae", "Staphylococcus simulans",
        "Staphylococcus cohnii", "Staphylococcus pasteuri",
        "Staphylococcus sciuri", "Staphylococcus petrasii",
        "Staphylococcus succinus", "Staphylococcus xylosus"]


def fetch(taxon, level="complete_genome"):
    out, token = [], None
    while True:
        q = ["filters.assembly_source=refseq", "filters.exclude_atypical=true",
             f"filters.assembly_level={level}", "page_size=500"]
        if token:
            q.append("page_token=" + token)
        url = API.format(urllib.parse.quote(taxon)) + "?" + "&".join(q)
        with urllib.request.urlopen(url, timeout=180) as fh:
            d = json.load(fh)
        out += d.get("reports", [])
        token = d.get("next_page_token")
        if not token:
            return out
        time.sleep(1)


if __name__ == "__main__":
    for t in CAND:
        r = fetch(t)
        src = sum(1 for x in r
                  if any(a["name"] == "isolation_source"
                         for a in x["assembly_info"].get("biosample", {}).get("attributes", [])))
        print(f"{t:34s} n={len(r):4d}  isolation_source={src}", flush=True)
        time.sleep(0.5)
