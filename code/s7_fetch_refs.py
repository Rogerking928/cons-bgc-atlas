#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""命名基因簇用的參考肽／酵素（UniProt）。

不硬寫 accession：用名稱＋分類（Staphylococcus）查 UniProt REST，把查到的
accession、蛋白名、物種、長度一起記進 data/refs.tsv，事後可核對。
查不到就留白並印出來，**不要用推測的 accession 填**。
"""
import json, os, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DST = os.path.join(HERE, "data", "refs.faa")
TSV = os.path.join(HERE, "data", "refs.tsv")
API = "https://rest.uniprot.org/uniprotkb/search"

QUERIES = [
    ("gallidermin", 'protein_name:gallidermin AND taxonomy_name:Staphylococcus'),
    ("epidermin_EpiA", 'gene:epiA AND taxonomy_name:Staphylococcus'),
    ("Pep5", 'protein_name:"Pep5" AND taxonomy_name:Staphylococcus'),
    ("epilancin_K7", 'protein_name:epilancin AND taxonomy_name:Staphylococcus'),
    ("epicidin_280", 'protein_name:epicidin AND taxonomy_name:Staphylococcus'),
    ("nukacin", 'protein_name:nukacin AND taxonomy_name:Staphylococcus'),
    ("staphylococcin_C55", '"staphylococcin C55" AND taxonomy_name:Staphylococcus'),
    ("aureocin_A53", 'protein_name:aureocin AND taxonomy_name:Staphylococcus'),
    ("lugdunin_LugA", 'protein_name:lugdunin AND taxonomy_name:Staphylococcus'),
    ("lugdunin_LugC", '"lugdunin" AND taxonomy_name:Staphylococcus AND length:[1500 TO 4000]'),
    ("lugdunin_LugD", 'gene:lugD AND taxonomy_name:Staphylococcus'),
    ("aureusimine_AusA", '"aureusimine" AND taxonomy_name:Staphylococcus AND length:[2000 TO 3000]'),
    ("epifadin", '"epifadin" AND taxonomy_name:Staphylococcus'),
    ("GrsA_reference", 'accession:P0C061'),
]


def search(q, size=3):
    url = (API + "?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "size": size,
         "fields": "accession,protein_name,organism_name,length,reviewed"}))
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return json.load(fh).get("results", [])
        except Exception:
            if attempt == 3:
                return []
            time.sleep(4 * (attempt + 1))


def seq_of(acc):
    u = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    s = urllib.request.urlopen(u, timeout=120).read().decode()
    return "".join(s.split("\n")[1:])


if __name__ == "__main__":
    rows, seen = [], set()
    with open(DST, "w") as out:
        for label, q in QUERIES:
            hits = search(q)
            if not hits:
                print("NOT FOUND:", label, flush=True)
                rows.append((label, "", "", "", ""))
                continue
            # 優先取 reviewed（SwissProt）
            hits.sort(key=lambda h: (not h.get("entryType", "").startswith("UniProtKB reviewed"),))
            h = hits[0]
            acc = h["primaryAccession"]
            if acc in seen:
                continue
            seen.add(acc)
            pname = (h.get("proteinDescription", {}).get("recommendedName", {})
                     .get("fullName", {}).get("value")
                     or h.get("proteinDescription", {}).get("submissionNames", [{}])[0]
                     .get("fullName", {}).get("value", ""))
            org = h.get("organism", {}).get("scientificName", "")
            s = seq_of(acc)
            out.write(f">{acc} {label} | {pname} | {org}\n{s}\n")
            rows.append((label, acc, pname, org, str(len(s))))
            print(f"{label:22s} {acc:10s} {len(s):5d}  {pname[:44]:46s} {org}", flush=True)
            time.sleep(0.4)
    with open(TSV, "w") as fh:
        fh.write("label\taccession\tprotein_name\torganism\tlength\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("->", DST)
