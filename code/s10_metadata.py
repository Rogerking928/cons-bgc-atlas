#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Flatten the assembly reports into one metadata table."""
import csv, glob, json, os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FIELDS = ["assembly", "species", "strain", "bioproject", "biosample",
          "submitter", "release_date", "collection_date", "geo_loc_name",
          "host", "host_disease", "isolation_source", "sequencing_tech",
          "genome_size", "gc_percent", "n_replicons", "total_genes"]

rows = []
for p in sorted(glob.glob(os.path.join(HERE, "data", "reports", "*.json"))):
    for r in json.load(open(p)):
        ai = r.get("assembly_info", {})
        bs = ai.get("biosample", {})
        at = {a["name"]: a.get("value", "") for a in bs.get("attributes", [])}
        st = r.get("assembly_stats", {})
        an = r.get("annotation_info", {}).get("stats", {}).get("gene_counts", {})
        rows.append({
            "assembly": r["accession"],
            "species": r["organism"]["organism_name"],
            "strain": r["organism"].get("infraspecific_names", {}).get("strain", ""),
            "bioproject": ai.get("bioproject_accession", ""),
            "biosample": bs.get("accession", ""),
            "submitter": ai.get("submitter", ""),
            "release_date": ai.get("release_date", ""),
            "collection_date": bs.get("collection_date", at.get("collection_date", "")),
            "geo_loc_name": bs.get("geo_loc_name", at.get("geo_loc_name", "")),
            "host": bs.get("host", at.get("host", "")),
            "host_disease": bs.get("host_disease", at.get("host_disease", "")),
            "isolation_source": bs.get("isolation_source", at.get("isolation_source", "")),
            "sequencing_tech": ai.get("sequencing_tech", ""),
            "genome_size": st.get("total_sequence_length", ""),
            "gc_percent": st.get("gc_percent", ""),
            "n_replicons": st.get("number_of_contigs", ""),
            "total_genes": an.get("total", ""),
        })
with open(os.path.join(OUT, "metadata.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    w.writeheader(); w.writerows(rows)
print(len(rows), "-> out/metadata.csv")
