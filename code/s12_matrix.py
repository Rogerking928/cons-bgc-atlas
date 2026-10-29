#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Per-genome occupancy matrix: which cluster families each assembly carries,
where they sit, and which lone structural/immunity genes are left over.
"""
import collections, csv, gzip, json, os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
STRUCT = {"Bacteriocin_II", "Bacteriocin_IIc", "Bacteriocin_IId",
          "Bacteriocin_IIi", "BacteriocIIc_cy", "Antimicrobial17", "MccV",
          "L_biotic_typeA", "Mersacidin", "EntA_Immun"}


def main():
    meta = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "metadata.csv")))}
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    fams = list(csv.DictReader(open(os.path.join(OUT, "families.tsv")), delimiter="\t"))
    keep = [f["family"] for f in fams if int(f["n_assemblies"]) >= 4]

    carry = collections.defaultdict(dict)      # assembly -> family -> location
    inside = collections.defaultdict(set)      # assembly -> proteins in a cluster
    for c in cl:
        f = fam_of[c["cluster"]]
        loc = "plasmid" if c["replicon_type"] == "plasmid" else "chromosome"
        prev = carry[c["assembly"]].get(f)
        carry[c["assembly"]][f] = loc if prev in (None, loc) else "both"
    for r in csv.DictReader(open(os.path.join(OUT, "cluster_cds.tsv")), delimiter="\t"):
        inside[r["assembly"]].add(r["protein"])

    orphan = collections.defaultdict(set)
    with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            a, prot, pf, lab = line.split("\t")[:4]
            if lab in STRUCT and prot not in inside.get(a, ()):
                orphan[a].add(lab)

    rows = []
    for a, m in sorted(meta.items()):
        row = {"assembly": a, "species": m["species"].split()[1],
               "strain": m["strain"], "bioproject": m["bioproject"],
               "isolation_source": m["isolation_source"],
               "host_disease": m["host_disease"],
               "collection_date": m["collection_date"][:4],
               "geo": m["geo_loc_name"].split(":")[0]}
        for f in keep:
            row[f] = carry.get(a, {}).get(f, "")
        row["n_families"] = sum(1 for f in keep if row[f])
        row["orphan_struct"] = ",".join(sorted(orphan.get(a, ())))
        rows.append(row)
    p = os.path.join(OUT, "genome_matrix.csv")
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(len(rows), "genomes,", len(keep), "families ->", p)
    for f in keep:
        n = sum(1 for r in rows if r[f])
        pl = sum(1 for r in rows if r[f] in ("plasmid", "both"))
        print(f"  {f}: {n:4d} genomes ({100*n/len(rows):5.1f}%), plasmid-borne in {pl}")
    print("  genomes with no family:", sum(1 for r in rows if r["n_families"] == 0))
    print("  of those, with orphan structural genes:",
          sum(1 for r in rows if r["n_families"] == 0 and r["orphan_struct"]))


if __name__ == "__main__":
    main()
