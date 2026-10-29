#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文表格的表身（TSV，之後由 make_doc.py 排進 Word）。"""
import collections, csv, json, os, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
CALL = 0.5


def w(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(header)
        wr.writerows(rows)
    print(f"  {name}: {len(rows)} rows", flush=True)


def main():
    md = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "metadata.csv")))}
    mat = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv")))}
    fams = list(csv.DictReader(open(os.path.join(OUT, "families.tsv")), delimiter="\t"))
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
    loc = list(csv.DictReader(open(os.path.join(OUT, "colocalisation.csv"))))
    amr = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv")))}
    mlst = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "mlst.csv")))} \
        if os.path.exists(os.path.join(OUT, "mlst.csv")) else {}
    keep = [f["family"] for f in fams if int(f["n_assemblies"]) >= 4]

    # Table 1 面板組成
    rows = []
    by_sp = collections.defaultdict(list)
    for a, r in mat.items():
        by_sp[r["species"]].append(a)
    for sp, aa in sorted(by_sp.items(), key=lambda kv: -len(kv[1])):
        src = collections.Counter(names.source_group(md[a]["isolation_source"]) for a in aa)
        sizes = [int(md[a]["genome_size"]) / 1e6 for a in aa]
        nrep = [int(md[a]["n_replicons"]) for a in aa]
        rows.append([f"S. {sp}", len(aa),
                     len({md[a]["bioproject"] for a in aa}),
                     len(aa) - src["not stated"],
                     src["blood"], src["device / implant"], src["skin / nose"],
                     f"{statistics.median(sizes):.2f}",
                     f"{statistics.median(nrep):.0f}",
                     sum(1 for a in aa if float(amr[a].get("mecA", 0) or 0) >= CALL)])
    rows.append(["all", len(mat), len({r["bioproject"] for r in md.values()}),
                 sum(1 for a in mat if names.source_group(md[a]["isolation_source"]) != "not stated"),
                 sum(r[4] for r in rows), sum(r[5] for r in rows), sum(r[6] for r in rows),
                 "", "", sum(1 for a in mat if float(amr[a].get("mecA", 0) or 0) >= CALL)])
    w("table1_panel.tsv",
      ["Species", "Genomes", "BioProjects", "With isolation source", "Blood",
       "Device/implant", "Skin/nose", "Median genome (Mb)", "Median replicons",
       "mecA+"], rows)

    # Table 2（正文，7 欄）＋ Table S1（被砍掉的那幾欄）
    rows, extra = [], []
    for f in fams:
        fam = f["family"]
        if fam not in keep:
            continue
        cs = [c for c in cl if fam_of[c["cluster"]] == fam]
        sp = collections.Counter(mat[a_]["species"] for a_ in {c["assembly"] for c in cs})
        lr = [r for r in loc if r["family"] == fam]
        shared = sum(1 for r in lr if int(r["n_acquired_on_replicon"]) > 0)
        near = [int(r["nearest_bp"]) for r in lr if r["nearest_bp"] != ""]
        rows.append([fam, names.fam_name(fam), f["kind"], len(cs), f["n_assemblies"],
                     "; ".join(f"{k} {v}" for k, v in sp.most_common(2)),
                     f"{int(f['median_span']) / 1000:.1f}"])
        extra.append([fam, names.fam_short(fam), f["on_plasmid"],
                      sum(1 for c in cs if int(c["mobile_cds"]) > 0),
                      f"{100 * shared / len(cs):.0f}%",
                      f"{statistics.median(near) / 1000:.0f}" if near else "\u2014",
                      f["median_cds"]])
    w("table2_families.tsv",
      ["Family", "Name", "Class", "Clusters", "Genomes", "Species (top 2, genomes)",
       "Median core span (kb)"], rows)
    w("tableS1_family_extra.tsv",
      ["Family", "Cluster", "On plasmid", "With mobility gene",
       "Sharing replicon with AMR gene",
       "Median distance to nearest AMR gene (kb)", "Median CDS in window"], extra)

    # Table 3 PCR
    if os.path.exists(os.path.join(OUT, "primers.csv")):
        best = {}
        for r in csv.DictReader(open(os.path.join(OUT, "primers.csv"))):
            k = r["target"]
            sc = (float(r["sensitivity"]), float(r["specificity"]), -int(r["product_bp"]))
            if k not in best or sc > best[k][0]:
                best[k] = (sc, r)
        rows = [[k, r["family"], names.fam_short(r["family"]), r["marker"],
                 r["forward"], r["reverse"], r["product_bp"],
                 f"{float(r['tm_f']):.1f}/{float(r['tm_r']):.1f}",
                 f"{100 * float(r['sensitivity']):.0f}%",
                 f"{100 * float(r['specificity']):.1f}%",
                 f"{r['amplified_carriers']}/{r['carriers_tested']}",
                 f"{r['amplified_non_carriers']}/{r['non_carriers_tested']}"]
                for k, (_, r) in sorted(best.items())]
        w("tableS4_primers.tsv",
          ["Target", "Family", "Cluster", "Marker domain", "Forward (5'-3')",
           "Reverse (5'-3')", "Product (bp)", "Tm F/R (C)", "Sensitivity",
           "Specificity", "Carriers amplified", "Non-carriers amplified"], rows)

    # Table 4 前驅肽與理論質量
    p = os.path.join(OUT, "precursor_masses.csv")
    if os.path.exists(p):
        rows, seen = [], set()
        for r in csv.DictReader(open(p)):
            key = (r["family"], r["core_seq"])
            if key in seen:
                continue
            seen.add(key)
            n = sum(1 for x in csv.DictReader(open(p))
                    if x["family"] == r["family"] and x["core_seq"] == r["core_seq"])
            # 完整梯子太長，排進 Word 會爆頁；表裡只放範圍，逐格數字留在
            # out/precursor_masses.csv 與 Source_Data.xlsx
            # 已發表的質量是單一同位素，比對一律同型；平均質量另給一欄，
            # 因為 MALDI 在這個質量區間兩種都用得到（2026-09-20 審閱意見）
            lad = [float(v.split(":")[1]) for v in r["ladder_mono"].split(";")]
            lad2 = ([float(v.split(":")[1]) for v in r["ladder_mono_decarb"].split(";")]
                    if r["ladder_mono_decarb"] else [])
            rows.append([r["family"], names.fam_short(r["family"]), n,
                         r["core_aa"], r["core_seq"],
                         r["mono_core"], r["avg_core"],
                         f"{min(lad):.2f}\u2013{max(lad):.2f}",
                         (f"{min(lad2):.2f}\u2013{max(lad2):.2f}" if lad2 else "\u2014"),
                         (f"{r['published']} {r['published_mass']} ({r['published_type']})"
                          if r["published"] else "\u2014"),
                         r["best_ladder_gap"] if r["best_ladder_gap"] != "" else "\u2014"])
        rows.sort(key=lambda x: (x[0], -int(x[2])))
        w("tableS5_masses.tsv",
          ["Family", "Cluster", "Clusters with this core", "Core (aa)",
           "Predicted core peptide",
           "Core monoisotopic (Da)", "Core average (Da)",
           "Monoisotopic dehydration ladder (Da)",
           "Ladder with C-terminal decarboxylation (Da)",
           "Published peptide and mass", "Gap, same mass type (Da)"], rows)

    # Table 5 抗藥共帶（株系分層）
    p = os.path.join(OUT, "amr_stratified.csv")
    if os.path.exists(p):
        # BH 校正後的 q 值在 associations.csv 裡（每個物種一個檢定家族）
        qof = {(r["test"], r["x"], r["y"]): r["q"] for r in csv.DictReader(
            open(os.path.join(OUT, "associations.csv")))}

        def qcell(r):
            q = qof.get((f"family~AMR({r['species']})", r["family"], r["determinant"]))
            return f"{float(q):.2f}" if q not in (None, "") else "\u2014"

        def cmh_cell(r, cut, what):
            v = r.get(f"{what}@{cut}")
            if not v:
                return "\u2014"
            if what == "cmh_odds" and (float(v) == 0.0 or float(v) == float("inf")):
                return "n.e."      # 退化估計：只有一層有資訊時 OR 會塌到 0 或 inf
            return f"{float(v):.2f}" if what == "cmh_odds" else f"{float(v):.1e}"

        rows = []
        for r in csv.DictReader(open(p)):
            if not r["p"] or float(r["p"]) > 0.05:
                continue
            # 有零格時 OR 是 0 或 inf，不可估計：報 Haldane-Anscombe 的點估計並標出來
            orr = (f"{float(r['odds']):.2f}" if not int(r["not_estimable"])
                   else f"n.e. ({float(r['odds_haldane_anscombe']):.1f})")
            flag = "single-project" if int(r["single_project_signal"]) else ""
            rows.append([f"S. {r['species']}", r["family"], names.fam_short(r["family"]),
                         r["determinant"],
                         f"{r['n_fam_det']}/{int(r['n_fam_det']) + int(r['n_fam_nodet'])}",
                         f"{r['n_nofam_det']}/{int(r['n_nofam_det']) + int(r['n_nofam_nodet'])}",
                         orr, f"{float(r['p']):.1e}",
                         qcell(r),
                         cmh_cell(r, "0.01", "cmh_odds"), cmh_cell(r, "0.01", "cmh_p"),
                         r.get("informative_strata@0.01", ""),
                         cmh_cell(r, "0.02", "cmh_odds"), cmh_cell(r, "0.02", "cmh_p"),
                         r.get("informative_strata@0.02", ""),
                         f"{float(r['worst_p_leave_one_project_out']):.1e}",
                         r["worst_project"], flag])
        rows.sort(key=lambda x: float(x[7]))
        w("tableS2_amr.tsv",
          ["Species", "Family", "Cluster", "Determinant",
           "Carriers with determinant", "Non-carriers with determinant",
           "Crude OR", "Crude p", "BH q", "CMH OR (cut 0.01)", "CMH p (cut 0.01)",
           "Informative strata (0.01)", "CMH OR (cut 0.02)", "CMH p (cut 0.02)",
           "Informative strata (0.02)", "Worst leave-one-project-out p",
           "Project removed", "Flag"], rows)

    # Table 6 定序型 × 家族
    if mlst:
        rows = []
        for sp in ("epidermidis", "lugdunensis", "haemolyticus", "hominis"):
            aa = [a for a, r in mat.items() if r["species"] == sp]
            sts = collections.Counter(mlst.get(a, {}).get("ST", "") for a in aa)
            for st, n in sts.most_common():
                if not st or n < 4:
                    continue
                g = [a for a in aa if mlst.get(a, {}).get("ST") == st]
                car = {f: sum(1 for a in g if mat[a].get(f)) for f in keep}
                car = {k: v for k, v in car.items() if v}
                rows.append([f"S. {sp}", f"ST{st}", n,
                             "; ".join(f"{k} {v}" for k, v in sorted(car.items())) or "none",
                             sum(1 for a in g if float(amr[a].get("mecA", 0) or 0) >= CALL),
                             sum(1 for a in g if float(amr[a].get("erm(A)", 0) or 0) >= CALL),
                             sum(1 for a in g if float(amr[a].get("blaZ", 0) or 0) >= CALL)])
        w("tableS3_mlst.tsv",
          ["Species", "Sequence type", "Genomes", "Cluster families carried",
           "mecA+", "erm(A)+", "blaZ+"], rows)


if __name__ == "__main__":
    main()
