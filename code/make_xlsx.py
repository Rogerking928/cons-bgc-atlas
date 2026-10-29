#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Source_Data.xlsx：**每一格圖一個分頁**，裡面就是圖上畫出來的那些數字。

之前的版本是「每個輸出檔一個分頁」，讀者要自己從檔案裡湊出某一格畫的是什麼。
這版改成一格一頁（F1A、F1B…、S9B），欄名和圖上的座標軸／圖例用同一套字，
另外每張表一頁。第一頁是索引，寫明哪個分頁對應哪一格、以及它是從 out/ 的哪個檔算出來的。

  python3 code/make_xlsx.py   ->  out/Source_Data.xlsx
"""
import collections, csv, json, os, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

import names

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
WB = Workbook()
WB.remove(WB.active)
INDEX = []


def read(name, delim=","):
    return list(csv.DictReader(open(os.path.join(OUT, name)), delimiter=delim))


def js(name):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else {}


def sheet(tag, what, source, header, rows):
    """一個分頁＝一格圖（或一張表）。第一列寫它是什麼、從哪來。"""
    ws = WB.create_sheet(tag[:31])
    ws.cell(1, 1, f"{tag} — {what}").font = Font(bold=True)
    ws.cell(2, 1, f"computed from out/{source}").font = Font(italic=True, size=9)
    for j, h in enumerate(header, start=1):
        c = ws.cell(4, j, h)
        c.font = Font(bold=True)
        c.alignment = Alignment(wrap_text=True, vertical="bottom")
    for i, r in enumerate(rows, start=5):
        for j, v in enumerate(r, start=1):
            ws.cell(i, j, v)
    for j in range(1, len(header) + 1):
        # 有些分頁的列長度不一（例如 F5B 底下接了第二段小表），取值時要防呆
        widths = [len(str(header[j - 1]))]
        widths += [len(str(r[j - 1])) for r in rows[:50] if len(r) >= j]
        ws.column_dimensions[get_column_letter(j)].width = max(12, min(40, max(widths) + 2))
    INDEX.append([tag, what, f"out/{source}", len(rows)])
    return ws


# ------------------------------------------------------------------ inputs
md = {r["assembly"]: r for r in read("metadata.csv")}
mat = list(read("genome_matrix.csv"))
matd = {r["assembly"]: r for r in mat}
cl = read("clusters.tsv", "\t")
cld = {r["cluster"]: r for r in cl}
ccds = read("cluster_cds.tsv", "\t")
fams = read("families.tsv", "\t")
fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
KEEP = [f["family"] for f in fams if int(f["n_assemblies"]) >= 4]
KIND = {f["family"]: f["kind"] for f in fams}
loc = read("colocalisation.csv")
amr = {r["assembly"]: r for r in read("amr_per_genome.csv")}
mlst = {r["assembly"]: r for r in read("mlst.csv")}
gc = read("cluster_gc.csv")
dn = read("distance_null.csv")
dnj = js("distance_null.json")
rule = js("calling_rule.json")
strat = read("amr_stratified.csv")
ext = read("external_pcr.csv")
extj = js("external_pcr.json")
deep = read("deepbgc_results.csv")
prim = read("primers.csv")
mass = read("precursor_masses.csv")
base = js("baselines.json")
paired = js("paired_tests.json")
hard = js("hard_negatives.json")
ladder = js("ladder.json")
rep = js("precursor_report.json")
gain = js("atlas_gain_summary.json")
salj = js("finetune_summary.json")
mutj = js("mutagenesis.json")
mutr = read("mutagenesis.csv")
calj = js("fdr_calibration.json")
adom = read("adomain_codes.csv")
rulebase = js("rule_baseline.json")
ftap = js("finetune_ap.json")
counts = collections.Counter(r["species"] for r in mat)
CARRIERS = {r["species"] for r in mat if any(r.get(f) for f in KEEP)}
ORDER = [s for s in ("epidermidis", "haemolyticus", "capitis", "lugdunensis",
                     "hominis", "saprophyticus", "simulans", "warneri")
         if s in CARRIERS] + sorted(CARRIERS - {"epidermidis", "haemolyticus",
                                                "capitis", "lugdunensis", "hominis",
                                                "saprophyticus", "simulans", "warneri"})
best_primer = {}
for r in prim:
    k = r["target"]
    sc = (float(r["sensitivity"]), float(r["specificity"]), -int(r["product_bp"]))
    if k not in best_primer or sc > best_primer[k][0]:
        best_primer[k] = (sc, r)
PRI = {k: v[1] for k, v in best_primer.items()}

# ------------------------------------------------------------------ Figure 1
kinds = collections.defaultdict(lambda: collections.defaultdict(set))
for c in cl:
    kinds[matd[c["assembly"]]["species"]][c["kind"]].add(c["assembly"])
# Figure 1A 畫的是「≥20 株的物種」，包含完全沒有基因簇的 S. saprophyticus
# （圖上標 none detected），所以不能用只含帶菌物種的 ORDER
allsp = [s for s in ("epidermidis", "haemolyticus", "capitis", "lugdunensis",
                     "hominis", "saprophyticus", "simulans", "warneri")
         if counts[s] >= 20]
sheet("F1A", "genomes of each species carrying a cluster of each class",
      "clusters.tsv + genome_matrix.csv",
      ["species", "genomes in panel", "with an NRPS/PKS cluster", "with a hybrid cluster",
       "with a RiPP cluster", "NRPS/PKS (%)", "hybrid (%)", "RiPP (%)"],
      [[s, counts[s]] + [len(kinds[s][k]) for k in ("NRPS/PKS", "hybrid", "RiPP")]
       + [round(100 * len(kinds[s][k]) / counts[s], 1)
          for k in ("NRPS/PKS", "hybrid", "RiPP")] for s in allsp])

rows = []
for fam in KEEP:
    for sp in ORDER:
        n = sum(1 for r in mat if r["species"] == sp and r.get(fam))
        if n:
            rows.append([fam, names.fam_short(fam), sp, n, counts[sp],
                         round(100 * n / counts[sp], 1)])
sheet("F1B", "family x species heat map: cell number and cell shade",
      "genome_matrix.csv", ["family", "cluster", "species", "genomes carrying",
                            "genomes of that species", "carried by (% of species)"], rows)

rows = []
for fam in KEEP:
    cs = [c for c in cl if fam_of[c["cluster"]] == fam]
    p = sum(1 for c in cs if c["replicon_type"] != "chromosome")
    rows.append([fam, names.fam_short(fam), len(cs) - p, p,
                 "plasmid" if p == len(cs) else "chromosome" if p == 0 else "both"])
sheet("F1B_squares", "replicon square beside each family label", "clusters.tsv",
      ["family", "cluster", "clusters on a chromosome", "clusters on a plasmid",
       "square colour"], rows)

rows = []
for fam in KEEP:
    v = [int(c["span"]) / 1000 for c in cl if fam_of[c["cluster"]] == fam]
    rows.append([fam, names.fam_short(fam), KIND[fam], len(v), round(min(v), 2),
                 round(statistics.median(v), 2), round(max(v), 2)])
sheet("F1C", "core cluster span per family (dot = median, line = range)",
      "clusters.tsv", ["family", "cluster", "class", "clusters", "min span (kb)",
                       "median span (kb)", "max span (kb)"], rows)

rows = []
for f in sorted({fam_of[c["cluster"]] for c in cl}, key=lambda x: int(x[1:])):
    cs = [c for c in cl if fam_of[c["cluster"]] == f]
    kk = collections.Counter(c["kind"] for c in cs)
    rows.append([f, names.fam_short(f) if f in KEEP else "(not drawn in Figure 1)",
                 KIND.get(f, ""), len(cs), kk["NRPS/PKS"], kk["hybrid"], kk["RiPP"],
                 "yes" if f in KEEP else "no"])
sheet("F1_reconcile",
      "how the 562 clusters map onto families: the class of each cluster against the "
      "class of its family, and which families Figure 1B and 1C draw",
      "clusters.tsv + families.tsv",
      ["family", "cluster", "family class (Table 2, F1C)", "clusters",
       "clusters classed NRPS/PKS", "clusters classed hybrid", "clusters classed RiPP",
       "drawn in Figure 1B and 1C"], rows)

# ------------------------------------------------------------------ Figure 2
for tag, cid in (("F2A", "C00115"), ("F2B", "C00077")):
    g = sorted([r for r in ccds if r["cluster"] == cid], key=lambda r: int(r["start"]))
    lo = int(g[0]["start"])
    sheet(tag, f"gene map of cluster {cid} ({fam_of[cid]}, "
               f"{cld[cid]['replicon_type']})", "cluster_cds.tsv",
          ["protein", "locus tag", "start (bp)", "end (bp)", "offset from left edge (bp)",
           "strand", "length (aa)", "domains", "product"],
          [[r["protein"], r["locus_tag"], r["start"], r["end"], int(r["start"]) - lo,
            r["strand"], r["aa"], r["domains"], r["product"]] for r in g])

sheet("F2C", "distance from each cluster to the nearest acquired resistance gene, "
             "observed and under the permutation null",
      "distance_null.csv",
      ["cluster", "family", "replicon type", "replicon length (bp)",
       "acquired resistance genes on the replicon", "observed distance (bp)",
       "null median distance (bp)", "one-sided P"],
      [[r["cluster"], r["family"], r["replicon_type"], r["replicon_len"], r["n_amr"],
        r["observed_bp"], round(float(r["null_median_bp"]), 1), r["p_one_sided"]]
       for r in dn])
def _ci(v, key):
    """區間取兩位有效的一位小數；沒有就留白（每簇那層本來就沒算區間）。"""
    c = v.get(key)
    return ("", "") if not c else (round(c[0], 1), round(c[1], 1))


sheet("F2C_tests", "the two tests drawn on F2C, with the family-level "
      "bootstrap intervals the text reports", "distance_null.json",
      ["level", "replicon type", "n", "observed median (kb)",
       "observed 95% CI lo", "observed 95% CI hi", "null median (kb)",
       "null 95% CI lo", "null 95% CI hi", "median difference (kb)",
       "difference 95% CI lo", "difference 95% CI hi", "closer than null",
       "Wilcoxon P"],
      [["per cluster", k, v["n_clusters"], round(v["observed_median_kb"], 1),
        "", "", round(v["null_median_kb"], 1), "", "", "", "", "",
        "", v["wilcoxon_p"]]
       for k, v in dnj.get("groups", {}).items()] +
      [["per family", k, v["n_families"], round(v["observed_median_kb"], 1),
        *_ci(v, "observed_ci_kb"), round(v["null_median_kb"], 1),
        *_ci(v, "null_ci_kb"),
        round(v["median_difference_kb"], 1) if "median_difference_kb" in v else "",
        *_ci(v, "difference_ci_kb"),
        f"{v['closer_than_null']}/{v['n_families']}", v["wilcoxon_p"]]
       for k, v in dnj.get("family_level", {}).items() if v.get("wilcoxon_p")])

rows = []
for fam in KEEP:
    rs = [r for r in loc if r["family"] == fam]
    if not rs:
        continue
    sh = sum(1 for r in rs if int(r["n_acquired_on_replicon"]) > 0)
    rows.append([fam, names.fam_short(fam), KIND[fam], len(rs), sh,
                 round(100 * sh / len(rs), 1)])
sheet("F2D", "clusters sharing a replicon with an acquired determinant",
      "colocalisation.csv", ["family", "cluster", "class", "clusters",
                             "sharing a replicon", "share (%)"], rows)

cnt = {k: collections.Counter() for k in ("NRPS/PKS", "hybrid", "RiPP")}
for c in cl:
    m = int(c["mobile_cds"])
    cnt[c["kind"]]["0" if m == 0 else "1" if m == 1 else ">=2"] += 1
sheet("F2E", "mobility genes inside the cluster window, by class", "clusters.tsv",
      ["class", "clusters", "0 mobility genes", "1", ">=2", "0 (%)", "1 (%)", ">=2 (%)"],
      [[k, sum(cnt[k].values())] + [cnt[k][x] for x in ("0", "1", ">=2")]
       + [round(100 * cnt[k][x] / max(1, sum(cnt[k].values())), 1)
          for x in ("0", "1", ">=2")] for k in cnt])

sheet("F2F", "cluster GC minus genome GC", "cluster_gc.csv",
      ["cluster", "family", "assembly", "replicon type", "span (bp)", "cluster GC (%)",
       "genome GC (%)", "difference (percentage points)"],
      [[r["cluster"], r["family"], r["assembly"], r["replicon_type"], r["span"],
        r["cluster_gc"], r["genome_gc"], r["delta"]] for r in gc
       if r["family"] in KEEP])

# ------------------------------------------------------------------ Figure 3
rows = []
for defn, d in rule.get("by_ST", {}).items():
    for st, v in d.items():
        rows.append([defn, st, v["n"], v["carriers"],
                     round(100 * v["carriers"] / v["n"], 1)])
sheet("F3A", "S. epidermidis F001 carriage by sequence type, under three definitions",
      "calling_rule.json", ["carriage definition", "sequence type", "genomes",
                            "carriers", "carriage (%)"], rows)

rows = []
for defn, d in rule.get("association", {}).items():
    for det, v in d.items():
        rows.append([defn, det, v["a"], v["b"], v["c"], v["d"], round(v["odds"], 4),
                     round(v["odds_haldane"], 4), round(v["ci"][0], 4),
                     round(v["ci"][1], 4), v["p"]])
sheet("F3B", "F001 against each determinant, under each carriage definition",
      "calling_rule.json",
      ["carriage definition", "determinant", "carriers with determinant",
       "carriers without", "non-carriers with", "non-carriers without", "odds ratio",
       "Haldane-Anscombe odds ratio", "95% CI lower", "95% CI upper", "Fisher P"], rows)

lug = [a for a, r in matd.items() if r["species"] == "lugdunensis"]
lst = collections.Counter(mlst.get(a, {}).get("ST", "") for a in lug)
lkeep = [s for s, n in lst.most_common() if s and n >= 3]
rows = []
for st in lkeep + ["other/untyped"]:
    g = ([a for a in lug if mlst.get(a, {}).get("ST") == st] if st != "other/untyped"
         else [a for a in lug if mlst.get(a, {}).get("ST") not in lkeep])
    rows.append([st, len(g)] + [sum(1 for a in g if matd[a].get(f))
                                for f in ("F003", "F006", "F005")])
sheet("F3C", "S. lugdunensis carriage by sequence type", "genome_matrix.csv + mlst.csv",
      ["sequence type", "genomes", "F003 aryl-NRPS", "F006 lugdunin",
       "F005 NRPS x YcaO hybrid"], rows)

# ------------------------------------------------------------------ Figure 4
targets = sorted(extj.get("per_target", {}))
rows = []
for t in targets:
    fam = t.split("-")[0]
    sp = collections.Counter(x["species"] for x in mat if x.get(fam)).most_common(1)
    sp = sp[0][0] if sp else ""
    tot = sum(1 for x in mat if x["species"] == sp)
    car = sum(1 for x in mat if x["species"] == sp and x.get(fam))
    d = [x for x in ext if x["group"] == "draft CoNS" and x["species"] == sp]
    card = sum(int(x[f"{t}_carrier"]) for x in d)
    p = PRI[t]
    rows.append([t, fam, sp, car, tot, round(100 * car / tot, 1), card, len(d),
                 round(100 * card / max(1, len(d)), 1),
                 p["amplified_carriers"], p["carriers_tested"],
                 p["amplified_non_carriers"], p["non_carriers_tested"]])
sheet("F4A", "inputs to the predictive values: carriage in the panel and in the drafts, "
             "and the in-silico performance the Beta sampling uses",
      "genome_matrix.csv + external_pcr.csv + primers.csv",
      ["target", "family", "main species",
       "panel carriers, main species only", "panel genomes of that species",
       "panel carriage (%)", "draft carriers, main species only",
       "draft genomes of that species", "draft carriage (%)",
       "carriers amplified (panel, all species)",
       "carriers tested (panel, all species)",
       "non-carriers amplified (panel, all species)",
       "non-carriers tested (panel, all species)"], rows)

rows = []
for t in targets:
    e = extj["per_target"][t]
    p = PRI[t]
    rows.append([t, p["amplified_carriers"], p["carriers_tested"],
                 round(100 * float(p["sensitivity"]), 1),
                 e["draft_carriers"],
                 round(100 * e["sensitivity"], 1) if e["sensitivity"] is not None else "",
                 round(100 * e["sensitivity_ci"][0], 1) if e["sensitivity_ci"] else "",
                 round(100 * e["sensitivity_ci"][1], 1) if e["sensitivity_ci"] else ""])
sheet("F4B", "sensitivity on the design panel and on the held-out drafts",
      "primers.csv + external_pcr.json",
      ["target", "carriers amplified (panel, all species)",
       "carriers (panel, all species)", "sensitivity, panel (%)",
       "draft carriers (all species)", "sensitivity, drafts (%)",
       "95% CI lower", "95% CI upper"], rows)

sheet("F4C", "genomes that should not amplify and did",
      "external_pcr.json",
      ["target", "draft non-carriers tested", "draft non-carriers amplified",
       "non-target genomes tested", "non-target genomes amplified"],
      [[t, extj["per_target"][t]["draft_non_carriers"],
        extj["per_target"][t]["draft_non_carriers"]
        - round(extj["per_target"][t]["specificity"]
                * extj["per_target"][t]["draft_non_carriers"]),
        extj["per_target"][t]["non_target_tested"],
        extj["per_target"][t]["non_target_amplified"]] for t in targets])

rows = []
for t in targets:
    fam = t.split("-")[0]
    sp = collections.Counter(x["species"] for x in mat if x.get(fam)).most_common(1)
    sp = sp[0][0] if sp else ""
    tot = sum(1 for x in mat if x["species"] == sp)
    car = sum(1 for x in mat if x["species"] == sp and x.get(fam))
    d = [x for x in ext if x["group"] == "draft CoNS" and x["species"] == sp]
    card = sum(int(x[f"{t}_carrier"]) for x in d)
    rows.append([fam, names.fam_short(fam), sp, car, tot, round(100 * car / max(1, tot), 1),
                 card, len(d), round(100 * card / max(1, len(d)), 1)])
sheet("F4D", "family carriage inside and outside the design panel",
      "genome_matrix.csv + external_pcr.csv",
      ["family", "cluster", "main species", "panel carriers, main species only",
       "panel genomes of that species", "panel carriage (%)",
       "draft carriers, main species only", "draft genomes of that species",
       "draft carriage (%)"], rows)

# ------------------------------------------------------------------ Figure 5
rows = []
for m, cur in base.get("curves", {}).items():
    pr, rc = cur[0], cur[1]
    for p_, r_ in zip(pr, rc):
        rows.append([m, round(r_, 5), round(p_, 5)])
sheet("F5A", "precision-recall curves, grouped cross-validation", "baselines.json",
      ["method", "recall", "precision"], rows)

rows = [[m, round(base["methods"][m]["cv_ap"], 4)] for m in base.get("methods", {})
        if base["methods"][m].get("cv_ap") is not None]
rows.append(["leader rules (this paper's implementation)", round(rulebase.get("ap", 0), 4)])
sheet("F5B", "average precision per method, and the paired differences",
      "baselines.json + rule_baseline.json + paired_tests.json",
      ["method", "average precision"], rows + [[], ["paired difference", "estimate",
                                                    "95% CI lower", "95% CI upper"]] +
      [[k, round(v["difference"], 4), round(v["ci"][0], 4), round(v["ci"][1], 4)]
       for k, v in paired.items() if isinstance(v, dict) and "difference" in v])

sheet("F5C", "negative sets of increasing difficulty", "hard_negatives.json",
      ["negative set", "negatives", "positives", "prevalence (AP of random guessing)",
       "average precision", "AUROC", "average precision, length only"],
      [[k, v["n_negatives"], hard.get("n_positives"), round(v["prevalence"], 4),
        round(v["ap"], 4), round(v["auroc"], 4), round(v["ap_length_only"], 4)]
       for k, v in hard.get("negative_sets", {}).items()])

sheet("F5D", "performance by 10-mer identity to the nearest other training peptide",
      "hard_negatives.json",
      ["identity bin", "training peptides", "average precision", "95% CI lower",
       "95% CI upper"],
      [[k, v["n_positives"], round(v["ap"], 4),
        round(v["ap_ci"][0], 4) if v.get("ap_ci") and v["ap_ci"][0] else "",
        round(v["ap_ci"][1], 4) if v.get("ap_ci") and v["ap_ci"][1] else ""]
       for k, v in hard.get("identity_strata", {}).items()])

sheet("F5E", "learning curve", "hard_negatives.json",
      ["fraction of positive clusters", "training peptides", "average precision"],
      [[d["fraction"], d["n_positives"], round(d["ap"], 4)]
       for d in hard.get("learning_curve", [])])

sheet("F5F", "ESM-2 size ladder", "ladder.json",
      ["model", "parameters (M)", "embedding dimension", "average precision",
       "seconds to embed"],
      [[k, {"8M": 8, "35M": 35, "150M": 150}.get(k), v["dim"], round(v["cv_ap"], 4),
        v["seconds"]] for k, v in ladder.items() if k in ("8M", "35M", "150M")])

# ------------------------------------------------------------------ Figure 6
en = rep.get("enrichment_near_machinery", {})
sheet("F6A", "candidates near cluster machinery, the label-free external check",
      "precursor_report.json",
      ["group", "within 10 kb of machinery", "farther", "share within 10 kb (%)"],
      [["high scorers, no profile hit", en.get("high_near"), en.get("high_far"),
        round(100 * en.get("high_near", 0) / max(1, en.get("high_near", 0) + en.get("high_far", 0)), 2)],
       ["all other short CDS", en.get("rest_near"), en.get("rest_far"),
        round(100 * en.get("rest_near", 0) / max(1, en.get("rest_near", 0) + en.get("rest_far", 0)), 2)],
       ["odds ratio", round(en.get("odds_ratio", 0), 1), "permutation P",
        en.get("permutation_p")]])

sheet("F6B", "clusters with a precursor, before and after the model",
      "atlas_gain_summary.json",
      ["family", "cluster", "class", "clusters", "with a profile precursor",
       "gaining a candidate (q<0.20)", "gaining a candidate (q<0.05)"],
      [[f, names.fam_short(f), KIND.get(f, ""), d["clusters"], d["with_known_precursor"],
        d["gaining_a_candidate"], d["gaining_a_candidate_q05"]]
       for f, d in sorted(gain.items())
       if f.startswith("F") and d["clusters"] >= 4 and KIND.get(f) in ("RiPP", "hybrid")])

sheet("F6C", "top-attributed residue next to a double-glycine motif",
      "finetune_summary.json",
      ["group", "sequences", "top-attributed residue next to a GG motif",
       "share (%)"],
      [[k, v[1], v[0], round(100 * v[0] / v[1], 1)]
       for k, v in salj.get("near_gg", {}).items()])
meta = salj.get("meta_profile", {})
rows = []
for g, prof in meta.items():
    for pos in sorted(prof, key=int):
        rows.append([g, int(pos), round(prof[pos][0], 4), prof[pos][1]])
sheet("F6C_inset", "median integrated-gradient attribution, aligned on the motif",
      "finetune_summary.json",
      ["group", "residue relative to the GG motif", "median attribution",
       "sequences contributing"], rows)

sheet("F6D", "change in detector log-odds after editing the motif",
      "mutagenesis.csv",
      ["group", "edit", "sequence", "log-odds before", "log-odds after", "change"],
      [[r["group"], r["kind"], r["seq"], r["base"], r["mutated"], r["delta"]]
       for r in mutr])
sheet("F6D_summary", "medians drawn on F6D", "mutagenesis.json",
      ["group", "edit", "sequences", "median log-odds before", "median change",
       "median change, control edit", "share whose score falls"],
      [[k, v["kind"], v["n"], round(v["baseline_median_logodds"], 3),
        round(v["delta_median"], 3),
        round(v["control_delta_median"], 3) if v.get("control_delta_median") is not None else "",
        round(v["fraction_dropping"], 3)] for k, v in mutj.items()])

sheet("F6E", "share of sequences above each score threshold", "fdr_calibration.json",
      ["score threshold", "real short CDS (%)", "shuffled background (%)",
       "shuffled precursors (%)"],
      [[round(c["threshold"], 4), round(100 * c["real_rate"], 4),
        round(100 * c["decoy_background_rate"], 4),
        round(100 * c["decoy_from_positives_rate"], 4)] for c in calj.get("curve", [])])
sheet("F6F", "the same at the two calling thresholds", "fdr_calibration.json",
      ["threshold", "score", "real above (%)", "shuffled background above (%)",
       "shuffled precursors above (%)", "real above (n)",
       "shuffled background above (n)", "shuffled precursors above (n)",
       "decoy-based estimator of the false discovery proportion"],
      [[k, round(v["score"], 4), round(100 * v["real_rate"], 3),
        round(100 * v["decoy_background_rate"], 3),
        round(100 * v["decoy_from_positives_rate"], 3), v["real_above"],
        v["decoy_background_above"], v["decoy_from_positives_above"],
        round(v["empirical_fdp_background"], 2)]
       for k, v in calj.get("thresholds", {}).items()])

# ------------------------------------------------------------ supplementary
sheet("S1A", "complete genomes per species", "metadata.csv",
      ["species", "genomes"], [[s, n] for s, n in counts.most_common()])
src = collections.Counter(names.source_group(md[a]["isolation_source"]) for a in matd)
sheet("S1B", "isolation source, grouped by keyword", "metadata.csv",
      ["isolation source group", "genomes"], [[k, v] for k, v in src.most_common()])

rows = []
for fam in ("F001", "F004", "F005", "F011"):
    cs = [c for c, f in fam_of.items() if f == fam]
    if not cs:
        continue
    cid = sorted(cs, key=lambda c: int(cld[c]["n_cds"]))[len(cs) // 2]
    for r in sorted([x for x in ccds if x["cluster"] == cid],
                    key=lambda r: int(r["start"])):
        rows.append([fam, cid, r["protein"], r["start"], r["end"], r["strand"], r["aa"],
                     r["domains"], r["product"]])
sheet("S2", "gene maps of four more representative clusters", "cluster_cds.tsv",
      ["family", "cluster", "protein", "start (bp)", "end (bp)", "strand", "length (aa)",
       "domains", "product"], rows)

rows = []
for sp, fam in (("epidermidis", "F001"), ("lugdunensis", "F006")):
    pool = [a for a, r in matd.items() if r["species"] == sp]
    for g in ["blood", "device / implant", "skin / nose", "other clinical"]:
        sub = [a for a in pool if names.source_group(md[a]["isolation_source"]) == g]
        rows.append([sp, fam, g, len(sub), sum(1 for a in sub if matd[a].get(fam)),
                     round(100 * sum(1 for a in sub if matd[a].get(fam)) / len(sub), 1)
                     if sub else ""])
sheet("S3", "carriage by isolation source (descriptive only)",
      "metadata.csv + genome_matrix.csv",
      ["species", "family", "isolation source group", "genomes", "carriers",
       "carriage (%)"], rows)

dom_rows = []
dom_by = collections.defaultdict(list)
import gzip
with gzip.open(os.path.join(OUT, "domains.tsv.gz"), "rt") as fh:
    next(fh)
    for line in fh:
        a_, prot, pf, lab, af, at, ev, sc, pl = line.rstrip("\n").split("\t")
        dom_by[(a_, prot)].append((int(af), int(at), lab, int(pl)))
for fam, title in (("F006", "lugdunin"), ("F001", "aureusimine (epidermidis)"),
                   ("F003", "aryl-capped NRPS (lugdunensis)"),
                   ("F005", "NRPS x YcaO hybrid")):
    cs = [c for c, f in fam_of.items() if f == fam]
    if not cs:
        continue
    cid = cs[0]
    prots = [r for r in ccds if r["cluster"] == cid and int(r["aa"] or 0) >= 400
             and r["domains"]]
    prots.sort(key=lambda r: -int(r["aa"]))
    for r in prots[:4]:
        for af, at, lab, pl in sorted(dom_by.get((r["assembly"], r["protein"]), [])):
            dom_rows.append([fam, cid, r["assembly"], r["protein"], r["product"],
                             int(r["aa"]), lab, af, at])
sheet("S4A", "domain architecture of the synthetases drawn in Figure S4A",
      "cluster_cds.tsv + domains.tsv.gz",
      ["family", "cluster", "assembly", "protein", "product", "protein length (aa)",
       "domain", "domain start (residue)", "domain end (residue)"], dom_rows)

sheet("S4B", "Stachelhaus specificity codes of the adenylation domains",
      "adomain_codes.csv",
      ["family", "cluster", "assembly", "protein", "product", "domain",
       "signature", "closest published consensus", "matching positions out of 10"],
      [[r["family"], r["cluster"], r["assembly"], r["protein"], r["product"],
        r["domain"], r["signature"], r["closest"], r["identity"]] for r in adom])
prof = js("gc_profiles.json")
pick = next((k for k, v in prof.items() if v["family"] == "F006"), None)
if pick:
    d = prof[pick]
    sheet("S4C", f"GC in 2 kb windows along {d['replicon']}", "gc_profiles.json",
          ["window start (bp)", "GC (%)", "genome mean GC (%)", "cluster start (bp)",
           "cluster end (bp)"],
          [[i * d["window"], v, d["genome_gc"], d["start"], d["end"]]
           for i, v in enumerate(d["gc"])])

tree = js("tree_meta.json")
sheet("S5", "genomes in the similarity tree, with the ring values",
      "tree_meta.json + genome_matrix.csv + amr_per_genome.csv",
      ["assembly", "species", "F001", "F004", "F006", "F009", "mecA"],
      [[a, matd[a]["species"]] + [1 if matd[a].get(f) else 0
                                  for f in ("F001", "F004", "F006", "F009")]
       + [1 if float(amr.get(a, {}).get("mecA", 0) or 0) >= 0.5 else 0]
       for a in tree.get("accs", [])])

scan_seen = {}
for r in read("precursor_scan.csv"):
    scan_seen.setdefault(r["seq"], r)
hold = {h["seq"] for h in read("precursor_holdout.csv")}
sheet("S6A", "detector score distribution, by distinct sequence",
      "precursor_scan.csv + precursor_holdout.csv",
      ["group", "sequences", "median score"],
      [["profile precursor", sum(1 for r in scan_seen.values() if r["known_domain"]),
        round(float(np.median([float(r["score"]) for r in scan_seen.values()
                               if r["known_domain"]])), 4)],
       ["RefSeq-named, no profile", len([s for s in scan_seen if s in hold]),
        round(float(np.median([float(scan_seen[s]["score"]) for s in scan_seen
                               if s in hold])), 4)],
       ["other distinct CDS",
        sum(1 for s, r in scan_seen.items() if not r["known_domain"] and s not in hold),
        round(float(np.median([float(r["score"]) for s, r in scan_seen.items()
                               if not r["known_domain"] and s not in hold])), 4)]])
sheet("S6B", "average precision on the held-out RefSeq-named set, with the two "
             "recalls quoted in the paper (they are different operating points)",
      "baselines.json + precursor_report.json",
      ["method", "held-out peptides", "average precision on the held-out set",
       "recall at a 0.1% false-positive rate (the baseline comparison)",
       "recovered at that point (of 64)",
       "recall at the calling threshold used for the candidate sets",
       "recovered at that point (of 64)"],
      [[m, base.get("n_holdout"),
        round(v["holdout_ap"], 4) if v.get("holdout_ap") is not None else "",
        round(v["holdout_recall"], 4), round(v["holdout_recall"] * 64),
        (round(rep["holdout_annotation"]["recall_at_threshold"], 4)
         if m == "ESM-2 + logistic" else ""),
        (round(rep["holdout_annotation"]["recall_at_threshold"] * 64)
         if m == "ESM-2 + logistic" else "")]
       for m, v in base.get("methods", {}).items()])

sheet("S7", "LoRA fine-tuning against the matched frozen baseline", "paired_tests.json",
      ["fold", "frozen ESM-2 35M + logistic (AP)", "LoRA fine-tuned (AP)"],
      [[i + 1, round(a, 4), round(b, 4)] for i, (a, b) in
       enumerate(zip(paired["finetuned_vs_frozen_same_folds"]["frozen_folds"],
                     paired["finetuned_vs_frozen_same_folds"]["finetuned_folds"]))] +
      [["mean", round(paired["finetuned_vs_frozen_same_folds"]["frozen_mean"], 4),
        round(paired["finetuned_vs_frozen_same_folds"]["finetuned_mean"], 4)],
       ["paired difference",
        round(paired["finetuned_vs_frozen_same_folds"]["paired_difference"], 4),
        f"95% CI {paired['finetuned_vs_frozen_same_folds']['ci_t'][0]:+.4f} to "
        f"{paired['finetuned_vs_frozen_same_folds']['ci_t'][1]:+.4f}"]])

seen = {}
for r in mass:
    seen.setdefault((r["family"], r["core_seq"]), []).append(r)
sheet("S8", "predicted core peptides and their mass ladders",
      "precursor_masses.csv",
      ["family", "cluster", "genomes with this core", "precursor (aa)", "core (aa)",
       "cleavage motif", "core sequence", "Ser+Thr", "monoisotopic core (Da)",
       "average core (Da)", "monoisotopic dehydration ladder (Da)",
       "ladder with C-terminal decarboxylation (Da)", "published peptide",
       "published mass (Da)", "gap, same mass type (Da)"],
      [[k[0], names.fam_short(k[0]), len(v), v[0]["precursor_aa"], v[0]["core_aa"],
        v[0]["cleavage"], k[1], v[0]["n_ser_thr"], v[0]["mono_core"], v[0]["avg_core"],
        v[0]["ladder_mono"], v[0]["ladder_mono_decarb"], v[0]["published"],
        v[0]["published_mass"], v[0]["best_ladder_gap"]] for k, v in seen.items()])

sheet("S9", "DeepBGC 0.1.10 on the same DNA", "deepbgc_results.csv",
      ["region", "assembly", "status", "BGC calls", "calls classed RiPP",
       "best DeepBGC score", "product classes"],
      [[r["region"], r["assembly"], r["status"], r["bgcs"], r["ripp"],
        r["best_score"], r["classes"]] for r in deep])

# ------------------------------------------------------------------ tables
for tag, fn, what in (("Table_1", "table1_panel.tsv", "Table 1. the panel"),
                      ("Table_2", "table2_families.tsv", "Table 2. the cluster families"),
                      ("Table_S1", "tableS1_family_extra.tsv",
                       "Table S1. family location and neighbours"),
                      ("Table_S2", "tableS2_amr.tsv", "Table S2. resistance associations"),
                      ("Table_S3", "tableS3_mlst.tsv", "Table S3. sequence types"),
                      ("Table_S4", "tableS4_primers.tsv", "Table S4. the PCR scheme"),
                      ("Table_S5", "tableS5_masses.tsv", "Table S5. predicted masses")):
    rows = list(csv.reader(open(os.path.join(OUT, fn)), delimiter="\t"))
    sheet(tag, what, fn, rows[0], rows[1:])

# ------------------------------------------------------------------ index
idx = WB.create_sheet("Index")
idx.cell(1, 1, "Source data: one sheet per figure panel and per table").font = Font(bold=True, size=12)
idx.cell(2, 1, "Every sheet holds the values that are actually drawn, recomputed from "
               "the files in out/ that the figure scripts read.").font = Font(italic=True, size=9)
hdr = ["sheet", "what it shows", "computed from", "rows"]
for j, h in enumerate(hdr, start=1):
    idx.cell(4, j, h).font = Font(bold=True)
for i, r in enumerate(INDEX, start=5):
    for j, v in enumerate(r, start=1):
        idx.cell(i, j, v).alignment = Alignment(vertical="top", wrap_text=(j == 2))
idx.column_dimensions["A"].width = 14
idx.column_dimensions["B"].width = 78
idx.column_dimensions["C"].width = 46
idx.column_dimensions["D"].width = 8
WB._sheets.insert(0, WB._sheets.pop(WB._sheets.index(idx)))
dest = os.path.join(OUT, "Source_Data.xlsx")
WB.save(dest)
print("->", dest, f"({len(WB.sheetnames)} sheets)")
