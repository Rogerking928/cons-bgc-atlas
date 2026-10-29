#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tables_and_Figures.docx（HRS 版式：橫向、一頁一項、表在前、圖在後）。

每個圖說裡的數字都從 out/ 讀出來算，不可以手打。
  python3 code/tables.py && python3 code/make_doc.py && python3 code/check_doc.py
"""
import collections, csv, json, os, re, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import docbuild as DB
import names
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")


def tsv(name):
    return [r for r in csv.reader(open(os.path.join(OUT, name)), delimiter="\t")]


def read(name, delim=","):
    return list(csv.DictReader(open(os.path.join(OUT, name)), delimiter=delim))


def js(name):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else {}


# ------------------------------------------------------------------ numbers
md = {r["assembly"]: r for r in read("metadata.csv")}
mat = {r["assembly"]: r for r in read("genome_matrix.csv")}
cl = read("clusters.tsv", "\t")
fams = read("families.tsv", "\t")
fam_of = json.load(open(os.path.join(OUT, "cluster_family.json")))
loc = read("colocalisation.csv")
amr = {r["assembly"]: r for r in read("amr_per_genome.csv")}
mlst = {r["assembly"]: r for r in read("mlst.csv")} if os.path.exists(os.path.join(OUT, "mlst.csv")) else {}
rep = js("precursor_report.json")
hardneg = js("hard_negatives.json")
mutagen = js("mutagenesis.json")
nulldist = js("distance_null.json")
stsens = js("st_sensitivity.json")
en = rep["enrichment_near_machinery"] if rep else {}


def carriage_pair(fam):
    """(panel %, draft %) for a family, rounded exactly as Figure 4A and 4D round."""
    import collections as _c
    sp = _c.Counter(x["species"] for x in mat.values() if x.get(fam)).most_common(1)
    sp = sp[0][0] if sp else ""
    tot = sum(1 for x in mat.values() if x["species"] == sp)
    car = sum(1 for x in mat.values() if x["species"] == sp and x.get(fam))
    d = [x for x in extrows if x["group"] == "draft CoNS" and x["species"] == sp]
    t = [k for k in extj.get("per_target", {}) if k.startswith(fam)]
    card = sum(int(x[f"{t[0]}_carrier"]) for x in d) if t else 0
    return (round(100 * car / max(1, tot)), round(100 * card / max(1, len(d))))
sal_sum = js("finetune_summary.json")
near = sal_sum.get("near_gg", {}) if sal_sum else {}
callrule = js("calling_rule.json")
extpcr = js("external_pcr.json")
extj = extpcr
extrows = read("external_pcr.csv") if os.path.exists(os.path.join(OUT, "external_pcr.csv")) else []
deepbgc = read("deepbgc_results.csv") if os.path.exists(os.path.join(OUT, "deepbgc_results.csv")) else []
fdr = js("fdr_calibration.json")
cal_win = ("%d\u2013%d residues" % tuple(fdr["length_window"])) if fdr else ""
rulebase = js("rule_baseline.json")
paired = js("paired_tests.json")
twins = js("determinant_twins.json")
masses = read("precursor_masses.csv") if os.path.exists(os.path.join(OUT, "precursor_masses.csv")) else []
GAL = next((r for r in masses if r.get("published") == "gallidermin"), {})

# 同重（isobaric）的家族：預測核心肽質量一模一樣，質譜分不開，要在表註講明
_by_mass = collections.defaultdict(set)
for _r in masses:
    _by_mass[_r["mono_core"]].add((_r["family"], _r["core_seq"]))
ISOBARIC = [sorted(v) for v in _by_mass.values()
            if len({f for f, _ in v}) > 1]
gain = js("atlas_gain_summary.json")
base = js("baselines.json")
ladder = js("ladder.json")

N = len(mat)
NSP = len({r["species"] for r in mat.values()})
NBP = len({r["bioproject"] for r in md.values()})
NSRC = sum(1 for a in mat if names.source_group(md[a]["isolation_source"]) != "not stated")
KEEP = [f["family"] for f in fams if int(f["n_assemblies"]) >= 4]
NCL = len(cl)
KINDS = collections.Counter(c["kind"] for c in cl)
CARRY = len({c["assembly"] for c in cl})
PLAS = collections.Counter(c["kind"] for c in cl if c["replicon_type"] != "chromosome")


def fam_n(f):
    return sum(1 for r in mat.values() if r.get(f))


def med_span(f):
    v = [int(c["span"]) for c in cl if fam_of[c["cluster"]] == f]
    return statistics.median(v) / 1000 if v else 0


def near_median(kind_plasmid):
    v = [int(r["nearest_bp"]) for r in loc if r["nearest_bp"] != ""
         and ((r["replicon_type"] != "chromosome") == kind_plasmid)]
    return statistics.median(v) / 1000 if v else float("nan")


def amr_pos(det):
    return sum(1 for a in mat if float(amr[a].get(det, 0) or 0) >= 0.5)


ST = collections.Counter(mlst.get(a, {}).get("ST", "") for a, r in mat.items()
                         if r["species"] == "epidermidis")
ST2 = [a for a, r in mat.items() if r["species"] == "epidermidis"
       and mlst.get(a, {}).get("ST") == "2"]
ST23 = [a for a, r in mat.items() if r["species"] == "epidermidis"
        and mlst.get(a, {}).get("ST") == "23"]
F001_ST2 = sum(1 for a in ST2 if mat[a].get("F001"))
F001_ST23 = sum(1 for a in ST23 if mat[a].get("F001"))

# 同一組呼叫在 0.3／0.7 下的數量（Methods 承諾的敏感度檢查，數字要在展示項裡）
_amr_dets = [k for k in next(iter(amr.values())) if k != "assembly"]
AMR_CALLS = {t: sum(1 for r in amr.values() for d in _amr_dets
                    if r[d] not in ("", None) and float(r[d]) >= t)
             for t in (0.3, 0.5, 0.7)}

def carriage_groups():
    """把八個目標分成不變／變低／變高三組，回傳可直接排進圖說的字串。"""
    fams = sorted({t.split("-")[0] for t in extpcr["per_target"]})
    same, down, up = [], [], []
    for f in fams:
        a_, b_ = carriage_pair(f)
        (same if a_ == b_ else down if b_ < a_ else up).append(
            f"{f} {a_} to {b_}")
    return same, down, up


# 進入距離檢定的染色體座裡，有多少是同一個 F001 座的複本（假重複的規模）
F001_IN_TEST = sum(1 for r in loc if r["family"] == "F001"
                   and r["nearest_bp"] != "" and r["replicon_type"] == "chromosome")
# 每個目標自己的 draft 非帶菌者分母（不是同一個數）
_neg = [int(v["draft_non_carriers"]) for v in extpcr["per_target"].values()]
DRAFT_NEG_LO, DRAFT_NEG_HI = min(_neg), max(_neg)


def fam_species(f):
    return {r["species"] for r in mat.values() if r.get(f)}


CONFINED = [f for f in KEEP if len(fam_species(f)) == 1]
CROSS = [f for f in KEEP if len(fam_species(f)) > 1]
# 兩個跨物種的家族都是質體型的，這是圖 1 標題的宣稱，跑不出來就不要寫
assert all(any(c["replicon_type"] != "chromosome" for c in cl
               if fam_of[c["cluster"]] == f) for f in CROSS), CROSS

items = []

# =============================================================== main items
items.append({
    "kind": "table", "rows": tsv("table1_panel.tsv"), "bold_rows": [0],
    "col_widths": [0.17, 0.09, 0.10, 0.13, 0.08, 0.11, 0.09, 0.11, 0.10, 0.08],
    "caption": f"Table 1. The panel: every complete RefSeq genome of the {NSP} "
               f"coagulase-negative staphylococcal species that colonise and infect "
               f"people, {N} in all.",
    "note": f"Assemblies were taken from NCBI Datasets on 20 September 2026, restricted "
            f"to RefSeq, assembly level complete genome, excluding assemblies flagged as "
            f"atypical; they come from {NBP} BioProjects. Mammaliicoccus sciuri and the "
            f"food and animal species S. xylosus and S. succinus were left out. Isolation "
            f"source is the free-text BioSample field grouped by keyword and is stated "
            f"for {NSRC} of {N} genomes ({100 * NSRC / N:.0f} per cent); the groups are "
            f"descriptive only and no result in this paper rests on them. Resistance "
            f"determinants were called by 10-residue k-mer containment against the "
            f"AMRFinderPlus reference proteins at a containment of 0.5; the program "
            f"itself was not run, so these are containment calls. At 0.3 and 0.7, run as "
            f"sensitivity checks, the panel carries {AMR_CALLS[0.3]:,} and "
            f"{AMR_CALLS[0.7]:,} genome-determinant calls against {AMR_CALLS[0.5]:,} at "
            f"the threshold used here. Figure S1 shows the same composition graphically.",
})

items.append({
    "kind": "table", "rows": tsv("table2_families.tsv"), "bold_rows": [0],
    "col_widths": [0.07, 0.32, 0.09, 0.09, 0.09, 0.22, 0.12],
    "caption": f"Table 2. The {len(KEEP)} cluster families carried by four or more "
               f"genomes.",
    "note": f"A cluster is a run of domain-carrying coding sequences on one replicon, "
            f"chained while neighbours stay within 10 kb, that contains at least one "
            f"maturation or condensation anchor which is not a primary-metabolism gene; "
            f"the reported span is that of the anchor genes, and genes within 5 kb on "
            f"either side are listed with the cluster but do not extend it. Families are "
            f"single-linkage groups at a 10-residue k-mer containment of 0.50. Names are "
            f"taken from the RefSeq product names of the member genes and describe the "
            f"machinery, not a characterised compound. Where each family sits, how often "
            f"it carries a mobility gene and how far it lies from the nearest resistance "
            f"gene are in Table S1 and in Figure 2. Families numbered above F013 are "
            f"carried by fewer than four genomes and are not tabulated here; the ones "
            f"that carry a predicted precursor peptide appear in Table S5, and every "
            f"cluster's family is in the source data.",
})

items.append({
    "kind": "figure", "path": os.path.join(OUT, "Figure1_atlas.png"),
    "caption": f"Figure 1. The atlas: {len(KEEP)} cluster families, {len(CONFINED)} of "
               f"them confined to a single species and the {len(CROSS)} that are not "
               f"carried on plasmids.",
    "note": f"(A) Percentage of each species' genomes carrying at least one cluster of "
            f"each class, for the eight species with 20 or more genomes; a species with "
            f"none is labelled rather than left blank. This panel classes each cluster "
            f"on its own contents and counts all {NCL} clusters. \u00a7\u00a7"
            f"(B) Genomes of each species carrying each family, for the {len(KEEP)} "
            f"families carried by four or more genomes; those families hold "
            f"{sum(1 for c in cl if fam_of[c['cluster']] in KEEP)} of the {NCL} "
            f"clusters, the rest falling in families of one to three genomes that are "
            f"not drawn here. The cell shade is the share of that species' genomes and "
            f"the number is the count, so a pale cell with a large number is a family "
            f"that is common in absolute terms but not universal. The square beside each "
            f"family name is where its clusters sit: chromosome, plasmid or both. A "
            f"family is labelled by the class of most of its clusters, which is not "
            f"always the class of every one: F005 carries a YcaO gene inside the cluster "
            f"window in 25 of its 38 clusters and outside it in the other 13, so panel A "
            f"counts 25 hybrid clusters while the family is called hybrid throughout. "
            f"Species are abbreviated to the first four letters of the epithet: epid. "
            f"S. epidermidis, haem. S. haemolyticus, capi. S. capitis, lugd. "
            f"S. lugdunensis, homi. S. hominis, simu. S. simulans, warn. S. warneri, "
            f"capr. S. caprae. \u00a7\u00a7"
            f"(C) Median and range of the core span; the marker is coloured by the "
            f"family's class on the same scale as panel A, navy for NRPS or PKS, amber "
            f"for RiPP and purple for hybrid. Panels B and C share the family "
            f"axis. The panel calls {NCL} clusters in {CARRY} of {N} genomes: "
            f"{KINDS['NRPS/PKS']} NRPS or PKS, {KINDS['hybrid']} hybrid and "
            f"{KINDS['RiPP']} RiPP.",
})

items.append({
    "kind": "figure", "path": os.path.join(OUT, "Figure2_linkage.png"),
    "caption": "Figure 2. The clusters share replicons with acquired resistance genes; "
               "on chromosomes they sit somewhat closer than chance and on plasmids "
               "farther, and neither survives counting each family once.",
    "note": f"(A, B) Two representative clusters with the 5 kb of neighbours on either "
            f"side, drawn to one common scale, one on a plasmid and one on a chromosome; "
            f"arrows are coloured by the role of the gene and acquired resistance genes "
            f"are black. Both are the member of their family carrying the most resistance "
            f"genes in the window, so they are chosen to show the linkage rather than to "
            f"be typical: the F006 copy drawn here spans "
            f"{int(next(c['span'] for c in cl if c['cluster'] == 'C00077')) / 1000:.1f} "
            f"kb against a family median of {med_span('F006'):.1f} kb in Table 2. The "
            f"remaining maps are in Figure S2. \u00a7\u00a7"
            f"(C) Distance from each cluster to the nearest acquired resistance gene on "
            f"the same replicon, with distances of zero drawn at a 0.1 kb floor; only "
            f"clusters with such a gene on their own replicon enter, "
            f"{nulldist['groups']['chromosome']['n_clusters']} on chromosomes and "
            f"{nulldist['groups']['plasmid']['n_clusters']} on plasmids. The "
            f"dashed line is the median of a null in which the resistance genes of that "
            f"replicon are placed at random {nulldist['n_permutations']} times, keeping "
            f"their number and the replicon length. On chromosomes the observed median "
            f"({nulldist['groups']['chromosome']['observed_median_kb']:.0f} kb) is a "
            f"little below the null "
            f"({nulldist['groups']['chromosome']['null_median_kb']:.0f} kb; Wilcoxon "
            f"signed-rank test on the paired observed and null distances, "
            f"P = {nulldist['groups']['chromosome']['wilcoxon_p']:.3f}); on plasmids it "
            f"is {nulldist['groups']['plasmid']['observed_median_kb']:.0f} kb against a "
            f"null of {nulldist['groups']['plasmid']['null_median_kb']:.0f} kb "
            f"(P = {nulldist['groups']['plasmid']['wilcoxon_p']:.3f}), which is farther "
            f"apart than chance, not closer. The short distances on plasmids are what a "
            f"short replicon gives you, and with "
            f"{nulldist['groups']['plasmid']['n_clusters']} plasmid clusters that P is "
            f"the smallest value the signed-rank test can return here, so the plasmid "
            f"side is a direction rather than a measurement. Those tests treat every "
            f"cluster as independent, and they are not: {F001_IN_TEST} of the "
            f"{nulldist['groups']['chromosome']['n_clusters']} chromosomal clusters in "
            f"the test are the same S. epidermidis locus in different genomes. Counted "
            f"once per family the difference disappears on both "
            f"replicon types (chromosome "
            f"{nulldist['family_level']['chromosome']['n_families']} families, "
            f"P = {nulldist['family_level']['chromosome']['wilcoxon_p']:.2f}; plasmid "
            f"{nulldist['family_level']['plasmid']['n_families']} families, "
            f"P = {nulldist['family_level']['plasmid']['wilcoxon_p']:.2f}), and that is "
            f"the estimate to believe. \u00a7\u00a7"
            f"(D) Share of each family's clusters that share a replicon with at least "
            f"one determinant. \u00a7\u00a7"
            f"(E) Mobility genes inside the cluster window, as a proportion of clusters "
            f"of each class. \u00a7\u00a7"
            f"(F) Cluster GC minus genome GC for every cluster of each family; the RiPP "
            f"families sit several points below their host genome, which is what an "
            f"acquired segment looks like.",
})

items.append({
    "kind": "figure", "path": os.path.join(OUT, "Figure3_lineage.png"),
    "caption": "Figure 3. One line in the cluster-calling rule manufactured a "
               "significant association.",
    "note": f"(A) Carriage of the aureusimine cluster in S. epidermidis by sequence "
            f"type, under three definitions. Under the first version of our calling "
            f"rule, which required at least two domain-carrying coding sequences in the "
            f"window, "
            f"{callrule['carriage']['strict (>=2 domain-carrying CDS)']} of "
            f"{callrule['n_genomes']} genomes counted as carriers and both ST23 and "
            f"ST20 counted as none, which is why two of the grey bars are labelled "
            f"zero; under the rule used in this paper, in which a lone multi-modular "
            f"synthetase is a cluster, it is "
            f"{callrule['carriage']['fixed (a lone megasynthetase counts)']}, which is "
            f"exactly the number of genomes RefSeq annotates with the synthetase at all. "
            f"The gene was always there: ausA is a single 2,387-residue protein and in "
            f"those genomes it happens to have no second domain-carrying gene within "
            f"10 kb. \u00a7\u00a7"
            f"(B) The same odds ratios under each definition. Against erm(A) the first "
            f"rule gives "
            f"{callrule['association']['strict (>=2 domain-carrying CDS)']['erm(A)']['odds']:.2f} "
            f"(P = {callrule['association']['strict (>=2 domain-carrying CDS)']['erm(A)']['p']:.0e}), "
            f"which reads as a strong finding about hospital lineages; the corrected "
            f"rule gives "
            f"{callrule['association']['fixed (a lone megasynthetase counts)']['erm(A)']['odds']:.2f} "
            f"(P = {callrule['association']['fixed (a lone megasynthetase counts)']['erm(A)']['p']:.2f}) "
            f"and there is nothing to report. The artefact was found by the external "
            f"PCR check of Figure 4, where the primer amplified genomes the atlas had "
            f"called non-carriers. \u00a7\u00a7"
            f"(C) Where family carriage does vary within a species: in S. lugdunensis "
            f"the aryl-capped NRPS is in every genome of every sequence type, while the "
            f"lugdunin locus is missing from {40 - fam_n('F006')} of the 40 genomes, "
            f"spread across three of the groups, and the NRPS-YcaO hybrid from "
            f"{40 - fam_n('F005')}, both of them in ST1. The differences are small and "
            f"the groups smaller, so the panel is descriptive. Isolation source, which "
            f"travels with the "
            f"submitting project, is in Figure S3 and is not used for inference.",
})

items.append({
    "kind": "figure", "path": os.path.join(OUT, "Figure4_assay.png"),
    "caption": f"Figure 4. The PCR scheme, and what it is worth once it leaves the "
               f"genomes it was designed on: {extpcr['n_draft']} draft genomes that were "
               f"never part of the design set and {extpcr['n_non_target']} genomes of "
               f"other species.",
    "note": f"(A) Positive predictive value of each target, computed twice: once at the "
            f"carriage seen in the design panel and once at the carriage seen in the "
            f"held-out drafts, with the two carriage figures written beside each row. "
            f"Predictive value depends on how common the family actually is, and "
            f"carriage is lower outside the design panel for "
            f"{len(carriage_groups()[1])} of the {len(extpcr['per_target'])} targets "
            f"({', '.join(carriage_groups()[1])} per cent), so for those the "
            f"panel-based value is the optimistic one; the three largest drops are the "
            f"gallidermin locus of S. capitis (F004), the lugdunin locus (F006) and the "
            f"plasmid lanthipeptide locus (F009). Intervals come from sampling "
            f"sensitivity, specificity "
            f"and prevalence from their Beta posteriors. Negative predictive value is "
            f"undefined where every genome of the species carries the family and is not "
            f"drawn. \u00a7\u00a7"
            f"(B) Sensitivity on the design panel and on the held-out drafts, with "
            f"Wilson intervals; the number beside each row is the number of draft "
            f"carriers. Carriage in the drafts is defined at the protein level, by "
            f"10-residue k-mer containment against the family's marker protein, because "
            f"a draft assembly breaks clusters across contigs. \u00a7\u00a7"
            f"(C) Where the assay should stay silent: draft genomes that are not "
            f"carriers, and genomes of S. aureus, Cutibacterium acnes, Enterococcus "
            f"faecalis, Escherichia coli, Streptococcus mitis and Micrococcus luteus. "
            f"No target amplified any of the {extpcr['n_non_target']} non-target "
            f"genomes, and at most two of the non-carrier drafts, of which each target "
            f"has its own number, {DRAFT_NEG_LO} to {DRAFT_NEG_HI}. "
            f"\u00a7\u00a7"
            f"(D) Family carriage in the drafts against carriage in the design panel, "
            f"within each family's main species, as percentages of that species' "
            f"genomes. Carriage is unchanged for "
            f"{', '.join(carriage_groups()[0])}; it is lower outside the panel for "
            f"{', '.join(carriage_groups()[1])}, so the complete-genome collection "
            f"over-represents those; and higher outside it for "
            f"{', '.join(carriage_groups()[2])}. F004 is the gallidermin locus of "
            f"S. capitis, F006 the lugdunin locus and F009 the plasmid lanthipeptide "
            f"locus. \u00a7\u00a7"
            f"DeepBGC on the same loci is in Figure S9. \u00a7\u00a7"
            f"Primer sequences, melting temperatures, amplicon sizes and raw counts are "
            f"in Table S4; the predicted peptide masses are in Figure S8 and Table S5. "
            f"These primers have not been run in a laboratory.",
})

items.append({
    "kind": "figure", "path": os.path.join(OUT, "Figure5_detector.png"),
    "caption": "Figure 5. A protein language model finds precursor peptides the profile "
               "library misses, and how hard that task really is.",
    "note": f"ESM-2 (150M) mean-pooled embeddings with logistic regression, trained on "
            f"{rep['n_pos']} precursor peptides and cross-validated by sequence cluster "
            f"so a near-duplicate of a training peptide cannot be scored as a held-out "
            f"success. \u00a7\u00a7"
            f"(A, B) Precision-recall against three classical baselines on the same "
            f"folds. The comparison that matters is against the strongest of them, 3-mer "
            f"counts: average precision "
            f"{base['methods']['ESM-2 + logistic']['cv_ap']:.2f} against "
            f"{base['methods']['3-mer counts']['cv_ap']:.2f}, paired difference "
            f"{paired['language_model_vs_3mer']['difference']:.2f} "
            f"(95% CI {paired['language_model_vs_3mer']['ci'][0]:.2f} to "
            f"{paired['language_model_vs_3mer']['ci'][1]:.2f}), bootstrapped by sequence "
            f"cluster. The leader-rule baseline is a rule set of the kind antiSMASH "
            f"and BAGEL use for precursor detection - length, a double-glycine or "
            f"equivalent cleavage motif, a core of 12 to 45 residues, cysteine, serine "
            f"and threonine content, and leader charge - scored 0 to 5; it reaches an "
            f"average precision of "
            f"{rulebase['ap']:.3f} here. It is our implementation of that kind of rule, "
            f"not those tools, because no published peptide-level RiPP predictor "
            f"installs and runs offline on this machine. \u00a7\u00a7"
            f"(C) The same detector against negative sets of increasing difficulty: the "
            f"easy set used for training, short coding sequences that lie inside a "
            f"cluster window but carry no precursor domain, composition-shuffled "
            f"positives, and length-matched short coding sequences. Average precision "
            f"rises with the harder sets only because they are smaller and richer in "
            f"positives, which is why the chance level is drawn for each and the "
            f"prevalence-free AUROC is drawn beside it: that falls from "
            f"{hardneg['negative_sets']['easy (far from machinery)']['auroc']:.3f} on the "
            f"easy set to "
            f"{hardneg['negative_sets']['hard: composition-shuffled positives']['auroc']:.3f} "
            f"against shuffled positives, so part of what the model reads is amino-acid "
            f"composition. That comparison is a peptide against its own shuffle and the "
            f"real peptide still wins; Figure 6E and F run the same shuffle against "
            f"ordinary short proteins instead, where the shuffles win, and the two "
            f"together are what place the weight on composition. A logistic regression "
            f"on sequence length alone is drawn as a triangle on each average-precision "
            f"bar: it reaches "
            f"{hardneg['negative_sets']['easy (far from machinery)']['ap_length_only']:.2f} "
            f"on the easy set, where it is far above that set's chance level because the "
            f"precursors are short, and it stays below the model on every set. \u00a7\u00a7"
            f"(D) Performance by how similar each training peptide is to the next "
            f"nearest one; peptides with less than 30 per cent 10-mer identity to any "
            f"other are still detected. \u00a7\u00a7"
            f"(E) Average precision against the number of training peptides: the curve "
            f"has not flattened at {rep['n_pos']}, which is the most likely reason "
            f"fine-tuning does not help (Figure S7). \u00a7\u00a7"
            f"(F) The same task over the ESM-2 size ladder; 35M is already within "
            f"{ladder['150M']['cv_ap'] - ladder['35M']['cv_ap']:.02f} of 150M.",
})

items.append({
    "kind": "figure", "path": os.path.join(OUT, "Figure6_evidence.png"),
    "caption": "Figure 6. The candidates stand up to a check that owes nothing to the "
               "labels; the explanation of why they do does not; and what the score "
               "actually reads is composition, which is why its q values are not false "
               "discovery rates.",
    "note": f"(A) High scorers with no profile hit sit within 10 kb of cluster machinery "
            f"{100 * en['high_near'] / (en['high_near'] + en['high_far']):.0f} per cent "
            f"of the time against "
            f"{100 * en['rest_near'] / (en['rest_near'] + en['rest_far']):.1f} per cent "
            f"for other short coding sequences (odds ratio {en['odds_ratio']:.0f}; "
            f"permutation P < {1 / en.get('n_permutations', 200):.4f} over "
            f"{en.get('n_permutations', 200):,} draws, in which the null count never "
            f"exceeded {en.get('permutation_null_max', 0)} of the "
            f"{en['high_near'] + en['high_far']} high scorers). Sequences are counted "
            f"once per distinct amino-acid sequence, not once per genome copy, and the "
            f"permutation shuffles positions within each replicon. Genomic position is "
            f"not available to the model, so this is an external check on the score; it "
            f"is not wholly free of the label definition, because the training negatives "
            f"were sampled more than 20 kb from any maturation or export gene. Nothing "
            f"scored here was in the training set, and Figure 5C scores the same "
            f"detector against negatives taken from inside cluster windows, where it "
            f"still separates them. It is a claim "
            f"about enrichment at the top of the ranking, which is not the same claim as "
            f"false discovery control: panels E and F show that the q labels do not "
            f"control the false discovery rate, and both statements can be true at once. "
            f"\u00a7\u00a7"
            f"(B) Clusters with a precursor before and after the model, for the RiPP and "
            f"hybrid families of four or more clusters: a purely non-ribosomal cluster "
            f"has no precursor peptide to find, so the question is undefined for it "
            f"rather than answered zero. Both thresholds are drawn because they do not "
            f"say the same thing: the {gain['F005']['gaining_a_candidate_q05']} hybrid "
            f"clusters of F005 gain a candidate at the strict threshold, while F007, "
            f"F010 and F011 gain theirs only at q < 0.20. \u00a7\u00a7"
            f"(C) What the model appears to rest on: the share of sequences whose single "
            f"most attributed residue lies next to a double-glycine motif, the site the "
            f"C39 peptidase cuts, "
            f"{100 * near['candidate, q < 0.05'][0] / near['candidate, q < 0.05'][1]:.0f} "
            f"per cent of the q < 0.05 candidates against "
            f"{100 * near['negative'][0] / near['negative'][1]:.0f} per cent of the "
            f"negatives. The inset shows why that share, and not the attribution profile "
            f"itself, is the statistic: every group peaks at the motif, negatives "
            f"included. \u00a7\u00a7"
            f"(D) The causal test of the same idea, and it fails. Editing the motif to "
            f"alanines changes the detector's log-odds by "
            f"{mutagen['candidate']['delta_median']:+.2f} in the candidates (median; a "
            f"control edit of two residues elsewhere gives "
            f"{mutagen['candidate']['control_delta_median']:+.2f}) and by "
            f"{mutagen['profile precursor']['delta_median']:+.2f} in the "
            f"profile-detected precursors; inserting a double-glycine into a negative "
            f"raises it by {mutagen['negative']['delta_median']:+.2f} from a median of "
            f"{mutagen['negative']['baseline_median_logodds']:.1f}, nowhere near the "
            f"calling threshold, and inserting two alanines instead raises it by "
            f"{mutagen['negative']['control_delta_median']:+.2f}, so most of even that "
            f"is the two extra residues. The left two groups have the motif removed and "
            f"the right one has it inserted, which is why they are separated. "
            f"Attribution points at the motif; editing the motif does not move the "
            f"score, so we do not claim the detector reads it. \u00a7\u00a7"
            f"(E, F) What it does read. Decoys are made by shuffling residues within a "
            f"sequence, which keeps the amino-acid composition and the length and "
            f"destroys any real precursor, in two strata so the result can be diagnosed "
            f"rather than only reported; real sequences and background decoys share the "
            f"same length window ({cal_win}), which is why the real set here is "
            f"{fdr['n_real_distinct']:,} sequences rather than the "
            f"{len({r['seq'] for r in read('precursor_scan.csv')}):,} of Figure S6. At "
            f"the q < 0.05 threshold "
            f"{100 * fdr['thresholds']['q<0.05']['real_rate']:.2f} per cent of real short "
            f"coding sequences clear it, against "
            f"{100 * fdr['thresholds']['q<0.05']['decoy_background_rate']:.2f} per cent "
            f"of background decoys, so the decoy-based estimator of the false discovery "
            f"proportion comes out at "
            f"{fdr['thresholds']['q<0.05']['empirical_fdp_background']:.1f}. A proportion "
            f"cannot exceed one: an estimate above one means this null is invalid, not "
            f"that the candidates are mostly false, and it is at best a loose upper "
            f"bound. What it does establish is that the q labels cannot be read as false "
            f"discovery rates, and the candidate sets should be treated as a ranked "
            f"list; the evidence that the top of that ranking is enriched for real "
            f"precursors is panel A. The reason the null fails is the third series: "
            f"{100 * fdr['thresholds']['q<0.05']['decoy_from_positives_rate']:.0f} per "
            f"cent of shuffled real precursors still clear the same threshold, so "
            f"shuffling does not destroy what the model reads, and shuffling an ordinary "
            f"protein removes whatever order-borne signal was keeping it low. Order "
            f"still helps at the margin - in Figure 5C the real precursors beat their "
            f"own shuffles, AUROC "
            f"{hardneg['negative_sets']['hard: composition-shuffled positives']['auroc']:.3f} "
            f"- but composition carries the call.",
})

# -------------------------------------------------------------- Figure 7
ascomp = js("antismash_compare.json")
pseudo = js("pseudogene_summary.json")
_as_out = ascomp["out_of_scope_products"]
_as_in = ascomp["antismash_regions_in_scope"]
_as_hit = ascomp["in_scope_regions_matched_by_atlas"]
_f1 = pseudo["families"]["F001"]
_f5 = pseudo["families"]["F005"]
_f6 = pseudo["families"]["F006"]
# 嚴格口徑：扣掉 antiSMASH 的兩個總稱類別後的範圍內 region
_mb_fam_tot = _mb_fam_hits = 0
for _r in csv.DictReader(open(os.path.join(OUT, "mibig_hits.csv"))):
    if _r["atlas_family"]:
        _mb_fam_tot += 1
        _mb_fam_hits += bool(_r["mibig_accession"])
_strict_tot, _strict_hit = 0, 0
for _r in csv.DictReader(open(os.path.join(OUT, "antismash_regions.csv"))):
    if _r["scope"] != "in":
        continue
    _ps = set(_r["product"].split(";"))
    if _ps == {"RiPP-like"}:
        continue
    if "opine-like-metallophore" in _ps and _ps <= {"lanthipeptide-class-v",
                                                    "opine-like-metallophore"}:
        continue
    _strict_tot += 1
    _strict_hit += int(_r["matched_by_atlas"])

items.append({
    "kind": "figure", "path": os.path.join(OUT, "Figure7_antismash.png"),
    "caption": "Figure 7. The same genomes read by antiSMASH: most of what it "
               "calls lies outside what this atlas sets out to find, the two "
               "agree on the loci that lie inside it, and the few "
               "disagreements led to synthetases the annotation had called "
               "pseudogenes.",
    "note": f"antiSMASH 8.0.4 was given the same RefSeq gene calls as the "
            f"atlas, so the two differ in the detector and not in the gene "
            f"prediction. \u00a7\u00a7"
            f"(A) The {ascomp['antismash_regions']:,} regions it calls, by "
            f"product. The four commonest are classes the 45 profiles of this "
            f"atlas do not cover: a terpene precursor and an agr autoinducer "
            f"locus in all {_as_out.get('terpene-precursor', 0)} and "
            f"{_as_out.get('cyclic-lactone-autoinducer', 0)} genomes "
            f"respectively, an NIS-type siderophore in "
            f"{_as_out.get('NI-siderophore', 0)} and a type III polyketide "
            f"synthase in {_as_out.get('T3PKS', 0)}. Shading marks whether a "
            f"product falls within the atlas's scope; the comparison in B uses "
            f"the shaded group only. \u00a7\u00a7"
            f"(B) Agreement in both directions. antiSMASH overlaps "
            f"{ascomp['atlas_clusters_overlapped_by_in_scope_region']} of the "
            f"atlas's {ascomp['atlas_clusters']} clusters; the atlas overlaps "
            f"{_as_hit} of the {_as_in} in-scope regions, or "
            f"{_strict_hit} of {_strict_tot} once antiSMASH's two catch-all "
            f"categories are set aside - a region called RiPP-like alone, and "
            f"a class-v lanthipeptide called together with staphylopine. Both "
            f"figures are given because the catch-alls are where the two "
            f"definitions disagree rather than where either is wrong. \u00a7\u00a7"
            f"MIBiG 4.0 assignment by knownclusterblast, supplied in full "
            f"with the supplementary data, returns a hit for "
            f"{_mb_fam_hits} of the {_mb_fam_tot} regions that match an atlas "
            f"cluster and nothing for the other {_mb_fam_tot - _mb_fam_hits}; "
            f"where it returns one it names the chemistry the family is named "
            f"for, including lugdunin at a median 93 and nukacin ISK-1 at 90 "
            f"per cent of the reference cluster's genes. \u00a7\u00a7"
            f"(C) The boundaries of the same loci. A region is a window around "
            f"the core genes; the atlas reports the span of the anchor genes "
            f"themselves, so the region is wider by construction. \u00a7\u00a7"
            f"(D) What the disagreements led to. For the families in which the "
            f"annotation calls an anchor gene a pseudogene, the number of "
            f"genomes in which the atlas calls a carrier and the number in "
            f"which the locus is there but the gene is pseudogenised. All "
            f"{_f1['genomes_with_pseudogene_but_not_called_carrier']} "
            f"S. epidermidis genomes that the corrected rule counts as "
            f"non-carriers of F001 carry ausA as a pseudogene, so the locus is "
            f"present in "
            f"{_f1['locus_present_including_pseudogenised']} genomes and "
            f"intact in {_f1['atlas_carriers']}; F005 is present in "
            f"{_f5['locus_present_including_pseudogenised']} and intact in "
            f"{_f5['atlas_carriers']}, and the lugdunin locus F006 present in "
            f"{_f6['locus_present_including_pseudogenised']} and intact in "
            f"{_f6['atlas_carriers']}. Across the panel there are "
            f"{pseudo['pseudogene_rows']} pseudogenised anchor genes in "
            f"{pseudo['genomes_affected']} genomes, "
            f"{pseudo['reasons'].get('frameshifted', 0)} of them frameshifted. "
            f"No protein entry is generated for a pseudogene, and this atlas "
            f"reads proteomes, so these are invisible to it.",
})

# ======================================================== supplementary items
items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS1_panel.png"),
    "caption": "Figure S1. The panel, graphically.",
    "note": "(A) Complete RefSeq genomes per species; the four species with fewer than "
            "20 genomes are pooled. (B) Isolation source, grouped by keyword from the "
            "free-text BioSample field. The same numbers are in Table 1.",
})
items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS2_maps.png"),
    "caption": "Figure S2. Gene maps of four more representative clusters, on one "
               "common scale.",
    "note": "All four are drawn at the same kilobases per inch, so a shorter cluster "
            "occupies less of the width. Arrows are coloured by the role of the gene, "
            "acquired resistance genes are black, and where a family has any cluster "
            "with a resistance gene in the window the representative shown is the one "
            "carrying the most of them.",
})
items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS3_source.png"),
    "caption": "Figure S3. Carriage by isolation source, which is descriptive only.",
    "note": "Isolation source travels with the submitting project, so these proportions "
            "cannot be read as an effect of the source; no result in the paper uses "
            "them.",
})
items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS4_architecture.png"),
    "caption": "Figure S4. Domain architecture, adenylation specificity codes and the "
               "GC profile across one locus.",
    "note": "(A) Every domain-carrying protein of at least 400 residues in one "
            "representative cluster of four families, drawn to a common residue scale. "
            "(B) Stachelhaus specificity codes of the adenylation domains, read at the "
            "ten pocket positions after global alignment to GrsA (UniProt P0C061) with "
            "BLOSUM62; the control, GrsA against itself, returns the published "
            "DAWTIAAICK. Each row is one adenylation domain of one synthetase, named by "
            "its protein, because a family can carry three synthetases and each has its "
            "own first domain; only the commonest signature of each is shown, which is "
            "why a row can count 287 of the 288 clusters of a family, the remaining one "
            "carrying a single-residue variant. Seven or eight matching positions out of "
            "ten is a suggestion of substrate, not an assignment. (C) GC in 2 kb "
            "windows around one cluster and the 40 kb on either side. Cluster GC minus "
            "genome GC by family was a fourth panel here and has been removed as a "
            "duplicate: it is Figure 2F.",
})
items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS5_tree.png"),
    "caption": f"Figure S5. Similarity tree of all {N} genomes with gene-content rings.",
    "note": "Average linkage on whole-proteome MinHash distances over 9-residue protein "
            "words, 2000 sketch values per genome. This is a similarity dendrogram, not "
            "a phylogeny; it is the source of the lineage strata used for the "
            "Cochran-Mantel-Haenszel estimates in Table S2. "
            "Branches are coloured by species where all their leaves agree. Rings, from "
            "the inside out: four cluster families, mecA, and candidate precursor "
            "peptides called by the detector.",
})
items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS6_scores.png"),
    "caption": "Figure S6. Detector scores, and the held-out set every method fails.",
    "note": f"(A) Score distribution by distinct sequence. The amber group is the "
            f"held-out annotation set: peptides RefSeq names a bacteriocin or "
            f"lantibiotic but no profile in the library recognises. "
            f"(B) Average precision on that set. Two recalls are quoted in this paper "
            f"and they are measured at different operating points: at the calling "
            f"threshold used for the candidate sets the detector recovers "
            f"{rep['holdout_annotation']['recall_at_threshold'] * 64:.0f} of the 64 "
            f"({100 * rep['holdout_annotation']['recall_at_threshold']:.0f} per cent), "
            f"and at the 0.1 per cent false-positive rate used for the baseline "
            f"comparison it recovers "
            f"{base['methods']['ESM-2 + logistic']['holdout_recall'] * 64:.0f} "
            f"({100 * base['methods']['ESM-2 + logistic']['holdout_recall']:.0f} per "
            f"cent). Every method tested does badly here, and this is reported as a "
            f"failure of the detector, not a success.",
})
items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS7_finetune.png"),
    "caption": "Figure S7. Fine-tuning against the matched frozen baseline.",
    "note": f"The same "
            f"{paired['finetuned_vs_frozen_same_folds']['n_sequences']:,} sequences "
            f"({paired['finetuned_vs_frozen_same_folds']['n_positives']} positives) and "
            f"the same five folds, ESM-2 35M frozen behind a logistic regression and then "
            f"LoRA fine-tuned: average precision "
            f"{paired['finetuned_vs_frozen_same_folds']['frozen_mean']:.3f} against "
            f"{paired['finetuned_vs_frozen_same_folds']['finetuned_mean']:.3f}, "
            f"across-fold paired difference "
            f"{paired['finetuned_vs_frozen_same_folds']['paired_difference']:+.3f} "
            f"(95% CI {paired['finetuned_vs_frozen_same_folds']['ci_t'][0]:+.3f} to "
            f"{paired['finetuned_vs_frozen_same_folds']['ci_t'][1]:+.3f}), which includes "
            f"zero. Lines join the same fold.",
})

items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS8_mass.png"),
    "caption": "Figure S8. Predicted core peptides and the masses a MALDI-TOF spectrum "
               "would show.",
    "note": "(A) Monoisotopic mass ladders of the modal core peptide of each family, "
            "with n the number of clusters carrying that core, "
            "with the published monoisotopic mass marked where one exists and the gap to "
            "the closest rung beside it. (B) Length and serine plus threonine content of "
            "every distinct predicted core peptide; marker area grows with the number of "
            "genomes carrying it. Five of the families here are single-genome families. "
            "Nothing in this figure has been measured; it is a lookup table for a "
            "spectrum, and the rule that generates it is checked against five peptides "
            "whose mature forms are annotated in UniProt (Table S5).",
})
def _spell(k):
    """A caption is prose: "each of eight families", not "each of 8"."""
    words = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight",
             "nine", "ten", "eleven", "twelve", "thirteen"]
    return words[k] if k < len(words) else str(k)


items.append({
    "kind": "figure", "path": os.path.join(OUT, "FigureS9_deepbgc.png"),
    "caption": "Figure S9. A published tool on the same DNA: DeepBGC 0.1.10.",
    "note": f"DeepBGC was run over a 10 kb window either side of one representative "
            f"cluster from each of "
            f"{_spell(sum(1 for r in deepbgc if not r['region'].startswith('control')))} "
            f"families and over one control window of the same length from each of "
            f"those genomes, which is the design the Methods state; those families "
            f"cover the distinct machinery types of the atlas. "
            f"(A) It calls a cluster on "
            f"{sum(1 for r in deepbgc if not r['region'].startswith('control') and int(r['bgcs']) > 0)} "
            f"of {sum(1 for r in deepbgc if not r['region'].startswith('control'))} of "
            f"our loci and on "
            f"{sum(1 for r in deepbgc if r['region'].startswith('control') and int(r['bgcs']) > 0)} "
            f"of {sum(1 for r in deepbgc if r['region'].startswith('control'))} controls. "
            f"(B) Per locus: which ones it calls and which it classes as RiPP. It misses "
            f"the C39 exporter locus of S. haemolyticus and classes only the two plasmid "
            f"lanthipeptide loci as RiPP. No published peptide-level RiPP predictor "
            f"installs and runs offline on this machine, so this comparison is at the "
            f"cluster level; that is a limitation, not a claim of superiority.",
})

items.append({
    "kind": "table", "rows": tsv("tableS4_primers.tsv"), "bold_rows": [0],
    "col_widths": [0.08, 0.05, 0.13, 0.07, 0.15, 0.15, 0.06, 0.08, 0.07, 0.07, 0.05, 0.04],
    "caption": "Table S4. The PCR typing scheme in full.",
    "note": "Primers are 18 to 25-mers exact in at least 90 per cent of the carriers' "
            "copies of the marker gene; melting temperatures are nearest-neighbour "
            "values at 50 mM sodium and 0.25 micromolar primer. Sensitivity and "
            "specificity come from exact-match in-silico PCR on both strands of every "
            "genome in the panel, scored against the cluster calls, so they are measured "
            "on the same genomes the primers were designed on and are optimistic. Exact "
            "matching is also stricter than a real reaction, which tolerates mismatches "
            "away from the 3' end. These primers have not been run in a laboratory.",
})
items.append({
    "kind": "table", "rows": tsv("tableS5_masses.tsv"), "bold_rows": [0],
    "col_widths": [0.05, 0.13, 0.06, 0.05, 0.21, 0.09, 0.09, 0.11, 0.11, 0.10],
    "caption": "Table S5. Predicted mature core peptides and the masses a MALDI-TOF "
               "spectrum would show.",
    "note": "The leader is cut at the last Gly-Gly, Gly-Ala, Gly-Ser, Pro-Gln, Pro-Arg, "
            "Pro-Lys, Ala-Arg, Val-Arg or Glu-Arg motif that leaves a core of 12 to 45 "
            "residues. That rule is checked on every run against five peptides whose "
            "mature forms are annotated in UniProt (gallidermin P21838, epidermin "
            "P08136, Pep5 P19578, epilancin 15X P86047 and nukacin Q9KWM4) and the "
            "pipeline stops if any of them moves. Both mass types are given, because the "
            "published masses of these peptides are monoisotopic and an average mass "
            "runs about 1.7 Da higher at this size: comparing the two would manufacture "
            "a discrepancy of exactly that order. Each lanthionine ring removes one "
            "water, so the ladder gives the core after 0 to n dehydrations, and the "
            "second ladder also subtracts the oxidative decarboxylation of a C-terminal "
            "cysteine. Compared like with like, the closest rung is "
            f"{float(GAL.get('best_ladder_gap', 0)):.2f} Da from the published "
            f"monoisotopic mass of gallidermin and of epidermin, against "
            f"{float(GAL.get('best_ladder_gap_if_average', 0)):.2f} Da if average masses "
            "are compared with monoisotopic ones. "
            + ("Mass alone does not separate every family: "
               + "; ".join(" and ".join(f"{f} ({q})" for f, q in grp)
                           + " differ only by leucine against isoleucine, which are "
                             "isomers, so their predicted cores have the same mass and "
                             "the same ladder"
                           for grp in ISOBARIC)
               + " - the same relation gallidermin and epidermin have to each other. "
                 "Telling those apart needs the sequence or the PCR scheme, not the "
                 "spectrum. " if ISOBARIC else "")
            + f"The count column counts clusters, not genomes: one genome can carry "
              f"two clusters of a family, which is why the F009 rows add up to "
              f"{sum(1 for r in masses if r['family'] == 'F009')} clusters in "
              f"{len({r['assembly'] for r in masses if r['family'] == 'F009'})} genomes. "
              f"A core is predicted for "
              f"{sum(1 for r in masses if r['family'] == 'F004')} of the "
              f"{fam_n('F004')} F004 clusters; the remaining precursor carries no "
              f"cleavage motif from the list above that leaves a core of 12 to 45 "
              f"residues, so it is left without a prediction rather than given a guessed "
              f"one. These are lookup windows for a spectrum, not measurements.",
})
items.append({
    "kind": "table", "rows": tsv("tableS2_amr.tsv"), "bold_rows": [0],
    "col_widths": [0.08, 0.04, 0.10, 0.07, 0.06, 0.06, 0.07, 0.05, 0.05, 0.05,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.05],
    "caption": "Table S2. Cluster families against acquired resistance determinants, "
               "crude and after stratification by lineage.",
    "note": "Tests are within species and only where at least 20 genomes carry the "
            "family and at least 10 carry the determinant. The table lists the four "
            "pairs whose crude test reaches P < 0.05; the BH q column is the "
            "Benjamini-Hochberg value across all tests of that species, and none of "
            "them reaches 0.05, so no association here survives correction for multiple "
            "testing. The two count columns give "
            "carriers of the family with the determinant over all carriers, and the same "
            "for non-carriers. Two kinds of cell are written n.e., and the bracket tells "
            "them apart. In the crude odds ratio column, n.e. followed by a figure in "
            "brackets is a table with a zero cell, where the odds ratio is not "
            "estimable and the bracket gives the Haldane-Anscombe estimate instead. "
            "In the stratified columns, n.e. stands alone. Lineage "
            "strata are average-linkage clusters of the whole-proteome MinHash distances "
            "at two cut-offs; the number of informative strata is reported because a "
            "stratum in which every genome carries the family, or none does, contributes "
            "nothing, and it is zero for the S. capitis rows; at the third cut-off "
            "scanned, 0.005, it is zero for those rows too and one for the S. epidermidis "
            "row, which is why only two cut-offs are tabulated. A stratified estimate "
            "written n.e., without a bracket, is one that a single informative stratum "
            "has driven to zero or infinity; it is a degenerate pooled estimate rather "
            "than a zero cell, and no Haldane-Anscombe substitute is offered for it "
            "because the problem is the number of strata and not the size of a count. "
            "The flag column marks rows "
            "whose crude association does not survive the removal of one BioProject. "
            "Determinants that are perfectly co-carried in this panel give identical rows: "
            + "; ".join(f"{a_} and {b_} in all {n_} carriers"
                        for a_, b_, n_ in twins.get("identical_determinant_pairs", [])
                        if a_ not in ("mecA",) and b_ not in ("qacB",)) +
            ". Reference entries for the same gene (" +
            ", ".join(twins.get("duplicate_reference_entries_dropped", [])) +
            ") were dropped so that one gene is not counted twice.",
})
items.append({
    "kind": "table", "rows": tsv("tableS3_mlst.tsv"), "bold_rows": [0],
    "col_widths": [0.14, 0.11, 0.08, 0.37, 0.10, 0.10, 0.10],
    "caption": "Table S3. Sequence types with four or more genomes, and what they carry.",
    "note": "Seven-locus schemes were taken from PubMLST (S. epidermidis, "
            "S. haemolyticus, S. hominis) and from the Institut Pasteur BIGSdb "
            "(S. lugdunensis) and each genome was typed by exact allele match on either "
            "strand. Unauthenticated access to both databases is capped at records "
            "submitted up to 31 December 2024, so a genome whose alleles are newer comes "
            "back untyped rather than mistyped; every allele of the S. capitis scheme is "
            "newer than that cap, so that species is not typed here at all.",
})
items.append({
    "kind": "table", "rows": tsv("tableS1_family_extra.tsv"), "bold_rows": [0],
    "col_widths": [0.08, 0.24, 0.10, 0.14, 0.18, 0.16, 0.10],
    "caption": "Table S1. Where each family sits and what it sits next to.",
    "note": "The columns cut from Table 2, kept here so the main table fits the page. "
            "Distance to the nearest acquired resistance gene is measured on the same "
            "replicon only, and Figure 2C shows what that distance is worth against a "
            "null in which the resistance genes are placed at random.",
})

# ------------------------------------------------------- house post-processing
# 沿用腸球菌那篇的版式：每項一頁、過長的表自己插分頁、被壓扁的圖給一整頁、
# 長圖說按分格切段、三線黑線、砍掉尾端空段。

ABBR_ALL = {
    "AP": "average precision",
    "AUPRC": "area under the precision-recall curve",
    "CDS": "coding sequence",
    "CI": "confidence interval",
    "CMH": "Cochran-Mantel-Haenszel",
    "ESM-2": "evolutionary scale modelling 2 protein language model",
    "GC": "guanine-plus-cytosine content",
    "IS": "insertion sequence",
    "LoRA": "low-rank adaptation",
    "MALDI-TOF": "matrix-assisted laser desorption/ionisation time-of-flight",
    "MLST": "multilocus sequence typing",
    "NRPS": "non-ribosomal peptide synthetase",
    "OR": "odds ratio",
    "PCP": "peptidyl carrier protein",
    "PCR": "polymerase chain reaction",
    "PKS": "polyketide synthase",
    "RiPP": "ribosomally synthesised and post-translationally modified peptide",
    "SCC": "staphylococcal cassette chromosome",
    "ST": "sequence type",
    "TE": "thioesterase",
}


def abbr(*blocks):
    """只列這一項自己用到的縮寫（表下面放沒用到的定義是雜訊）。"""
    hay = " ".join(str(b) for b in blocks)
    used = [f"{k}, {v}" for k, v in ABBR_ALL.items()
            if re.search(r"(?<![A-Za-z0-9])" + re.escape(k) + r"(?![A-Za-z0-9])", hay)]
    return ("Abbreviations: " + "; ".join(used) + ".") if used else ""


def mark_panels(note, mark="\u00a7\u00a7"):
    """長圖說在分格處斷段，讀者才找得到 (D) 在哪。"""
    return re.sub(r"(?<=[.]) (\([A-H]\) )", mark + r"\1", note)


def page_break_items(path):
    """過長的表會把圖說塞進表頭那一列，分頁記號跟著一起被吃掉；
    在表前面插一段 1 pt 的空段落帶分頁，Word 才會從新的一頁開始。"""
    doc = Document(str(path))
    body = doc.element.body
    added = 0
    for tbl in body.findall(qn("w:tbl")):
        first = tbl.find(qn("w:tr"))
        if first is None or not "".join(first.itertext()).startswith(("Table ", "Figure ")):
            continue
        prev = tbl.getprevious()
        if prev is not None and prev.tag == qn("w:p"):
            if prev.findall(f".//{qn('w:br')}[@{qn('w:type')}='page']"):
                continue
            pPr_prev = prev.find(qn("w:pPr"))
            if pPr_prev is not None and pPr_prev.find(qn("w:sectPr")) is not None:
                continue
        par = OxmlElement("w:p")
        pPr = OxmlElement("w:pPr")
        rPr = OxmlElement("w:rPr")
        sz = OxmlElement("w:sz"); sz.set(qn("w:val"), "2")
        rPr.append(sz); pPr.append(rPr)
        spacing = OxmlElement("w:spacing")
        spacing.set(qn("w:after"), "0"); spacing.set(qn("w:before"), "0")
        pPr.insert(0, spacing)
        par.append(pPr)
        r = OxmlElement("w:r")
        rPr2 = OxmlElement("w:rPr")
        sz2 = OxmlElement("w:sz"); sz2.set(qn("w:val"), "2"); rPr2.append(sz2)
        r.append(rPr2)
        br = OxmlElement("w:br"); br.set(qn("w:type"), "page"); r.append(br)
        par.append(r)
        tbl.addprevious(par)
        added += 1
    doc.save(str(path))
    return added


def squeezed_figures(min_width=7.0):
    """圖說太長把圖壓到某個寬度以下的，挑出來給一整頁。"""
    out = []
    for it in items:
        if it["kind"] != "figure":
            continue
        if DB.fit_figure_width(it) < min_width:
            out.append(it["caption"].split(".")[0] + ".")
    return tuple(out)


def full_page_figures(path, captions):
    """指定的圖：圖說留在前一頁，圖自己佔一頁，放到版面容得下的最大尺寸。"""
    from docx.shape import InlineShape
    doc = Document(str(path))
    pars = doc.paragraphs
    done = []
    for i, par in enumerate(pars):
        if not any(par.text.startswith(c) for c in captions):
            continue
        pic = next((p for p in pars[i + 1:i + 5]
                    if p._p.findall(".//" + qn("wp:inline"))), None)
        if pic is None:
            continue
        inline = pic._p.findall(".//" + qn("wp:inline"))[0]
        shape = InlineShape(inline)
        ratio = shape.height / shape.width
        width = min(int(DB.IMG_MAX_W * 914400),
                    int((DB.TEXT_H - 0.12) * 914400 / ratio))
        shape.width, shape.height = width, int(width * ratio)
        DB._page_break_before(pic)
        done.append((par.text.split(".")[0], round(width / 914400, 2)))
    doc.save(str(path))
    return done


def split_notes(path, mark="\u00a7\u00a7"):
    doc = Document(str(path))
    n = 0
    for par in list(doc.paragraphs):
        if mark not in par.text:
            continue
        parts = [p.strip() for p in par.text.split(mark)]
        src = par.runs[0]
        par.runs[0].text = parts[0]
        for extra in par.runs[1:]:
            extra.text = ""
        prev = par
        for chunk in parts[1:]:
            new = doc.add_paragraph()
            r = new.add_run(chunk)
            r.font.size = src.font.size
            r.font.italic = src.font.italic
            r.font.name = src.font.name
            new.paragraph_format.space_after = par.paragraph_format.space_after
            new.paragraph_format.space_before = Pt(3)
            prev._p.addnext(new._p)
            prev = new
            n += 1
    doc.save(str(path))
    return n


def rule_tables(path):
    """三線黑線：表頭上下各一條，最後一列下面一條（docbuild 產的是無框線表）。"""
    def borders(*edges):
        b = OxmlElement("w:tcBorders")
        for e in edges:
            el = OxmlElement(f"w:{e}")
            el.set(qn("w:val"), "single"); el.set(qn("w:sz"), "8")
            el.set(qn("w:space"), "0"); el.set(qn("w:color"), "000000")
            b.append(el)
        return b

    def rule(tr, *edges):
        for tc in tr.findall(qn("w:tc")):
            tcPr = tc.find(qn("w:tcPr"))
            if tcPr is None:
                tcPr = OxmlElement("w:tcPr"); tc.insert(0, tcPr)
            for old in tcPr.findall(qn("w:tcBorders")):
                tcPr.remove(old)
            tcPr.append(borders(*edges))

    d = Document(str(path))
    n = 0
    for tbl in d.element.body.findall(qn("w:tbl")):
        rows = tbl.findall(qn("w:tr"))
        if not rows:
            continue
        head = 1 if (len(rows) > 1 and len(rows[0].findall(qn("w:tc"))) == 1
                     and len(rows[1].findall(qn("w:tc"))) > 1) else 0
        last = len(rows) - 1
        rule(rows[head], "top", "bottom")
        if head + 1 <= last:
            rule(rows[head + 1], *(("top", "bottom") if head + 1 == last else ("top",)))
        if last > head + 1:
            rule(rows[last], "bottom")
        n += 1
    d.save(str(path))
    return n


def trim_trailing_blanks(path):
    doc = Document(str(path))
    body = doc.element.body
    n = 0
    for el in reversed(list(body)):
        if el.tag != qn("w:p"):
            break
        if "".join(el.itertext()).strip():
            break
        pPr = el.find(qn("w:pPr"))
        if pPr is not None and pPr.find(qn("w:sectPr")) is not None:
            break
        body.remove(el)
        n += 1
    if n:
        doc.save(str(path))
    return n



TITLE = ("How a public-genome atlas manufactures findings: the checks that overturned "
         "our own results in the biosynthetic gene clusters of the coagulase-negative "
         "staphylococci")


def front_matter():
    from docx import Document as _Doc
    from docx.enum.text import WD_ALIGN_PARAGRAPH

    doc = _Doc()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(11)

    def p(text, size=11, bold=False, after=6):
        par = doc.add_paragraph()
        par.alignment = WD_ALIGN_PARAGRAPH.LEFT
        par.paragraph_format.space_after = Pt(after)
        r = par.add_run(text)
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.name = "Times New Roman"
        return par

    p("Tables and Figures", size=15, bold=True, after=10)
    p(TITLE, size=12, bold=True, after=10)
    p(f"Public complete genomes only: {N} RefSeq assemblies of {NSP} "
      f"coagulase-negative staphylococcal species, retrieved on 20 September 2026, "
      f"with {extpcr['n_draft']} draft genomes and {extpcr['n_non_target']} genomes of "
      f"other species held out for validation. No isolates, sequencing or phenotypes of "
      f"our own.", size=11, after=8)
    p(f"What held up: {len(CONFINED)} of the {len(KEEP)} cluster families are confined "
      "to a single species and the two that are not are plasmid-borne; two families "
      "- the aureusimine locus of S. epidermidis (F001) and the aryl-capped NRPS "
      "of S. lugdunensis (F003) - carry at the same rate in 480 draft genomes that took "
      "no part in the design (Figure 4D); the PCR scheme keeps its sensitivity outside "
      "the design set, amplifies none of 72 genomes of other species and at most two of "
      "about 440 non-carrier drafts (Figure 4B, C); and the detector's candidates are "
      "enriched near cluster machinery, 68 per cent against 1.2 per cent, on a check "
      "that uses no labels (Figure 6A).", size=10.5, after=8)
    p("Four of our own results did not survive their own checks: the clusters do not "
      "sit closer to resistance genes than chance once each family is counted once "
      "(Figure 2C); the one significant association with resistance was made by a line "
      "in our cluster-calling rule (Figure 3); the carriage the design panel reports is "
      "higher than carriage outside it, and the predictive values follow (Figure 4A, D); "
      "and the detector's q values are not false discovery rates, because what it reads "
      "is amino-acid composition (Figure 6D-F).", size=10.5, after=14)
    p("Contents", size=12, bold=True, after=6)
    for it in items:
        head, rest = it["caption"].split(".", 1)
        p(" " + head + ".  " + rest.strip(), size=9.5, after=2)
    return doc


def insert_divider(path, first_supp="Figure S1."):
    """補充材料的分隔頁：DB.build 不認得分隔，所以事後插進去，並把分頁記號
    從原本那一項移到分隔頁上（兩邊都留會多一張空白頁）。"""
    doc = Document(str(path))
    target = next((p for p in doc.paragraphs if p.text.startswith(first_supp)), None)
    if target is None:
        for tbl in doc.tables:
            for c in tbl.rows[0].cells:
                if c.text.startswith(first_supp):
                    target = tbl
                    break
    if target is None:
        raise SystemExit(f"divider anchor not found: {first_supp}")
    el = target._p if hasattr(target, "_p") else target._tbl
    if hasattr(target, "_p"):
        pPr = target._p.find(qn("w:pPr"))
        if pPr is not None:
            for br in pPr.findall(qn("w:pageBreakBefore")):
                pPr.remove(br)

    def before(text, size, bold):
        par = doc.add_paragraph()
        r = par.add_run(text)
        r.font.size = Pt(size); r.font.bold = bold
        r.font.name = "Times New Roman"
        par.paragraph_format.space_after = Pt(6)
        el.addprevious(par._p)
        return par

    head = before("Supplementary Material", 14, True)
    DB._page_break_before(head)
    DB._keep_next(head)
    sub = before("The panel described graphically, the remaining gene maps, carriage by "
                 "isolation source, the domain architecture and specificity codes, the "
                 "similarity tree, the detector's score distribution and held-out set, "
                 "the fine-tuning control, and the full versions of the tables that were "
                 "cut to fit the page.", 11, False)
    DB._keep_next(sub)
    doc.save(path)


def keep_caption_with_item(path):
    """圖說到圖之間的每一段都掛 keepNext，Word 才不會從中間斷頁。

    split_notes() 把長圖說拆成一 panel 一段之後，只有第一段留著 docbuild 設的
    keepNext，後面幾段與圖之間就可以自由斷頁——Figure 5 的圖說因此落在前一頁、
    圖落在下一頁。這裡補上，從圖說起算到圖的前一段為止。
    """
    BLIP = "{http://schemas.openxmlformats.org/drawingml/2006/main}blip"
    doc = Document(str(path))
    body = list(doc.element.body)
    cap = re.compile(r"^(Figure|Table) S?\d+\.")
    n, run_start = 0, None
    for i, el in enumerate(body):
        if el.tag.endswith("}tbl"):
            run_start = None
            continue
        if not el.tag.endswith("}p"):
            continue
        from docx.text.paragraph import Paragraph
        if el.find(".//" + BLIP) is not None:      # 圖本身：這一串到此為止
            run_start = None
            continue
        text = Paragraph(el, doc).text.strip()
        if cap.match(text):
            run_start = i
        if run_start is not None and text:
            pPr = el.find(qn("w:pPr"))
            if pPr is None or pPr.find(qn("w:keepNext")) is None:
                DB._keep_next(Paragraph(el, doc))
                n += 1
    doc.save(str(path))
    return n


def break_before_first_item(path, first="Table 1."):
    """投稿檔裡第一項前面補一個分頁：DB.build 只從第二項開始加，
    否則 Table 1 會接在 Figure legends 那一節後面、跟文字共用一頁。"""
    doc = Document(str(path))
    for p_ in doc.paragraphs:
        if p_.text.strip().startswith(first):
            DB._page_break_before(p_)
            doc.save(str(path))
            return True
    for t in doc.tables:
        for c in t.rows[0].cells:
            if c.text.strip().startswith(first):
                return False
    return False


def fits_on_one_page(items):
    """回報哪些項目的「圖說＋圖」高過版心，兩者一定會被 Word 拆到兩頁。"""
    out = []
    for it in items:
        if DB.estimate_height(it) > DB.TEXT_H:
            out.append(it["caption"].split(".")[0])
    return out


def sort_key(it):
    m = re.match(r"(Table|Figure) (S?)(\d+)\.", it["caption"])
    kind, supp, num = m.group(1), m.group(2), int(m.group(3))
    if supp:                      # 補充：先圖後表
        return (1, 0 if kind == "Figure" else 1, num)
    return (0, 0 if kind == "Table" else 1, num)


if __name__ == "__main__":
    items.sort(key=sort_key)
    for it in items:
        cells = " ".join(" ".join(str(c) for c in row) for row in it.get("rows", []))
        line = abbr(it["caption"], it.get("note", ""), cells)
        if line:
            it["note"] = (it.get("note", "").rstrip() + " " + line).strip()

    dest = os.path.join(OUT, "Tables_and_Figures.docx")
    DB.build(items, dest, doc=front_matter())
    print(f"page breaks added before {page_break_items(dest)} over-long tables")
    # full_page_figures() 會把圖說留在前一頁、圖自己佔一頁。作者要的是圖說與圖
    # 同一頁，而且 22 項裡只有 Figure 6 真的放不下（圖說本身 5.2 in／版心 6.7 in），
    # 其餘 21 項本來就容得下，所以不再呼叫它；Word 的 keepNext 會把兩者綁在一起。
    over = fits_on_one_page(items)
    if over:
        print("圖說與圖無法同頁（圖說本身太長）: " + ", ".join(over))
    print(f"legends split into {split_notes(dest)} extra paragraphs")
    insert_divider(dest)
    print(f"three-line rules added to {rule_tables(dest)} tables")
    import italics
    # legends roman, italics only where they mean something
    print(f"legend runs set roman: {italics.roman_everything(dest)}; "
          f"gene symbols and binomials italicised in "
          f"{italics.italicise_docx(dest)} runs")
    # 一定要在 split_notes() 之後：拆出來的圖說段落是那一步才產生的，
    # 先補 keepNext 等於沒補（2026-09-21 踩過）
    print(f"keepNext 補在 {keep_caption_with_item(dest)} 段圖說上")
    trimmed = trim_trailing_blanks(dest)
    if trimmed:
        print(f"removed {trimmed} empty paragraphs from the end")
    print("->", dest, f"({len(items)} items)")
