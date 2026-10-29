#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""What the detector adds to the atlas.

For every called cluster, and in particular for the families that were defined
only by an export module, ask whether a high-scoring candidate precursor -
annotated or unannotated - lies inside the cluster or within 5 kb of it.
"""
import collections, csv, json, os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
FLANK = 5000


def main():
    rep = json.load(open(os.path.join(OUT, "precursor_report.json")))
    thr = rep["threshold_999th_percentile_of_negatives"]
    # a second, looser tier: candidates are reported at two false-discovery
    # rates rather than at one, because at q < 0.05 the detector calls almost
    # nothing new and at q < 0.20 it calls a set worth following up
    thr20 = rep.get("threshold_q20", thr)
    fam = json.load(open(os.path.join(OUT, "cluster_family.json")))
    cl = list(csv.DictReader(open(os.path.join(OUT, "clusters.tsv")), delimiter="\t"))
    # a profile-detected precursor counts as present whatever the model scores
    # it; only the novel candidates are subject to the calling threshold
    scan = [r for r in csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv")))
            if r["known_domain"] or float(r["score"]) >= thr20]
    orf = []
    p = os.path.join(OUT, "orf_candidates.csv")
    if os.path.exists(p):
        orf = [r for r in csv.DictReader(open(p)) if float(r["score"]) >= thr]

    by = collections.defaultdict(list)
    for r in scan:
        by[(r["assembly"], r["replicon"])].append(
            (int(r["start"]), int(r["end"]), "annotated", r["protein"],
             float(r["score"]), r["seq"], r["known_domain"]))
    # unannotated ORFs are kept in a separate layer: they are reported, but
    # they are not counted as a gain to the atlas, because the positional test
    # in s21 does not support them (see out/atlas_gain_summary.json, _orf)
    orf_by = collections.defaultdict(list)
    for r in orf:
        orf_by[(r["assembly"], r["replicon"])].append(
            (int(r["start"]), int(r["end"]), "unannotated ORF",
             f"{r['replicon']}:{r['start']}-{r['end']}({r['strand']})",
             float(r["score"]), r["seq"], ""))

    rows = []
    for c in cl:
        lo, hi = int(c["start"]) - FLANK, int(c["end"]) + FLANK
        hits = [h for h in by.get((c["assembly"], c["replicon"]), [])
                if h[0] <= hi and h[1] >= lo]
        orf_hits = [h for h in orf_by.get((c["assembly"], c["replicon"]), [])
                    if h[0] <= hi and h[1] >= lo]
        known = [h for h in hits if h[6]]
        novel = [h for h in hits if not h[6]]
        novel05 = [h for h in novel if h[4] >= thr]
        best = max(novel, key=lambda h: h[4]) if novel else None
        rows.append(dict(
            cluster=c["cluster"], family=fam[c["cluster"]], assembly=c["assembly"],
            replicon_type=c["replicon_type"],
            had_known_precursor=int(bool(known)),
            candidates=len(novel),
            candidates_q05=len(novel05),
            best_score=best[4] if best else "",
            best_kind=best[2] if best else "",
            best_id=best[3] if best else "",
            best_seq=best[5] if best else "",
            unannotated_orf_candidates=len(orf_hits)))
    with open(os.path.join(OUT, "atlas_gain.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    summary = {}
    for f in sorted({r["family"] for r in rows}):
        sub = [r for r in rows if r["family"] == f]
        summary[f] = {
            "clusters": len(sub),
            "with_known_precursor": sum(r["had_known_precursor"] for r in sub),
            "gaining_a_candidate": sum(1 for r in sub
                                       if not r["had_known_precursor"] and r["candidates"]),
            "gaining_a_candidate_q05": sum(1 for r in sub if not r["had_known_precursor"]
                                           and r["candidates_q05"]),
            "unannotated_orf_candidates_nearby": sum(r["unannotated_orf_candidates"] for r in sub),
        }
    # an independent corroboration that never entered the model: what does the
    # RefSeq product name say about the high-scoring proteins the profiles missed?
    import re as _re
    BACT = _re.compile(r"bacteriocin|lantibiotic|lanthipeptide|cin [A-Z]|"
                       r"microcin|autoinducer|pheromone", _re.I)
    allscan = list(csv.DictReader(open(os.path.join(OUT, "precursor_scan.csv"))))
    hi_u, bg_u = {}, {}
    for r in allscan:
        if r["known_domain"]:
            continue
        (hi_u if float(r["score"]) >= thr else bg_u).setdefault(r["seq"], r)
    named = [r for r in hi_u.values() if BACT.search(r["product"])]
    unnamed = [r for r in hi_u.values() if not BACT.search(r["product"])]
    bg_named = sum(1 for r in bg_u.values() if BACT.search(r["product"]))

    def gg(seq):
        return bool(_re.search(r"G[GA]", seq[5:45]))
    summary["_annotated"] = {
        "distinct_above_threshold": len(hi_u),
        "named_as_bacteriocin_or_lantibiotic": len(named),
        "not_named": len(unnamed),
        "named_share_of_the_rest": round(bg_named / max(1, len(bg_u)), 5),
        "double_glycine_motif_hits": sum(1 for r in hi_u.values() if gg(r["seq"])),
        "double_glycine_motif_rest": round(
            sum(1 for r in bg_u.values() if gg(r["seq"])) / max(1, len(bg_u)), 4),
    }
    print(f"  of {len(hi_u)} distinct high scorers with no profile hit, "
          f"{len(named)} are already named a bacteriocin or lantibiotic by RefSeq "
          f"(against {100*bg_named/max(1,len(bg_u)):.2f}% of the rest) and "
          f"{len(unnamed)} are unnamed")

    # how the unannotated candidates sit relative to every ORF in the same
    # windows, and how many the threshold is expected to throw up by chance
    if os.path.exists(p):
        allorf = list(csv.DictReader(open(p)))
        from scipy import stats
        hi_d = [int(r["nearest_anchor_bp"]) for r in allorf if float(r["score"]) >= thr]
        bg_d = [int(r["nearest_anchor_bp"]) for r in allorf if float(r["score"]) < thr]
        u, pu = stats.mannwhitneyu(hi_d, bg_d, alternative="less") if hi_d else (0, 1)
        uniq = len({r["seq"] for r in allorf})
        summary["_orf"] = {
            "orfs": len(allorf), "distinct": uniq,
            "above_threshold": len(hi_d),
            "distinct_above_threshold": len({r["seq"] for r in allorf
                                             if float(r["score"]) >= thr}),
            "expected_by_chance_at_threshold": round(uniq * 0.001, 1),
            "median_bp_to_machinery_hits": float(np.median(hi_d)) if hi_d else None,
            "median_bp_to_machinery_rest": float(np.median(bg_d)) if bg_d else None,
            "mannwhitney_p": float(pu)}
        print(f"  unannotated ORFs: {len(allorf):,} ({uniq:,} distinct); "
              f"{len(hi_d)} above threshold, about "
              f"{uniq*0.001:.0f} expected by chance; median distance to machinery "
              f"{np.median(hi_d):,.0f} bp vs {np.median(bg_d):,.0f} bp, P = {pu:.1e}")

    json.dump(summary, open(os.path.join(OUT, "atlas_gain_summary.json"), "w"), indent=1)
    for f, d in summary.items():
        if f.startswith("_"):
            continue
        if d["clusters"] >= 5:
            print(f"  {f}: {d['clusters']:4d} clusters, {d['with_known_precursor']:4d} already had "
                  f"a precursor, {d['gaining_a_candidate']:4d} gain one at q<0.2 "
                  f"({d['gaining_a_candidate_q05']} at q<0.05) "
                  f"({d['unannotated_orf_candidates_nearby']} unannotated ORFs nearby)")
    print("-> out/atlas_gain.csv")


if __name__ == "__main__":
    main()
