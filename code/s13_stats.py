#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Association tests, all of them run twice: once on the whole panel and once
stratified by lineage, because in RefSeq the sampling structure (clone and
submitting project) is confounded with almost every label.

  1. family x species
  2. family x family, within species
  3. family x acquired-resistance determinant, within species
     + lineage-stratified Cochran-Mantel-Haenszel
     + leave-one-BioProject-out sensitivity
"""
import collections, csv, json, os, sys
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
AMR_CALL = 0.5
LINEAGE_CUTS = (0.005, 0.010, 0.020)


def bh(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    q = np.empty_like(ps)
    m = len(ps)
    prev = 1.0
    for rank, i in enumerate(o[::-1]):
        k = m - rank
        prev = min(prev, ps[i] * m / k)
        q[i] = prev
    return q


def fisher(a, b, c, d):
    return stats.fisher_exact([[a, b], [c, d]])


def woolf_ci(a, b, c, d):
    """粗 OR 的 95% 區間（Woolf，四格各加 0.5 的 Haldane 修正）。"""
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    lor = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return float(np.exp(lor - 1.96 * se)), float(np.exp(lor + 1.96 * se))


def cmh_ci(tables):
    """CMH 的 95% 區間，Robins-Breslow-Greenland 變異數。"""
    tabs = [t for t in tables if (t[0] + t[1]) and (t[2] + t[3])
            and (t[0] + t[2]) and (t[1] + t[3])]
    R = S = 0.0
    PR = PS_QR = QS = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        if n == 0:
            continue
        r = a * d / n
        s_ = b * c / n
        P = (a + d) / n
        Q = (b + c) / n
        R += r; S += s_
        PR += P * r
        PS_QR += P * s_ + Q * r
        QS += Q * s_
    if R == 0 or S == 0:
        return None, None
    var = PR / (2 * R * R) + PS_QR / (2 * R * S) + QS / (2 * S * S)
    lor = np.log(R / S)
    se = np.sqrt(var)
    return float(np.exp(lor - 1.96 * se)), float(np.exp(lor + 1.96 * se))


def cmh(tables):
    """Cochran-Mantel-Haenszel, strata with no information drop out."""
    tabs = [t for t in tables if (t[0] + t[1]) and (t[2] + t[3])
            and (t[0] + t[2]) and (t[1] + t[3])]
    if not tabs:
        return None, None, 0
    num = den = 0.0
    E = V = S = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        num += a * d / n
        den += b * c / n
        E += (a + b) * (a + c) / n
        V += ((a + b) * (c + d) * (a + c) * (b + d)) / (n * n * (n - 1)) if n > 1 else 0
        S += a
    if den == 0 or V == 0:
        return None, None, len(tabs)
    chi = (abs(S - E) - 0.5) ** 2 / V
    return num / den, float(stats.chi2.sf(chi, 1)), len(tabs)


def lineages(accs_keep, cut):
    D = np.load(os.path.join(OUT, "tree_D.npy"))
    meta = json.load(open(os.path.join(OUT, "tree_meta.json")))
    idx = {a: i for i, a in enumerate(meta["accs"])}
    out = {}
    for sp in sorted(set(accs_keep.values())):
        sub = [a for a in accs_keep if accs_keep[a] == sp and a in idx]
        if len(sub) < 3:
            continue
        ii = [idx[a] for a in sub]
        d = D[np.ix_(ii, ii)]
        d = (d + d.T) / 2
        np.fill_diagonal(d, 0)
        Z = linkage(squareform(d, checks=False), method="average")
        lab = fcluster(Z, cut, criterion="distance")
        for a, l in zip(sub, lab):
            out[a] = f"{sp}{l}"
    return out


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, "genome_matrix.csv"))))
    fams = [f for f in rows[0] if f.startswith("F") and len(f) == 4]
    amr = {r["assembly"]: r for r in csv.DictReader(open(os.path.join(OUT, "amr_per_genome.csv")))}
    amr_fams = [k for k in next(iter(amr.values())) if k != "assembly"]
    for r in rows:
        a = amr.get(r["assembly"], {})
        for k in amr_fams:
            r["AMR:" + k] = 1 if float(a.get(k, 0) or 0) >= AMR_CALL else 0
    # AMRFinderPlus 對同一個基因有多筆參考序列（mecA 與 mecA-ceftar、qacA 與 qacB
    # 在這個面板上逐株完全一致），兩筆都報等於把同一件事數兩次。
    # aph(2'')-Ia 在葡萄球菌就是雙功能酶 aac(6')-Ie-aph(2'')-Ia 的另一半，
    # 逐株完全共現（181/181），兩個都列等於同一個基因報兩次。
    DUPLICATE_REFS = {"mecA-ceftar", "qacB", "aph(2'')-Ia"}
    keep_amr = [k for k in amr_fams
                if k not in DUPLICATE_REFS
                and 20 <= sum(r["AMR:" + k] for r in rows) <= len(rows) - 20]

    # 逐株完全共現的決定子成對記下來：表裡兩列數字一模一樣不是複製貼上的錯，
    # 是這兩個基因在本面板上分不開，要在註腳講清楚。
    twins = []
    for i, k1 in enumerate(amr_fams):
        v1 = [r["AMR:" + k1] for r in rows]
        if sum(v1) < 5:
            continue
        for k2 in amr_fams[i + 1:]:
            v2 = [r["AMR:" + k2] for r in rows]
            if v1 == v2:
                twins.append([k1, k2, sum(v1)])
    json.dump({"identical_determinant_pairs": twins,
               "duplicate_reference_entries_dropped": sorted(DUPLICATE_REFS)},
              open(os.path.join(OUT, "determinant_twins.json"), "w"), indent=1)
    print("identical determinant pairs:", twins, flush=True)
    sp_of = {r["assembly"]: r["species"] for r in rows}
    lins = {c: lineages(sp_of, c) for c in LINEAGE_CUTS} \
        if os.path.exists(os.path.join(OUT, "tree_D.npy")) else {c: {} for c in LINEAGE_CUTS}
    print(f"{len(rows)} genomes; {len(fams)} families; "
          f"{len(keep_amr)} resistance determinants with >=20 carriers; "
          + "; ".join(f"cut {c}: {len(set(l.values()))} lineages"
                      for c, l in lins.items()), flush=True)

    # 物種名單：有 >=20 株的物種才進檢定（12 種裡的 8 種）
    SPECIES = sorted({r["species"] for r in rows}
                     & {sp for sp in {r["species"] for r in rows}
                        if sum(1 for r in rows if r["species"] == sp) >= 20})
    print("species tested:", SPECIES, flush=True)

    res = []
    # 1. family x species（每個物種 one-vs-rest）
    for f in fams:
        for sp in SPECIES:
            a = sum(1 for r in rows if r["species"] == sp and r[f])
            b = sum(1 for r in rows if r["species"] == sp and not r[f])
            c = sum(1 for r in rows if r["species"] != sp and r[f])
            d = sum(1 for r in rows if r["species"] != sp and not r[f])
            if a + c == 0:
                continue
            orr, p = fisher(a, b, c, d)
            res.append(dict(test="family~species", x=f, y=sp + " vs rest",
                            n11=a, n10=b, n01=c, n00=d, odds=orr, p=p))

    # 2. family x family, within species
    for sp in SPECIES:
        sub = [r for r in rows if r["species"] == sp]
        present = [f for f in fams if sum(1 for r in sub if r[f]) >= 10]
        for i, f1 in enumerate(present):
            for f2 in present[i + 1:]:
                a = sum(1 for r in sub if r[f1] and r[f2])
                b = sum(1 for r in sub if r[f1] and not r[f2])
                c = sum(1 for r in sub if not r[f1] and r[f2])
                d = sum(1 for r in sub if not r[f1] and not r[f2])
                orr, p = fisher(a, b, c, d)
                res.append(dict(test=f"family~family({sp})", x=f1, y=f2,
                                n11=a, n10=b, n01=c, n00=d, odds=orr, p=p))

    # 3. family x resistance determinant, within species
    strat_rows = []
    for sp in SPECIES:
        sub = [r for r in rows if r["species"] == sp]
        present = [f for f in fams if sum(1 for r in sub if r[f]) >= 20]
        for f in present:
            for k in keep_amr:
                key = "AMR:" + k
                a = sum(1 for r in sub if r[f] and r[key])
                b = sum(1 for r in sub if r[f] and not r[key])
                c = sum(1 for r in sub if not r[f] and r[key])
                d = sum(1 for r in sub if not r[f] and not r[key])
                if a + c < 10:
                    continue
                orr, p = fisher(a, b, c, d)
                res.append(dict(test=f"family~AMR({sp})", x=f, y=k,
                                n11=a, n10=b, n01=c, n00=d, odds=orr, p=p))
                # lineage-stratified, at three clustering cut-offs
                cmhres = {}
                for cut, lin in lins.items():
                    tabs, mixed, n_mixed = [], 0, 0
                    by_l = collections.defaultdict(list)
                    for r in sub:
                        by_l[lin.get(r["assembly"], "NA")].append(r)
                    for l, rr in by_l.items():
                        t = (sum(1 for r in rr if r[f] and r[key]),
                             sum(1 for r in rr if r[f] and not r[key]),
                             sum(1 for r in rr if not r[f] and r[key]),
                             sum(1 for r in rr if not r[f] and not r[key]))
                        if (t[0] + t[1]) and (t[2] + t[3]):
                            mixed += 1; n_mixed += len(rr)
                        tabs.append(t)
                    cmhres[cut] = cmh(tabs) + (mixed, n_mixed) + cmh_ci(tabs)
                mor, mp, nstr, mixed, n_mixed, mlo, mhi = cmhres[0.020]
                # leave-one-BioProject-out
                worst_p, worst_bp = p, ""
                # 排序：worst_p 打平時，誰被記成 worst_project 不可以取決於
                # set 的走訪順序，否則同一份輸入會長出兩種結果（2026-10-29
                # 乾淨目錄重跑時實際撞到：同樣 p=0.1205，兩個 BioProject 互換）。
                bps = sorted({r["bioproject"] for r in sub})
                for bp in bps:
                    s2 = [r for r in sub if r["bioproject"] != bp]
                    aa = sum(1 for r in s2 if r[f] and r[key])
                    bb = sum(1 for r in s2 if r[f] and not r[key])
                    cc = sum(1 for r in s2 if not r[f] and r[key])
                    dd = sum(1 for r in s2 if not r[f] and not r[key])
                    _, pp = fisher(aa, bb, cc, dd)
                    if pp > worst_p:
                        worst_p, worst_bp = pp, bp
                # 2×2 的四格一定要寫出來：沒有分母的 odds ratio 沒辦法評估。
                # 有零格時 OR 會是 0 或 inf，改報 Haldane-Anscombe（四格各加 0.5）
                # 的點估計，並把「不可估計」旗標打開。
                zero = 0 in (a, b, c, d)
                ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                # 訊號是不是只來自單一 BioProject：剔掉最不利的那個計畫之後就不顯著
                single = worst_p > 0.05 and p <= 0.05
                clo, chi = woolf_ci(a, b, c, d)
                row = dict(species=sp, family=f, determinant=k,
                           odds_lo=round(clo, 3), odds_hi=round(chi, 3),
                           n_fam_det=a, n_fam_nodet=b, n_nofam_det=c, n_nofam_nodet=d,
                           n=a + b + c + d,
                           odds=orr, odds_haldane_anscombe=round(ha, 3),
                           not_estimable=int(zero), p=p,
                           single_project_signal=int(single),
                           worst_p_leave_one_project_out=worst_p,
                           worst_project=worst_bp)
                for cut, (mo, mpv, ns, mx, nmx, lo_, hi_) in cmhres.items():
                    row[f"cmh_odds@{cut}"] = mo
                    row[f"cmh_lo@{cut}"] = lo_
                    row[f"cmh_hi@{cut}"] = hi_
                    row[f"cmh_p@{cut}"] = mpv
                    row[f"informative_strata@{cut}"] = mx
                    row[f"n_in_informative_strata@{cut}"] = nmx
                strat_rows.append(row)

    for test in sorted({r["test"] for r in res}):
        idx = [i for i, r in enumerate(res) if r["test"] == test]
        q = bh([res[i]["p"] for i in idx])
        for i, qq in zip(idx, q):
            res[i]["q"] = qq
    with open(os.path.join(OUT, "associations.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0]))
        w.writeheader(); w.writerows(res)
    if strat_rows:
        with open(os.path.join(OUT, "amr_stratified.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(strat_rows[0]))
            w.writeheader(); w.writerows(strat_rows)
    print("wrote out/associations.csv", len(res), "tests", flush=True)
    sig = [r for r in res if r["q"] < 0.05]
    print("significant after BH:", len(sig))
    for r in sorted(sig, key=lambda r: r["p"])[:15]:
        print(f"  {r['test']:24s} {r['x']:6s} {r['y']:18s} "
              f"OR={r['odds']:.2f} p={r['p']:.2e} q={r['q']:.2e} "
              f"({r['n11']}/{r['n11']+r['n10']} vs {r['n01']}/{r['n01']+r['n00']})")


if __name__ == "__main__":
    main()
