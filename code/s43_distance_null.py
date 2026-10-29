#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""基因簇與抗藥基因的距離：跟隨機擺放比，還是比較近嗎？（2026-09-20 審閱意見）

「質體上中位 29 kb、染色體上 182 kb」本身不是證據——質體短，隨便擺也會很近。
這裡在**同一條 replicon 上**把抗藥基因的位置隨機重擺 1,000 次（保留 replicon 長度、
抗藥基因數目與基因簇位置），得到距離的虛無分布，再問觀察到的距離是不是顯著偏小。

輸出 out/distance_null.json。
"""
import collections, csv, json, os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "out")
B = 1000
rng = np.random.default_rng(20260920)


def circular_gap(pos, lo, hi, L):
    """基因簇 [lo, hi] 與點 pos 在環狀 replicon 上的最短距離。"""
    if lo <= pos <= hi:
        return 0
    d1 = (lo - pos) % L
    d2 = (pos - hi) % L
    return int(min(d1, d2))


def main():
    loc = list(csv.DictReader(open(os.path.join(OUT, "colocalisation.csv"))))
    amr = list(csv.DictReader(open(os.path.join(OUT, "amr_located.csv"))))
    cl = {r["cluster"]: r for r in csv.DictReader(
        open(os.path.join(OUT, "clusters.tsv")), delimiter="\t")}
    by_rep = collections.defaultdict(list)
    for r in amr:
        if r["start"] and r["end"]:
            by_rep[(r["assembly"], r["replicon"])].append(
                (int(r["start"]), int(r["end"]), r["determinant"]))

    obs = {"chromosome": [], "plasmid": []}
    null = {"chromosome": [], "plasmid": []}
    per_cluster = []
    for r in loc:
        if r["nearest_bp"] == "":
            continue
        kind = "plasmid" if r["replicon_type"] != "chromosome" else "chromosome"
        c = cl[r["cluster"]]
        L = int(c["replicon_len"])
        lo, hi = int(c["start"]), int(c["end"])
        genes = by_rep.get((r["assembly"], r["replicon"]), [])
        genes = [g for g in genes if g[2] not in ("fusD", "mecA-ceftar", "qacB")]
        if not genes:
            continue
        d_obs = int(r["nearest_bp"])
        # 隨機重擺：每個抗藥基因在 replicon 上均勻取一個起點，長度不變
        lens = np.array([g[1] - g[0] for g in genes])
        draws = rng.integers(0, L, size=(B, len(genes)))
        d_null = np.empty(B)
        for b in range(B):
            d_null[b] = min(circular_gap(int(p + l // 2), lo, hi, L)
                            for p, l in zip(draws[b], lens))
        obs[kind].append(d_obs)
        null[kind].append(float(np.median(d_null)))
        per_cluster.append(dict(cluster=r["cluster"], family=r["family"],
                                replicon_type=kind, replicon_len=L,
                                n_amr=len(genes), observed_bp=d_obs,
                                null_median_bp=float(np.median(d_null)),
                                p_one_sided=float((np.sum(d_null <= d_obs) + 1) / (B + 1))))

    res = {"n_permutations": B, "groups": {}, "family_level": {}}
    for kind in ("chromosome", "plasmid"):
        if not obs[kind]:
            continue
        o = np.array(obs[kind], dtype=float)
        n = np.array(null[kind], dtype=float)
        # 逐個基因簇配對比較：觀察值 vs 自己那條 replicon 的虛無中位數
        ratio = o / np.maximum(n, 1)
        from scipy import stats
        try:
            w = stats.wilcoxon(o, n)
            p = float(w.pvalue)
        except Exception:
            p = float("nan")
        res["groups"][kind] = {
            "n_clusters": len(o),
            "observed_median_kb": float(np.median(o) / 1000),
            "null_median_kb": float(np.median(n) / 1000),
            "median_ratio_observed_over_null": float(np.median(ratio)),
            "wilcoxon_p": p,
            "share_closer_than_null": float(np.mean(o < n))}
        g = res["groups"][kind]
        print(f"  {kind:11s} n={g['n_clusters']:4d}  observed {g['observed_median_kb']:7.1f} kb"
              f"   null {g['null_median_kb']:7.1f} kb   ratio {g['median_ratio_observed_over_null']:.2f}"
              f"   closer than null in {100 * g['share_closer_than_null']:.0f}%"
              f"   Wilcoxon P = {g['wilcoxon_p']:.3g}", flush=True)

    # 以「基因簇」為單位做檢定是假重複：染色體那組 396 個座裡有 288 個是同一個
    # F001 座在 288 株 S. epidermidis 裡的複本。改成每個家族一個中位數再配對檢定，
    # 有效樣本數是家族數而不是座數（2026-09-20 審閱意見）。
    from scipy import stats as _st
    byfam = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in per_cluster:
        byfam[r["replicon_type"]][r["family"]].append((r["observed_bp"], r["null_median_bp"]))
    for kind, fams in byfam.items():
        o = np.array([np.median([x[0] for x in v]) for v in fams.values()], dtype=float)
        n = np.array([np.median([x[1] for x in v]) for v in fams.values()], dtype=float)
        if len(o) < 3:
            res["family_level"][kind] = {"n_families": len(o), "wilcoxon_p": None}
            continue
        try:
            p_ = float(_st.wilcoxon(o, n).pvalue)
        except Exception:
            p_ = float("nan")
        # 距離才是要報的量（審稿人 1），所以中位數配 bootstrap 區間；
        # 重抽的單位是家族，和上面的檢定一致。
        brng = np.random.default_rng(20261029)   # 名字不可叫 rng：會遮住排列用的那個
        idx = brng.integers(0, len(o), size=(10000, len(o)))
        bo = np.median(o[idx], axis=1) / 1000
        bn = np.median(n[idx], axis=1) / 1000
        bd = np.median((o - n)[idx], axis=1) / 1000
        res["family_level"][kind] = {
            "n_families": len(o),
            "observed_median_kb": float(np.median(o) / 1000),
            "observed_ci_kb": [float(np.percentile(bo, 2.5)),
                               float(np.percentile(bo, 97.5))],
            "null_median_kb": float(np.median(n) / 1000),
            "null_ci_kb": [float(np.percentile(bn, 2.5)),
                           float(np.percentile(bn, 97.5))],
            "median_difference_kb": float(np.median(o - n) / 1000),
            "difference_ci_kb": [float(np.percentile(bd, 2.5)),
                                 float(np.percentile(bd, 97.5))],
            "closer_than_null": int(np.sum(o < n)),
            "wilcoxon_p": p_}
        g = res["family_level"][kind]
        print(f"  family-level {kind:11s} {g['n_families']:3d} families  "
              f"observed {g['observed_median_kb']:7.1f} kb  null {g['null_median_kb']:7.1f} kb  "
              f"closer in {g['closer_than_null']}/{g['n_families']}  P = {g['wilcoxon_p']:.3g}",
              flush=True)

    json.dump(res, open(os.path.join(OUT, "distance_null.json"), "w"), indent=1)
    with open(os.path.join(OUT, "distance_null.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_cluster[0]))
        w.writeheader(); w.writerows(per_cluster)
    print("-> out/distance_null.json, out/distance_null.csv", flush=True)


if __name__ == "__main__":
    main()
