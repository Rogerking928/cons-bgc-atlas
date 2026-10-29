#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fetch Pfam HMM profiles from InterPro (gzipped HMM per accession)."""
import gzip, io, os, sys, time, urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "data", "hmm")
os.makedirs(OUT, exist_ok=True)
URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{}/?annotation=hmm"

PFAMS = {
    # NRPS / PKS catalytic domains (same set as the S. mutans island paper)
    "PF00109": "KS_N", "PF02801": "KS_C", "PF16197": "KS_assoc",
    "PF00501": "A", "PF13193": "A_C",
    "PF00550": "PCP_ACP", "PF00668": "C",
    "PF00698": "AT", "PF08659": "KR", "PF00975": "TE",
    "PF07993": "NAD_binding_4", "PF08242": "MT", "PF00202": "Aminotran_3",
    # RiPP / bacteriocin machinery (葡萄球菌版：羊毛硫肽、sactipeptide、環狀細菌素)
    "PF05147": "LanC", "PF04738": "LanB_N", "PF14028": "LanB_C",
    "PF13575": "LanM_cyc", "PF16934": "Mersacidin", "PF04604": "L_biotic_typeA",
    "PF03412": "Peptidase_C39", "PF02624": "YcaO", "PF05402": "PqqD",
    "PF01721": "Bacteriocin_II", "PF10439": "Bacteriocin_IIc",
    "PF09221": "Bacteriocin_IId", "PF12173": "BacteriocIIc_cy",
    "PF11758": "Bacteriocin_IIi", "PF08129": "Antimicrobial17",
    "PF17508": "MccV", "PF08951": "EntA_Immun", "PF00082": "Peptidase_S8",
    "PF04055": "Radical_SAM",
    # 葡萄球菌特有的前驅肽家族（腸球菌那套沒有）
    "PF02052": "Gallidermin", "PF08130": "Antimicrobial18",
    "PF04369": "Lactococcin_972", "PF19402": "Lanthipeptide_prec",
    # mobility markers
    "PF00589": "Phage_integrase", "PF13495": "Phage_int_SAM_4",
    "PF01609": "DDE_Tnp_1", "PF00665": "rve", "PF13683": "rve_3",
    "PF02914": "Recombinase_serine", "PF07508": "Recombinase",
    "PF03050": "DDE_Tnp_IS66", "PF01548": "DEDD_Tnp_IS110",
}

for pf, name in PFAMS.items():
    dst = os.path.join(OUT, f"{pf}.hmm")
    if os.path.exists(dst) and os.path.getsize(dst) > 100:
        continue
    for attempt in range(4):
        try:
            with urllib.request.urlopen(URL.format(pf), timeout=180) as fh:
                raw = fh.read()
            try:
                raw = gzip.decompress(raw)
            except Exception:
                pass
            open(dst, "wb").write(raw)
            print(pf, name, len(raw), flush=True)
            break
        except Exception as e:
            if attempt == 3:
                print("FAIL", pf, e, flush=True)
            time.sleep(4 * (attempt + 1))
    time.sleep(0.5)

cat = os.path.join(OUT, "library.hmm")
with open(cat, "w") as out:
    for pf in PFAMS:
        p = os.path.join(OUT, f"{pf}.hmm")
        if os.path.exists(p):
            out.write(open(p).read())
print("concatenated ->", cat, os.path.getsize(cat), flush=True)
