#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""家族編號 → 論文裡用的名字，以及分離部位的歸類。

名字只照 RefSeq 產物註記與功能域組成取，**不宣稱化學結構**：
lugdunin 與 aureusimine 是註記裡就有的名字，其餘用「做什麼的機器」描述。
改家族編號（重跑 s11）之後一定要回來同步這個檔。
"""

FAMILY = {
    "F001": ("aureusimine NRPS (epidermidis)", "NRPS/PKS"),
    "F002": ("aureusimine NRPS (capitis)", "NRPS/PKS"),
    "F003": ("aryl-capped NRPS (lugdunensis)", "NRPS/PKS"),
    "F004": ("gallidermin-family lantibiotic (capitis)", "RiPP"),
    "F005": ("aureusimine NRPS x YcaO RiPP hybrid", "hybrid"),
    "F006": ("lugdunin NRPS", "NRPS/PKS"),
    "F007": ("C39 exporter locus (haemolyticus)", "RiPP"),
    "F008": ("gallidermin-family lantibiotic, second type", "RiPP"),
    "F009": ("lacticin 481-family lantibiotic (plasmid)", "RiPP"),
    "F010": ("bacteriocin + C39 exporter (plasmid)", "RiPP"),
    "F011": ("class III lanthipeptide, LanKC (plasmid)", "RiPP"),
    "F012": ("aryl-capped NRPS (caprae)", "NRPS/PKS"),
    "F013": ("aureusimine NRPS, divergent locus", "NRPS/PKS"),
}

SHORT = {
    "F001": "aureusimine (epi)", "F002": "aureusimine (cap)",
    "F003": "aryl-NRPS (lug)", "F004": "gallidermin-family",
    "F005": "NRPS x YcaO hybrid", "F006": "lugdunin",
    "F007": "C39 exporter", "F008": "gallidermin-family II",
    "F009": "lacticin 481-family", "F010": "bacteriocin + C39",
    "F011": "class III lanthipeptide", "F012": "aryl-NRPS (cap)",
    "F013": "aureusimine (divergent)",
}

# 分離部位歸類（只用於描述面板，不做部位關聯的主結論）
SOURCE_GROUPS = [
    ("blood", ("blood", "bacteremia", "bacteraemia", "blood culture",
               "peripheral blood", "platelet", "red blood cell")),
    ("device / implant", ("driveline", "hip", "knee", "prosthe", "implant",
                          "catheter", "shunt", "pacemaker", "graft",
                          "synovial", "joint", "device", "lvad")),
    ("skin / nose", ("skin", "forehead", "axilla", "nose", "nasal", "nares",
                     "scalp", "back", "arm", "toe", "foot", "ear", "cheek")),
    ("urine", ("urine", "urinary", "bladder")),
    ("other clinical", ("wound", "pus", "abscess", "bone", "osteo", "csf",
                        "sputum", "lung", "pleur", "peritone", "eye",
                        "conjunctiv", "tissue", "swab", "clinical", "hospital",
                        "infection", "milk", "breast")),
]


def source_group(s):
    s = (s or "").lower()
    if not s or s in ("missing", "not applicable", "not collected", "na",
                      "unknown", "none"):
        return "not stated"
    for label, keys in SOURCE_GROUPS:
        if any(k in s for k in keys):
            return label
    return "other / environmental"


def fam_name(f):
    return FAMILY.get(f, (f, ""))[0]


def fam_short(f):
    return SHORT.get(f, f)
