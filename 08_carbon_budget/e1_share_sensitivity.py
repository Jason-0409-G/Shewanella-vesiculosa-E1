"""Share of the terminal cellobiose-hydrolysis step attributable to S. vesiculosa E1, and the corresponding
carbon flux (Mt C yr-1), with one-at-a-time sensitivity and a scenario table. Deterministic.

  F_hyd = C_ing x phi x eps,  C_ing = (Q/PP) x NPP                (carbon_load_v32 model)
  p_E1  = s_E1*m / (s_E1*m + f)                                    (E1 share among cellobiose hydrolysers)
  F_E1  = p_E1 x F_hyd

s_E1: relative abundance of E1 in the gut bacterial community; the 1 % is a stipulated upper ceiling on
     abundance (not measured), not E1's share of the terminal hydrolysis step
f  : relative abundance of the other beta-glucosidase-carrying bacteria -- 37.0 %: the five of 12 MAGs
     annotated with K05349/K05350 (Wei et al. 2025 mSystems, Table S5; read relative abundances in text,
     Fig. 3C: 18.6 + 8.9 + 6.5 + 1.9 + 1.1); normalised within MAG-mapped reads, not the whole community
m  : per-cell hydrolysis rate of E1 relative to the other hydrolysers -- not measured; base 1 (abundance-
     weighted apportioning; Stone et al. 2021 qSIP ~1:1; Abu-Ali et al. 2018), range 0.1-10 (Stone et al. 2021,
     deviations "an order of magnitude or more")
"""
import json
from pathlib import Path

NPP = 1949.0
BASE = dict(QPP=0.18, phi=0.005, eps=1.0, s_E1=0.01, f=0.370, m=1.0)   # symbols as in Supplementary S3.2
OAT = [  # (label, parameter, value, basis)
    ("Q/PP low", "QPP", 0.09, "Cavan et al. 2024"),
    ("Q/PP high", "QPP", 0.30, "Cavan et al. 2024"),
    ("phi low", "phi", 0.002, "assumed range"),
    ("phi high", "phi", 0.015, "Schmidt et al. 2006 x Kwok & Wong 2003, ~1.4 % rounded to 1.5 %"),
    ("eps 0.4", "eps", 0.4, "crustacean cellulose digestibility"),
    ("eps 0.1", "eps", 0.1, "gut-passage-time inference"),
    ("s_E1 (E1 abundance) 0.1 %", "s_E1", 0.001, "scenario below the 1 % ceiling"),
    ("m 0.1", "m", 0.1, "Stone et al. 2021"),
    ("m 3", "m", 3.0, "three cellulose-relevant beta-glucosidase genes in E1 (Supplementary S4.5) vs one per carrier MAG"),
    ("m 10", "m", 10.0, "Stone et al. 2021"),
]


def run(QPP, phi, eps, s_E1, f, m):
    F_hyd = NPP * QPP * phi * eps
    p = s_E1 * m / (s_E1 * m + f)
    return {"F_hyd_Mt": F_hyd, "p_E1_pct": 100 * p, "F_E1_t": 1e6 * p * F_hyd}


def rnd(d): return {k: (round(v, 4) if k != "F_E1_t" else round(v, 0)) for k, v in d.items()}


base = run(**BASE)
res = {"base (eps = 1, upper bound for conversion)": rnd(base),
       "base with eps 0.1-1": {"F_hyd_Mt": [round(run(**{**BASE, 'eps': 0.1})["F_hyd_Mt"], 3), round(base["F_hyd_Mt"], 3)],
                               "F_E1_t": [round(run(**{**BASE, 'eps': 0.1})["F_E1_t"]), round(base["F_E1_t"])]},
       "one_at_a_time": {}}
for lab, p, v, basis in OAT:
    out = run(**{**BASE, p: v})
    res["one_at_a_time"][lab] = {**rnd(out), "value": v, "basis": basis,
                                  "F_E1_rel_to_base": round(out["F_E1_t"] / base["F_E1_t"], 3)}
# E1 share and flux as a two-dimensional scenario table (no combined E1 envelope: stipulated values must not be
# multiplied into a single headline number). f = beta-glucosidase carriers among ALL gut bacteria: the 12 MAGs are
# mostly from culture-enriched metagenomes and represent "a small fraction of the gut microbiota" (Wei et al. 2025),
# so f is not constrained; 37 % assumes the unbinned majority carries beta-glucosidases in the same proportion.
# m = 3: gene-number weighting (E1: three cellulose-relevant beta-glucosidase genes, bglB + GH3_e108 + GH3_e227,
# Supplementary S4.5; each carrier MAG: 1 annotated copy, Wei et al. 2025b Table S5), motivated by
# Abu-Ali et al. 2018 (transcripts track genes).
F_nom, F_low = base["F_hyd_Mt"], run(**{**BASE, "eps": 0.1})["F_hyd_Mt"]
table = {}
for f in (0.05, 0.10, 0.20, 0.37):
    for m in (0.1, 1.0, 3.0, 10.0):
        sh = BASE["s_E1"] * m / (BASE["s_E1"] * m + f)
        table[f"f={int(f*100)}% m={m}"] = {"p_E1_pct": round(100 * sh, 2),
                                          "F_E1_t (eps 0.1-1, nominal phi)": [round(1e6 * sh * F_low), round(1e6 * sh * F_nom)]}
res["scenario_table (s_E1 = 1 % stipulated ceiling; f and m not measured)"] = table
Path(__file__).with_suffix(".json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print(json.dumps(res, indent=1, ensure_ascii=False))
