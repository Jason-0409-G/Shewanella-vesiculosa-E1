"""
Final first-order carbon model (v3, 2026-09-27) for the krill-gut terminal
cellobiose-hydrolysis step.  Deterministic; no probabilistic sampling
(following the best-estimate + one-at-a-time bounds approach of
Belcher et al. 2019 and Cavan et al. 2024).

Chain (per productive season, 181 d, Nov-Apr; Cavan et al. 2024 SI):
  Q        krill food consumption          = (Q/PP) * PP
  C_cell   ingested cellulose carbon       = Q * phi
  C_cb     cellobiose carbon (host-cleaved) = C_cell * eta          (eta = 1: upper bound)
  C_hyd    hydrolysed by gut microbes and
           assimilated by the krill holobiont = C_cb * h             (h = 1: upper bound)
  split of assimilated carbon: growth = K2 * C_hyd, respiration (+excretion) = (1 - K2) * C_hyd,
  K2 = GGE / AE (net growth efficiency)
  E1 share: p_E1 = s_E1*m/(s_E1*m + f); s_E1 = stipulated 1 % ceiling on E1 relative abundance in the gut
  community (not its share of the step), f = other beta-glucosidase carriers, m = relative per-cell activity
"""
import json
from pathlib import Path

# ---- parameters: (best, low, high) -----------------------------------------
PP = 1949.0                    # Mt C yr-1, Arrigo et al. 2008 (denominator of Cavan's Q/PP only)
P = {
    "QPP":  (0.18, 0.09, 0.30),   # Cavan et al. 2024 SI Tables S1-S3, 40-mm krill; indicative only
    "phi":  (0.005, 0.002, 0.015),# assumed (nominal 0.5 %); upper 1.5 %: Schmidt et al. 2006 station dominated by
                                  # armoured dinoflagellates (~60 % of identified volume) x cellulose C 0.6-2.3 % of cell C (Kwok & Wong 2003)
    "eta":  (1.0, 0.1, 1.0),      # upper bound; sensitivity 0.4 (crustacean digestibility) and 0.1-0.35 (retention-time inference)
    "h":    (1.0, 1.0, 1.0),      # upper bound (all host-derived cellobiose hydrolysed by microbes)
    "AE":   (0.80, 0.75, 0.85),   # Clarke et al. 1988 (75/80/85 %); full literature range 0.42-0.94 (Belcher 2019)
    "GGE":  (0.25, 0.20, 0.30),   # Straile 1997 zooplankton range, as used for krill by Belcher 2019 / Cavan 2024
}
ETA_SENS = {"eta 0.4 (crustacean cellulose digestibility 16-53 %)": 0.4,
            "eta 0.35 (retention-time inference, upper)": 0.35,
            "eta 0.1 (retention-time inference, lower)": 0.1}
WINTER = (0.14, 0.20)           # Meyer et al. 2010: winter/autumn feeding relative to late spring
SEASON_D, WINTER_D = 181, 184
N_KRILL = 7.8e14                # Atkinson et al. 2009
A_E1 = 0.01                     # stipulated ceiling on E1 relative abundance in the gut community
F_OTH = 0.37                    # other beta-glucosidase carriers (5 of 12 MAGs, Wei et al. 2025b mSystems); see e1_share_sensitivity.py
R_ACT = 1.0                     # relative per-cell hydrolysis activity (Stone et al. 2021; Abu-Ali et al. 2018)
S_E1 = A_E1 * R_ACT / (A_E1 * R_ACT + F_OTH)   # = 2.63 %
KRILL_WET_TO_C = 0.1075         # Belcher et al. 2019

def run(QPP, phi, eta, h, AE, GGE, winter=0.0):
    Q = PP * QPP * (1 + WINTER_D / SEASON_D * winter)
    C_cell = Q * phi
    C_cb = C_cell * eta
    C_hyd = C_cb * h
    K2 = GGE / AE
    return {"Q": Q, "assimilated_total": Q * AE, "egested_total": Q * (1 - AE),
            "C_cell": C_cell, "C_cb": C_cb, "C_hyd": C_hyd,
            "dAE_pp": C_hyd / Q * 100, "share_of_assimilated_pct": C_hyd / (Q * AE) * 100,
            "K2": K2, "growth": C_hyd * K2, "respiration": C_hyd * (1 - K2),
            "E1_hyd_t": C_hyd * S_E1 * 1e6, "E1_growth_t": C_hyd * S_E1 * K2 * 1e6}

best_kw = {k: v[0] for k, v in P.items()}
best = run(**best_kw)
r = lambda x, n=3: round(x, n)
out = {"description": __doc__.strip().splitlines()[0],
       "parameters": {k: {"best": v[0], "low": v[1], "high": v[2]} for k, v in P.items()},
       "best": {k: r(v, 4) for k, v in best.items()}}
b = best
out["best_units"] = {
    "C_cell_Mt": r(b["C_cell"]),
    "C_cb_Mt_upper": r(b["C_cb"]),
    "growth_as_krill_wet_Mt": r(b["growth"] / KRILL_WET_TO_C, 2),
    "per_krill_ug_C_per_day": r(b["C_cb"] * 1e12 / N_KRILL / SEASON_D * 1e6, 1),
}

# ---- one-at-a-time sensitivity ---------------------------------------------
oat = {}
for k, (bv, lo, hi) in P.items():
    if lo == hi == bv:
        continue
    for tag, val in (("low", lo), ("high", hi)):
        if val == bv:
            continue
        kw = dict(best_kw, **{k: val}); res = run(**kw)
        oat[f"{k} {tag} ({val})"] = {"C_hyd_Mt": r(res["C_hyd"]), "rel_pct": round(res["C_hyd"] / b["C_hyd"] * 100),
                                     "growth_Mt": r(res["growth"]), "respiration_Mt": r(res["respiration"])}
for tag, val in ETA_SENS.items():
    res = run(**dict(best_kw, eta=val))
    oat[tag] = {"C_hyd_Mt": r(res["C_hyd"]), "rel_pct": round(res["C_hyd"] / b["C_hyd"] * 100),
                "growth_Mt": r(res["growth"]), "respiration_Mt": r(res["respiration"])}
for w in WINTER:
    res = run(**best_kw, winter=w)
    oat[f"+ winter feeding ({int(w*100)}% of summer, {WINTER_D} d)"] = {
        "C_hyd_Mt": r(res["C_hyd"]), "rel_pct": round(res["C_hyd"] / b["C_hyd"] * 100),
        "growth_Mt": r(res["growth"]), "respiration_Mt": r(res["respiration"])}
# AE upper value 0.94 is within the measured E. superba range (72-94 %; Clarke et al. 1988, citing Kato et al. 1982).
# The literature minimum AE 0.42 is excluded: Cavan et al. 2024 SI reject low AE as implying implausible
# circumpolar consumption.
res = run(**dict(best_kw, AE=0.94))
oat["AE 0.94 (upper of measured E. superba range)"] = {
    "C_hyd_Mt": r(res["C_hyd"]), "rel_pct": round(res["C_hyd"] / b["C_hyd"] * 100),
    "K2": r(res["K2"]), "growth_Mt": r(res["growth"]), "respiration_Mt": r(res["respiration"])}
out["OAT"] = oat

# ---- combined-extremes envelope (C_hyd) ------------------------------------
lo_kw = dict(best_kw, QPP=P["QPP"][1], phi=P["phi"][1], eta=min(ETA_SENS.values()))
hi_kw = dict(best_kw, QPP=P["QPP"][2], phi=P["phi"][2], eta=1.0)
out["envelope_C_hyd_Mt"] = {"all_low (incl. eta 0.1)": r(run(**lo_kw)["C_hyd"], 4),
                            "all_low (eta fixed 1)": r(run(**dict(lo_kw, eta=1.0))["C_hyd"], 4),
                            "all_high": r(run(**hi_kw)["C_hyd"]),
                            "all_high + winter 20%": r(run(**hi_kw, winter=max(WINTER))["C_hyd"])}

# ---- phi x eta two-way table (C_hyd, Mt C) ---------------------------------
out["phi_x_eta_Mt"] = {f"phi {p*100:.1f}%": {f"eta {e}": r(run(**dict(best_kw, phi=p, eta=e))["C_hyd"])
                                               for e in (0.1, 0.35, 0.4, 1.0)}
                       for p in (0.002, 0.005, 0.010, 0.015)}

# ---- break-even and conversions --------------------------------------------
out["per_0.1pct_phi_Mt"] = r(run(**dict(best_kw, phi=0.001))["C_hyd"])
out["phi_needed_pct"] = {f"{t} Mt C, Q/PP {int(q*100)}%": r(t / (PP * q) * 100, 3)
                         for t in (1.0,) for q in (0.09, 0.18, 0.30)}
# Testable prediction: pellet glucan-C / pellet-C = phi*(1 - eta*h)/(1 - AE); phi unmeasured, so an
# absolute pellet value cannot discriminate hypotheses. Use paired food/pellet ratios (Conover 1966 ratio
# method, adapted): apparent cellulose digestibility = 1 - (G/C)_pellet * (1 - AE) / (G/C)_food.
out["pellet_glucan_pct_of_pellet_C_if_no_digestion"] = {
    f"phi {p*100:.1f}%": [r(p / (1 - ae) * 100, 2) for ae in (0.75, 0.85)] for p in (0.002, 0.005, 0.010, 0.015)}
out["pellet_to_food_glucanC_ratio_expected"] = {"eta*h = 1 (complete)": 0.0,
    "eta*h = 0.4": r((1 - 0.4) / (1 - best_kw["AE"]), 2), "eta*h = 0 (none)": r(1 / (1 - best_kw["AE"]), 2)}
# consistency check: implied egestion rate vs Cavan 2024 median 0.46 mg C ind-1 d-1 (5-95 %: 0.11-1.23)
out["implied_egestion_mgC_per_ind_per_day"] = r(b["egested_total"] * 1e12 / N_KRILL / SEASON_D * 1e3, 3)
out["E1"] = {"s_E1": A_E1, "f": F_OTH, "m": R_ACT, "p_E1_pct": r(S_E1 * 100, 2),   # symbols as in Supplementary Note 3: s_E1 = abundance ceiling, p_E1 = share
             "F_E1_t_nominal_phi_eps_0.1_1": [round(run(**dict(best_kw, eta=0.1))["E1_hyd_t"]), round(b["E1_hyd_t"])],
             "F_E1_t_phi_1.5pct_eps_0.1_1": [round(run(**dict(best_kw, phi=0.015, eta=0.1))["E1_hyd_t"]),
                                             round(run(**dict(best_kw, phi=0.015))["E1_hyd_t"])],
             "E1_growth_t_best": round(b["E1_growth_t"])}

path = Path(__file__).with_name("carbon_load_v32_results.json")
path.write_text(json.dumps(out, indent=2, ensure_ascii=False))
print(json.dumps(out, indent=2, ensure_ascii=False))
