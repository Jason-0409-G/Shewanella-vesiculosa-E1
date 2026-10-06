#!/usr/bin/env python3
"""Fig. 5a - annual carbon flux of terminal beta-1,4-glucan processing (scenario violin + sensitivity).

Left: dot-and-range chart on a log axis, one row each for krill carbon intake,
cellobiose carbon hydrolysed by gut bacteria, and its partition into growth and
respiration. Scenario ensemble: Q/PP, phi and epsilon are drawn log-uniformly and
AE and GGE uniformly, independently, within the ranges of Supplementary Data 7
(fixed seed). The ensemble only visualises the range and is not a probability
distribution. Symbols are the nominal values (or the epsilon = 1 upper bound),
thick bars the epsilon = 0.1-1 range at nominal parameters, thin whiskers the
combined envelope of the ranges. Inset: one-at-a-time sensitivity as fold change
from the nominal value.

Usage: python Fig5a_violin_sensitivity.py [N_VIOLIN] [PALETTE]

The published panel uses the defaults, N_VIOLIN=1 (violin on the cellobiose-carbon
row only) and PALETTE=B_route (colours matched to panel b). The other options
are earlier drafts that were not used in the manuscript: N_VIOLIN=4 draws a
violin on every row; the palettes "gold", "A_orig" and "C_blue" are alternative
colour schemes.

Output (current directory): Fig5a_violin_sensitivity_{N_VIOLIN}violin[s]_{PALETTE}.{svg,pdf,png}
Canvas 12.6 x 6.6 in, Arial 12 pt, editable SVG text.
"""
import sys
N_VIOLIN = int(sys.argv[1]) if len(sys.argv) > 1 else 1
PAL = sys.argv[2] if len(sys.argv) > 2 else "B_route"
import matplotlib as mpl
import matplotlib.pyplot as plt

# text, muted text, grid, slate; colours sampled from panel b
DK = "#1B2E4D"; MUT = "#566675"; GRID = "#E3E8F0"; SLATE = "#566675"
# row colours: intake, cellobiose (fill, stroke), growth, respiration, band
PALETTES = {
    "gold":  ("#3D67AD", "#DEB053", "#936830", "#518851", "#566675", "#FBF3E1"),
    "A_orig": ("#486CA8", "#B03050", "#B03050", "#549054", "#486CA8", "#FBF1F4"),
    "B_route": ("#3D67AD", "#AE1844", "#AE1844", "#518851", "#2E7A76", "#FBEFF2"),
    "C_blue": ("#3D67AD", "#AE1844", "#AE1844", "#6691C3", "#9DB6D8", "#FBEFF2"),
}
C_IN, C_CB, C_CB_D, C_GR, C_RE, C_BAND = PALETTES[PAL]
mpl.rcParams.update({"font.family": ["Arial"], "font.size": 12,
                     "axes.linewidth": 1.2, "pdf.fonttype": 42, "svg.fonttype": "none",
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "text.color": DK, "axes.edgecolor": DK, "axes.labelcolor": DK,
                     "xtick.color": DK, "ytick.color": DK,
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})
FS_SUB, FS_VAL, FS_SMALL = 10.5, 11.5, 10

fig = plt.figure(figsize=(12.6, 6.6))

# ---- left: scenario-ensemble violins on a log axis ----
import numpy as np
rng = np.random.default_rng(20260929)
NS = 200_000
lu = lambda lo, hi: np.exp(rng.uniform(np.log(lo), np.log(hi), NS))
NPP = 1949.0
QPP, PHI, EPS = lu(0.09, 0.30), lu(0.002, 0.015), lu(0.1, 1.0)
AE, GGE = rng.uniform(0.75, 0.85, NS), rng.uniform(0.20, 0.30, NS)
K2 = GGE / AE
intake = NPP * QPP
hyd = intake * PHI * EPS
grow, resp = hyd * K2, hyd * (1 - K2)
nominal = {"intake": NPP * 0.18, "hyd": NPP * 0.18 * 0.005 * 1.0}
nominal["grow"], nominal["resp"] = nominal["hyd"] * 0.3125, nominal["hyd"] * 0.6875
env = {"intake": (NPP * 0.09, NPP * 0.30),
       "hyd": (NPP * 0.09 * 0.002 * 0.1, NPP * 0.30 * 0.015 * 1.0)}
env["grow"] = (env["hyd"][0] * 0.20 / 0.85, env["hyd"][1] * 0.30 / 0.75)
env["resp"] = (env["hyd"][0] * (1 - 0.30 / 0.75), env["hyd"][1] * (1 - 0.20 / 0.85))

ax = fig.add_axes([0.225, 0.15, 0.745, 0.80])
rows = [  # y, label, sub-label, samples, key, colour
    (4.0, "Krill carbon intake", "× Q/PP", intake, "intake", C_IN),
    (2.8, "Cellobiose carbon hydrolysed", "by gut bacteria (e.g., E1)\n× φ × ε", hyd, "hyd", C_CB),
    (1.6, "→ Growth", "× GGE/AE", grow, "grow", C_GR),
    (0.6, "→ Respiration", "× (1 − GGE/AE)", resp, "resp", C_RE),
]
ax.plot(NPP, 5.2, "o", ms=10, mfc="white", mec=SLATE, mew=2.0, zorder=3)
ax.text(NPP / 1.45, 5.2, "1,949", fontsize=FS_VAL, va="center", ha="right", color=DK)
ax.text(-0.03, 5.22, "Southern Ocean NPP", transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=12, fontweight="bold")
grid = np.linspace(np.log10(0.005), np.log10(4000), 800)
for y, lab, sub, s, key, col in rows:
    ls = np.log10(s)
    bw = 0.06
    hist, edges = np.histogram(ls, bins=400, range=(grid[0], grid[-1]), density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    kern = np.exp(-0.5 * ((grid[:, None] - centers[None, :]) / bw) ** 2)
    dens = kern @ hist; dens /= dens.max()
    half = 0.36 * dens
    xs = 10 ** grid
    mask = dens > 0.01
    line = C_CB_D if key == "hyd" else col
    if N_VIOLIN == 4 or key == "hyd":
        ax.fill_between(xs[mask], y - half[mask], y + half[mask], color=col, alpha=(0.40 if PAL == "gold" else 0.28) if key == "hyd" else 0.28, lw=0, zorder=2)
        ax.plot(xs[mask], y + half[mask], color=line, lw=1.6, zorder=2); ax.plot(xs[mask], y - half[mask], color=line, lw=1.6, zorder=2)
    lo, hi = env[key]
    ax.plot([lo, hi], [y, y], color=line, lw=1.6, zorder=1)
    for c in (lo, hi): ax.plot([c, c], [y - 0.07, y + 0.07], color=line, lw=1.6)
    v = nominal[key]
    if key != "intake":          # thick bar: epsilon 0.1-1 at nominal phi; symbol: epsilon = 1
        ax.plot([v * 0.1, v], [y, y], color=line, lw=5.0, solid_capstyle="butt", zorder=3)
    ax.plot(v, y, "D" if key == "hyd" else "o", ms=9.5, mfc=col, mec=C_CB_D if (key == "hyd" and PAL == "gold") else "white", mew=1.8, zorder=4)
    tag = f"{v:,.0f}" if key == "intake" else f"{v * 0.1:.2f}–{v:.2f}"
    ax.text(max(hi, xs[mask].max()) * 1.3, y, tag, fontsize=FS_VAL, va="center", ha="left",
            color=C_CB_D if key == "hyd" else DK, fontweight="bold" if key == "hyd" else "normal")
    ax.text(-0.03, y + 0.03, lab, transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=12, fontweight="bold")
    ax.text(-0.03, y - 0.05, sub, transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=FS_SUB, color=MUT)
ax.axhspan(2.8 - 0.42, 2.8 + 0.42, color=C_BAND, zorder=0, lw=0)
ax.set_xscale("log"); ax.set_xlim(0.005, 4000); ax.set_ylim(0.1, 5.6); ax.set_yticks([])
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
for d in (0.01, 0.1, 1, 10, 100, 1000): ax.axvline(d, color=GRID, lw=0.6, zorder=0)
ax.set_xticks([0.01, 0.1, 1, 10, 100, 1000]); ax.set_xticklabels(["0.01", "0.1", "1", "10", "100", "1,000"])
ax.set_xlabel("Carbon (Mt C yr$^{-1}$, log scale)", fontsize=12)

# ---- inset: one-at-a-time sensitivity (fold change from the nominal value) ----
ins = fig.add_axes([0.40, 0.67, 0.20, 0.23])
sens = [  # label, low fold, high fold, colour, affected row
    ("Q/PP 9–30%", 0.5, 30 / 18, C_IN),
    ("φ 0.2–1.5%", 0.4, 3.0, C_CB),
    ("ε 0.1–1", 0.1, 1.0, C_CB),
    ("AE 0.75–0.85 (growth)", (0.25 / 0.85) / 0.3125, (0.25 / 0.75) / 0.3125, C_GR),
    ("GGE 0.20–0.30 (growth)", 0.8, 1.2, C_GR),
]
for i, (lab, lo, hi, col) in enumerate(sens[::-1]):
    ins.barh(i, hi - lo, left=lo, height=0.55, color=col, alpha=0.85 if PAL == "gold" else 0.75, lw=1.0 if (col == C_CB and PAL == "gold") else 0, edgecolor=C_CB_D if (col == C_CB and PAL == "gold") else "none")
    ins.text(0.075, i, lab, ha="right", va="center", fontsize=10.5)
ins.axvline(1, color=DK, lw=1.2)
ins.set_xscale("log"); ins.set_xlim(0.08, 5); ins.set_ylim(-0.6, len(sens) - 0.4)
ins.set_yticks([]); ins.set_xticks([0.1, 0.3, 1, 3]); ins.set_xticklabels(["×0.1", "×0.3", "×1", "×3"], fontsize=10.5)
ins.xaxis.set_minor_locator(mpl.ticker.NullLocator())
for sp in ("top", "right", "left"): ins.spines[sp].set_visible(False)
ins.set_title("One-at-a-time sensitivity (fold change from nominal)", fontsize=11.5, fontweight="bold",
              loc="left", x=-0.85, y=1.04)


for ext in ("svg", "pdf", "png"):
    fig.savefig(f"Fig5a_violin_sensitivity_{N_VIOLIN}violin{'s' if N_VIOLIN > 1 else ''}_{PAL}.{ext}", dpi=300 if ext == "png" else None)
print("written: Fig5a_violin_sensitivity_*.{svg,pdf,png}")
