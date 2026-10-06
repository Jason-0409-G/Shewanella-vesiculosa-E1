#!/usr/bin/env python3
"""Fig. 4b - rank-abundance of E1 CDS in the Kin-1 metatranscriptome.

Every CDS is ranked by transcript abundance; the six genes of the
beta-glucoside-processing island are highlighted. Gene symbols and colours
follow the published panel.

Input   L5_TPM_matrix_eggNOG_annotated.tsv  CoverM TPM joined to the Prokka table,
        written by 07_metatranscriptome/annotate_tpm_with_eggnog.py
        (07_metatranscriptome/); optional first argument
Output  Fig4b_island_tpm_rank.svg, Fig4b_island_tpm_rank_editable.svg (clip paths
        removed so that text stays editable in a vector editor), Fig4b_island_tpm_rank.png
"""

import re
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

matplotlib.rcParams.update({
    "font.family":      "Arial",
    "font.size":        7,
    "axes.linewidth":   0.8,
    "xtick.major.width":0.8,
    "ytick.major.width":0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype":     42,
    "svg.fonttype":     "none",   # keep text editable in the vector editor
})

TPM_TABLE = Path(sys.argv[1] if len(sys.argv) > 1 else "L5_TPM_matrix_eggNOG_annotated.tsv")
OUTDIR = Path(".")

# -- read data -------------------------------------------------------------
df = pd.read_csv(TPM_TABLE, sep="\t")
cds = df[df["ftype"] == "CDS"].copy()
cds = cds.dropna(subset=["Kin-1"])

# rank descending: highest TPM is rank 1
cds["rank"] = cds["Kin-1"].rank(ascending=False, method="first").astype(int)
cds = cds.sort_values("rank").reset_index(drop=True)

n_total     = len(cds)
n_expressed = int((cds["Kin-1"] > 1).sum())

# -- genes of the beta-glucoside-processing island -------------------------
# locus -> (display name, colour), sampled from the published panel.
ISLAND = {
    "PB002_01811": ("GH3_e227",        "#B6392B"),   # island-borne beta-glucosidase
    "PB002_01800": ("GH3_e108",        "#266C9A"),
    "PB002_01809": ("bglB",            "#1C8458"),   # GH1 beta-glucosidase
    "PB002_01808": ("bglT",            "#D78C3C"),   # beta-glucoside transporter
    "PB002_01807": ("galR",            "#6A3906"),   # HTH-type transcriptional regulator
    "PB002_01810": ("MFS transporter\n(TPM = 0)", "#763A8E"),
}

island_rows = {}
for locus, (name, col) in ISLAND.items():
    island_rows[locus] = cds[cds["locus_tag"] == locus].iloc[0]

# -- plot ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 3.5))

# background rank curve over all CDS, including those at TPM = 0
# note: TPM = 0 gives log(0+1) = 0, which a log axis would place misleadingly at y = 1,
# so the y axis is linear in log10(TPM+1) and the ticks are labelled with that value
L = lambda t: np.log10(np.asarray(t, dtype=float) + 1.0)
ax.plot(
    cds["rank"].values,
    L(cds["Kin-1"].values),
    color="#B0C4D8", lw=1.5, zorder=1, label="_nolegend_"
)

# divider at TPM = 1
ax.axvline(x=n_expressed, color="#5D6D7E", lw=0.8, ls="--", zorder=2)
ax.text(
    n_expressed - 30, 2.75,
    f"TPM > 1\nn = {n_expressed}",
    ha="right", va="top", fontsize=6, color="#5D6D7E"
)

# span of the island, across the expressed genes only (the MFS transporter is at TPM = 0)
expressed_island = [island_rows[l] for l in
                    ("PB002_01808", "PB002_01811", "PB002_01800",
                     "PB002_01809", "PB002_01807")]
r_left  = min(int(r["rank"]) for r in expressed_island)
r_right = max(int(r["rank"]) for r in expressed_island)
y_arrow = 1.0
ax.annotate(
    "", xy=(r_right, y_arrow), xytext=(r_left, y_arrow),
    arrowprops=dict(arrowstyle="<->", color="#C0392B", lw=1.0)
)
ax.text(
    (r_left + r_right) / 2, y_arrow + 0.15,
    "\u03b2-glucoside-processing island",
    ha="center", va="bottom", fontsize=6, color="#C0392B", fontstyle="italic"
)

# highlighted island genes, with labels
ANNOT_OFFSET = {   # (dx_rank, dy_log); y is the log value, so the offset is additive
    "PB002_01811": (-100,  0.55),
    "PB002_01800": ( 120,  0.45),
    "PB002_01809": ( 100,  0.55),
    "PB002_01808": (-120,  0.45),  # bglT: label up-left, clear of GH3_e227
    "PB002_01807": ( 140,  0.30),  # galR
    "PB002_01810": (-150,  0.60),
}
for locus, (name, col) in ISLAND.items():
    row = island_rows[locus]
    y_plot = L(row["Kin-1"])
    ax.scatter(
        row["rank"], y_plot,
        color=col, s=35, marker="o",
        zorder=5, edgecolors="white", linewidths=0.5
    )
    dx, dy = ANNOT_OFFSET[locus]
    ax.annotate(
        name,
        xy=(row["rank"], y_plot),
        xytext=(row["rank"] + dx, y_plot + dy),
        fontsize=5.5, color=col, fontweight="bold",
        arrowprops=dict(arrowstyle="-", color=col, lw=0.7),
        va="center"
    )

# -- axes --------------------------------------------------------------------
ax.set_xlim(-30, n_total + 80)
ax.set_ylim(-0.1, 3.35)   # log10(TPM+1): 0 is TPM 0, 3 is about TPM 1000; the top CDS runs off scale

ax.set_xlabel("CDS rank (Kin-1 TPM, descending)", fontsize=7.5)
ax.set_ylabel("log$_{10}$(TPM + 1)", fontsize=7.5)

# tick positions are the log values themselves, so bglT reads 2.32 and not 206
ax.set_yticks([0, 1, 2, 3])
ax.set_yticklabels(["0", "1", "2", "3"], fontsize=6.5)

# note that the most abundant CDS lies beyond the axis
ax.text(
    500, 3.22,
    "highest CDS off scale (max ≈ 194,000 TPM)",
    ha="left", va="bottom", fontsize=5.5, color="#888888", fontstyle="italic"
)

ax.tick_params(labelsize=6.5, length=3, width=0.8)

# summary statistics, top right

# -- legend ----------------------------------------------------------------
legend_handles = [Line2D([0], [0], color="#B0C4D8", lw=1.5, label="All CDS")]
for locus in ("PB002_01808", "PB002_01811", "PB002_01800",
              "PB002_01807", "PB002_01809", "PB002_01810"):
    name, col = ISLAND[locus]
    legend_handles.append(
        plt.scatter([], [], marker="o", color=col, s=35,
                    edgecolors="white", lw=0.5, label=name.replace("\n", " "))
    )
ax.legend(
    handles=legend_handles,
    loc="lower left", bbox_to_anchor=(0.01, 0.03),
    fontsize=5.5, framealpha=0.85, edgecolor="#CCCCCC",
    handlelength=1.2, handletextpad=0.4, borderpad=0.5
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout(pad=0.8)

# -- save ------------------------------------------------------------------
out_svg = OUTDIR / "Fig4b_island_tpm_rank.svg"
fig.savefig(out_svg, format="svg", dpi=300, bbox_inches="tight")
print(f"written: {out_svg}")
fig.savefig(OUTDIR / "Fig4b_island_tpm_rank.png", dpi=300, bbox_inches="tight")

# strip clip paths
with open(out_svg, "r", encoding="utf-8") as f:
    content = f.read()
cleaned = re.sub(r'\s*clip-path="url\([^)]+\)"', '', content)
cleaned = re.sub(r'<clipPath[^>]*>.*?</clipPath>', '', cleaned, flags=re.DOTALL)
out_edit = OUTDIR / "Fig4b_island_tpm_rank_editable.svg"
with open(out_edit, "w", encoding="utf-8") as f:
    f.write(cleaned)
print(f"written: {out_edit}")
