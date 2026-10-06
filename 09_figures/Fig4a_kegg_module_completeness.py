#!/usr/bin/env python3
"""Fig. 4a - KEGG module completeness of Antarctic krill, E1 and M7, with E1 expression.

Layout: heatmap of module completeness (3 columns: Krill, E1, M7), a bubble column
with the mean Kin-1 TPM of the E1 genes assigned to each module (bubble size =
log10(TPM + 1)), and a colour bar naming the four module categories (amino acid
synthesis, antioxidation, cofactors and vitamins, central carbon).

Completeness is the KEGG module completeness computed by kegg-pathways-completeness
v1.4.3 (06_cazy_kegg_transporters/run_kegg_module_completeness.sh). Non-zero partial
completeness in krill mostly reflects shared upstream steps or alternative animal
routes, not a genuine capacity for synthesis.

Inputs
------
data/modcomp/out/{host,E1,M7}/{org}_pathways.tsv   per-organism tool output
L5_TPM_matrix_v_L5reasearch.tsv   E1 TPM table with eggNOG KO column (07_metatranscriptome/);
                                  optional first argument

Outputs (current directory)
---------------------------
Fig4a_kegg_module_completeness.{svg,png}
Fig4a_kegg_module_completeness_DATA.tsv   completeness and mean TPM behind each row
"""
from pathlib import Path
import csv, math, re
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.colorbar import ColorbarBase

mpl.rcParams.update({"font.family": "Arial", "font.size": 7,
                     "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42})

WORK = Path("data/modcomp")
TPMF = Path(sys.argv[1] if len(sys.argv) > 1 else "L5_TPM_matrix_v_L5reasearch.tsv")
OUT  = Path(".")

# tool output directory -> column key (L5 is the internal key of E1)
ORG2COL = {"host":"krill","E1":"L5","M7":"M7"}
ORGS = list(ORG2COL.keys())
COLS = ["krill","L5","M7"]
COL_LABELS = {"krill":"Krill","L5":"E1","M7":"M7"}

# Curated modules, in four categories
CAT_ORDER  = ["AA Synthesis","Antioxidation","Cofactors & Vitamins","Central Carbon"]
CAT_COLORS = {"AA Synthesis":"#1565C0","Antioxidation":"#AD1457",
              "Cofactors & Vitamins":"#2E7D32","Central Carbon":"#E65100"}
MODULES = [
 ("AA Synthesis","M00020","Serine"),
 ("AA Synthesis","M00018","Threonine*"),
 ("AA Synthesis","M00021","Cysteine"),
 ("AA Synthesis","M00017","Methionine*"),
 ("AA Synthesis","M00019","Val/Ile*"),
 ("AA Synthesis","M00432","Leucine*"),
 ("AA Synthesis","M00570","Isoleucine*"),
 ("AA Synthesis","M00016","Lysine*"),
 ("AA Synthesis","M00028","Ornithine"),
 ("AA Synthesis","M00844","Arginine"),
 ("AA Synthesis","M00022","Shikimate"),
 ("AA Synthesis","M00023","Tryptophan*"),
 ("AA Synthesis","M00024","Phenylalanine*"),
 ("AA Synthesis","M00025","Tyrosine"),
 ("AA Synthesis","M00026","Histidine*"),
 ("Antioxidation","M00176","Sulfate reduction"),
 ("Antioxidation","M00118","Glutathione"),
 ("Antioxidation","M00121","Heme"),
 ("Antioxidation","M00846","Siroheme"),
 ("Cofactors & Vitamins","M00127","Thiamine (B1)"),
 ("Cofactors & Vitamins","M00125","Riboflavin (B2)"),
 ("Cofactors & Vitamins","M00124","Vit B6 (PLP)"),
 ("Cofactors & Vitamins","M00115","NAD (B3)"),
 ("Cofactors & Vitamins","M00120","CoA (B5)"),
 ("Cofactors & Vitamins","M00123","Biotin (B7)"),
 ("Cofactors & Vitamins","M00126","Folate (B9)"),
 ("Cofactors & Vitamins","M00122","Cobalamin (B12)"),
 ("Cofactors & Vitamins","M00116","Menaquinone (K2)"),
 ("Cofactors & Vitamins","M00117","Ubiquinone"),
 ("Cofactors & Vitamins","M00880","Mo cofactor"),
 ("Central Carbon","M00002","Glycolysis (core)"),
 ("Central Carbon","M00009","TCA cycle"),
 ("Central Carbon","M00004","Pentose-P"),
 ("Central Carbon","M00008","Entner-Doudoroff"),
 ("Central Carbon","M00012","Glyoxylate"),
 ("Central Carbon","M00982","Methylcitrate"),
 ("Central Carbon","M00579","Acetate (Pta-AckA)"),
]

# Per-organism completeness and the KOs matched in E1
comp = {}          # module -> {column: completeness %}
e1_ko = {}         # module -> KOs of E1 in that module
for org in ORGS:
    col = ORG2COL[org]
    with open(WORK/f"out/{org}/{org}_pathways.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = r["module_accession"]
            comp.setdefault(m, {})[col] = float(r["completeness"])
            if org == "E1":
                e1_ko[m] = set(re.findall(r"K\d{5}", r["matching_ko"] or ""))
def C(m,col): return comp.get(m,{}).get(col,0.0)   # missing = 0 %

# Mean Kin-1 TPM of the E1 genes of each module
ko2tpm = {}
with open(TPMF) as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        try: t = float(r["Kin-1"])
        except (ValueError, TypeError): continue
        for k in re.findall(r"K\d{5}", r.get("eggNOG_KO","") or ""):
            ko2tpm.setdefault(k, []).append(t)
def mod_tpm(m):
    ks = e1_ko.get(m, set())
    vals = [v for k in ks for v in ko2tpm.get(k, [])]
    return sum(vals)/len(vals) if vals else 0.0

rows = [(cat,mid,lab,mod_tpm(mid)) for cat,mid,lab in MODULES]
LV  = [math.log10(tp+1.0) for *_,tp in rows]          # bubble size = log10(TPM+1)
lvmax = max(LV) or 1.0
n_row, n_col = len(rows), len(COLS)

# data behind the panel
with open(OUT/"Fig4a_kegg_module_completeness_DATA.tsv","w") as out:
    out.write("category\tmodule\tlabel\t"+"\t".join(COLS)+"\te1_meanTPM_Kin1\n")
    for cat,mid,lab,tp in rows:
        out.write(f"{cat}\t{mid}\t{lab}\t"+"\t".join(f"{C(mid,c):.1f}" for c in COLS)+f"\t{tp:.1f}\n")

# Colours and geometry; ColorBrewer RdBu reversed (blue-white-red)
CMAP = LinearSegmentedColormap.from_list(
    "RdBu_r7", ["#2166AC","#4393C3","#92C5DE","#F7F7F7",
                "#FDDBC7","#F4A582","#D6604D","#B2182B"], N=256)
norm = Normalize(0.0, 100.0)
C_E1, C_EXPR = "#C0504D", "#E07B2E"

ML = 1.35            # left margin for row labels
CW_h = 0.72; W_h = n_col*CW_h
GAP = 0.32; W_b = 0.70
MR = 1.15
RH = 0.155; H = n_row*RH
MT, MB = 0.95, 1.65
fig_w = ML+W_h+GAP+W_b+MR
fig_h = MT+H+MB

fig = plt.figure(figsize=(fig_w, fig_h))
ax_h = fig.add_axes([ML/fig_w, MB/fig_h, W_h/fig_w, H/fig_h])
ax_b = fig.add_axes([(ML+W_h+GAP)/fig_w, MB/fig_h, W_b/fig_w, H/fig_h])
ax_h.set_xlim(-0.5,n_col-0.5); ax_h.set_ylim(-0.5,n_row-0.5); ax_h.invert_yaxis(); ax_h.axis("off")
ax_b.set_xlim(-0.5,0.5);       ax_b.set_ylim(-0.5,n_row-0.5); ax_b.invert_yaxis(); ax_b.axis("off")
j_e1 = COLS.index("L5")

# Completeness heatmap
for i,(cat,mid,lab,tp) in enumerate(rows):
    for j,c in enumerate(COLS):
        ax_h.add_patch(mpatches.Rectangle((j-0.5,i-0.5),1,1,
                       facecolor=CMAP(norm(C(mid,c))), edgecolor="white", lw=0.5, zorder=1))
    # row labels, left of the heatmap
    ax_h.text(-0.85, i, lab, ha="right", va="center", fontsize=5.8, color="#1a1a1a",
              zorder=7, clip_on=False)
for j,c in enumerate(COLS):
    style  = "normal" if c=="krill" else "italic"
    weight = "bold" if c=="L5" else "normal"
    color  = C_E1 if c=="L5" else "#4472C4" if c=="krill" else "#333333"
    ax_h.text(j, n_row-0.5+0.4, COL_LABELS[c], ha="right", va="top", rotation=45,
              fontsize=6.3, fontstyle=style, fontweight=weight, color=color, clip_on=False)
ax_h.add_patch(mpatches.Rectangle((j_e1-0.5,-0.5),1,n_row, fill=False,
               edgecolor=C_E1, lw=1.5, zorder=6, clip_on=False))
ax_h.text((n_col-1)/2, -1.15, "KEGG module completeness  (cell color)",
          ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#333333", clip_on=False)

# E1 expression bubbles
SCALE = 95.0/lvmax
ax_b.scatter([0]*n_row, range(n_row), s=[max(v*SCALE,3.0) for v in LV],
             c=C_EXPR, edgecolors="white", linewidths=0.4, zorder=4)
ax_b.text(0, -0.72, "E1", ha="center", va="bottom", fontsize=6.3,
          fontstyle="italic", fontweight="bold", color=C_E1, clip_on=False)
ax_b.text(0, -1.15, "E1 expression\nlog(TPM+1), Kin-1  (bubble size)",
          ha="center", va="bottom", fontsize=6.5, fontweight="bold", color=C_EXPR, clip_on=False)

# Category colour bar and names
MOD_X = 0.72
cat_bounds = {}
for i,(cat,*_ ) in enumerate(rows):
    cat_bounds.setdefault(cat,[i,i]); cat_bounds[cat][1]=i
for cat in CAT_ORDER:
    r0,r1 = cat_bounds[cat]; col = CAT_COLORS[cat]
    ax_b.add_patch(mpatches.Rectangle((MOD_X, r0-0.5), 0.11, (r1-r0)+1,
                   facecolor=col, alpha=0.9, lw=0, clip_on=False, zorder=5))
    ax_b.text(MOD_X+0.34, (r0+r1)/2, cat, ha="center", va="center", rotation=90,
              fontsize=5.6, weight="bold", color=col, clip_on=False)

# Legends
gx0 = (ML+0.35)/fig_w
cb_ax = fig.add_axes([gx0, 0.050, 1.7/fig_w, 0.012])
cb = ColorbarBase(cb_ax, cmap=CMAP, norm=norm, orientation="horizontal")
cb.set_label("Module completeness (%)  ·  cell color", fontsize=6); cb.ax.tick_params(labelsize=5.5)
cb.set_ticks([0,50,100])

# two size-legend bubbles (lower tertile and maximum of log10(TPM+1))
nz = sorted(v for v in LV if v > 0)
lv_lo = round(nz[len(nz)//3], 1); lv_hi = round(nz[-1], 1)
tvals = sorted({lv_lo, lv_hi})
size_handles = [plt.scatter([],[], s=v*SCALE, facecolor=C_EXPR,
                edgecolor="white", linewidth=0.4, label=f"{v:.1f}") for v in tvals]
fig.legend(handles=size_handles, loc="lower left", bbox_to_anchor=(gx0+1.9/fig_w, 0.015),
           ncol=2, frameon=False, fontsize=5.8,
           title="E1 log(TPM+1) (Kin-1)  ·  bubble size", title_fontsize=6.0,
           handletextpad=0.4, columnspacing=1.4, borderpad=0.2, scatterpoints=1)

stem = "Fig4a_kegg_module_completeness"
fig.savefig(OUT/f"{stem}.svg", bbox_inches="tight", facecolor="white")
fig.savefig(OUT/f"{stem}.png", bbox_inches="tight", facecolor="white", dpi=300)
print(f"written: {OUT/stem}.svg / .png ({fig_w:.2f} x {fig_h:.2f} in), {n_row} rows")
plt.close(fig)
