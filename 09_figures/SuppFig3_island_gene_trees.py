"""Fig. S3 (data-driven version) - protein ML trees for the four island genes.

Reads the IQ-TREE .treefile outputs of 05_selection_pressure (A bglB/GH1,
B GH3_e108, C GH3_e227, D MFS transporter), roots on S. livingstonensis,
labels nodes with SH-aLRT/UFBoot when either value reaches SUP_MIN (85, as in the
published panel), draws E1 in red
and one scale bar per tree. Text is kept editable (svg.fonttype = none).
SuppFig3_island_gene_trees_layout.py is the self-contained version with the
newicks inlined and the leaf order of the submitted figure.
Run from 09_figures/; writes
svgout/Fig. S3_from_treefiles.{svg,png}.
"""
from pathlib import Path
from Bio import Phylo
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

matplotlib.rcParams.update({
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
})


SUP_MIN = 85.0   # support threshold, matching the published panel
RES = Path("04_selection_pressure/paml/results")   # where 04_iqtree_gene_trees.sh and 04b write;
                                                   # see 05_selection_pressure/README.md on this path
OUT = Path("svgout")
E1_RED = "#A03A2C"

SP = {"L5": "S. vesiculosa E1", "M7": "S. vesiculosa M7", "frig": "S. frigidimarina",
      "sp02": "S. sp002836315", "sp14": "S. sp014164505", "pola": "S. polaris",
      "psyc": "S. psychromarinicola", "livi": "S. livingstonensis"}
LOCUS_TO_SHORT = {"PB002": "L5", "KDH10": "M7", "FRIG": "frig", "SP02": "sp02",
                  "SP14": "sp14", "POLA": "pola", "PSYC": "psyc", "LIVI": "livi"}
# the 8-taxon GH3_e227 tree uses full species names as leaf labels
FULL_TO_SHORT = {"S_vesiculosa_E1": "L5", "S_vesiculosa_M7": "M7",
                 "S_frigidimarina": "frig", "S_sp002836315": "sp02",
                 "S_sp014164505": "sp14", "S_polaris": "pola",
                 "S_psychromarinicola": "psyc", "S_livingstonensis": "livi"}

TREES = [  # (panel, treefile, title, scale-bar length)
 ("A", RES / "OG0001451_GH1/OG0001451_GH1_tree.treefile",
      "GH1  ($\\it{bglB}$, OG0001451)", 0.1),
 ("B", RES / "OG0002997_GH3_e108/OG0002997_GH3_e108_tree.treefile",
      "GH3_e108  (OG0002997)", 0.1),
 ("C", RES / "OG0001452_GH3_e227_8tax/e227_8tax.treefile",
      "GH3_e227  (OG0001452)", 0.05),
 ("D", RES / "OG0000506_MFS/OG0000506_MFS_tree.treefile",
      "MFS transporter  (OG0000506)", 0.05),
]


def prefix(name):
    if name in FULL_TO_SHORT:
        return FULL_TO_SHORT[name]
    return LOCUS_TO_SHORT.get(name.split("_", 1)[0], name.split("_", 1)[0])


def draw_tree(ax, treefile, title, scale_len, panel):
    t = Phylo.read(str(treefile), "newick")
    livi = [l for l in t.get_terminals() if prefix(l.name) == "livi"]
    if livi:
        t.root_with_outgroup(livi[0])

    def setx(cl, x0):
        cl._x = x0 + (cl.branch_length or 0.0)
        for c in cl.clades:
            setx(c, cl._x)
    setx(t.root, 0.0)
    leaves = t.get_terminals()
    for i, lf in enumerate(leaves):
        lf._y = i

    def sety(cl):
        if cl.is_terminal():
            return cl._y
        cl._y = sum(sety(c) for c in cl.clades) / len(cl.clades)
        return cl._y
    sety(t.root)
    n = len(leaves)
    xmax = max(lf._x for lf in leaves)

    # all branches as one compound path (one <path> per tree in the SVG)
    segs = []
    def collect(cl):
        for c in cl.clades:
            segs.append(((cl._x, cl._y), (cl._x, c._y)))   # vertical connector
            segs.append(((cl._x, c._y), (c._x, c._y)))     # horizontal branch
            collect(c)
    collect(t.root)
    segs.append(((-xmax * 0.02, t.root._y), (t.root._x, t.root._y)))  # root stub
    verts, codes = [], []
    for a, b in segs:
        verts += [a, b]; codes += [MplPath.MOVETO, MplPath.LINETO]
    ax.add_patch(PathPatch(MplPath(verts, codes), edgecolor="#333333", facecolor="none", lw=1.2))
    by = -0.85
    ax.add_patch(PathPatch(MplPath([(0, by), (scale_len, by)], [MplPath.MOVETO, MplPath.LINETO]),
                           edgecolor="black", facecolor="none", lw=1.3))  # scale bar

    for lf in leaves:
        pre = prefix(lf.name); sp = SP.get(pre, pre); is_e1 = (pre == "L5")
        ax.text(lf._x + xmax * 0.02, lf._y, sp, va="center", ha="left", fontsize=8,
                style="italic", color=E1_RED if is_e1 else "#222222",
                fontweight="bold" if is_e1 else "normal", clip_on=False)

    for cl in t.get_nonterminals():
        lab = cl.name or ""
        if "/" in lab:
            try:
                if any(float(p) >= SUP_MIN for p in lab.split("/")):
                    ax.text(cl._x - xmax * 0.006, cl._y + 0.12, lab, fontsize=6,
                            ha="right", va="bottom", color="#777777", clip_on=False)
            except ValueError:
                pass

    ax.text(scale_len / 2, by - 0.38, str(scale_len), ha="center", va="top", fontsize=7)

    ax.set_xlim(-xmax * 0.03, xmax * 2.4)
    ax.set_ylim(-1.6, n - 0.4)
    ax.axis("off")
    ax.set_title(title, fontsize=8.5, loc="left", pad=6, color="#222222")
    ax.text(-0.04, 1.04, panel, transform=ax.transAxes, fontsize=14,
            fontweight="bold", ha="right", va="bottom", color="#111111")


fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
for (panel, tf, title, sc), ax in zip(TREES, axs.flat):
    draw_tree(ax, tf, title, sc, panel)
plt.tight_layout(w_pad=5, h_pad=4)

stem = "Fig. S3_from_treefiles"
fig.savefig(OUT / f"{stem}.svg", format="svg", bbox_inches="tight", facecolor="white")
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
print(f"saved {OUT / stem}.svg / .png")
plt.close(fig)
