#!/usr/bin/env python3
"""Fig. 2b - maximum-likelihood species tree of the eight Shewanella genomes.

Reads the IQ-TREE output for the eight-strain supermatrix (100 single-copy
orthogroups sampled with seed 42; Supplementary Note 2), re-roots it on
S. livingstonensis as published, and draws a rectangular cladogram with branch
lengths to scale and node support labels.

IQ-TREE was run with -m MFP -bb 1000 -alrt 1000, so each internal node carries a
pair "SH-aLRT / ultrafast bootstrap". As in the published panel, only the
bootstrap half is drawn (--support both labels both) and nodes with bootstrap
below 80 are left unlabelled (--min-support 0 labels all; the unlabelled node
is the S. frigidimarina + S. sp002836315 pair, bootstrap 78).

Input
-----
Sv8_supermatrix.treefile   newick with "aLRT/UFboot" node labels

Output
------
Fig2b_species_tree_8genomes.{svg,png}
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

matplotlib.rcParams.update({"svg.fonttype": "none", "font.family": "Arial",
                            "font.size": 8})

# tip label in the tree file -> name used in the manuscript
TIP_NAMES = {
    "L5":   "S. vesiculosa E1",
    "M7":   "S. vesiculosa M7",
    "frig": "S. frigidimarina",
    "livi": "S. livingstonensis",
    "pola": "S. polaris",
    "psyc": "S. psychromarinicola",
    "sp02": "S. sp002836315",
    "sp14": "S. sp014164505",
}
ROOT_AT = "livi"
FOCAL = "S. vesiculosa E1"


def support_text(raw, mode):
    """IQ-TREE writes internal labels as 'aLRT/UFboot'."""
    if not raw:
        return ""
    txt = str(raw)
    if "/" not in txt:
        return txt
    alrt, boot = txt.split("/", 1)
    if mode == "both":
        return f"{alrt}/{boot}"
    if mode == "alrt":
        return alrt
    return boot


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("treefile", nargs="?", default="Sv8_supermatrix.treefile")
    ap.add_argument("--support", choices=("boot", "alrt", "both"), default="boot")
    ap.add_argument("--min-support", type=float, default=80,
                    help="label only nodes whose bootstrap is at least this value")
    args = ap.parse_args()

    tree = Phylo.read(args.treefile, "newick")
    tree.root_with_outgroup({"name": ROOT_AT})
    tree.ladderize(reverse=False)

    for clade in tree.find_clades():
        if clade.name in TIP_NAMES:
            clade.name = TIP_NAMES[clade.name]

    labels = {}
    for clade in tree.get_nonterminals():
        raw = clade.confidence if clade.confidence is not None else clade.name
        txt = support_text(raw, args.support)
        if txt and float(str(raw).split("/")[-1]) >= args.min_support:
            labels[clade] = txt
        clade.name = None            # keep Bio.Phylo from printing it as a tip

    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    Phylo.draw(
        tree, axes=ax, do_show=False,
        branch_labels=lambda c: labels.get(c, ""),
        label_colors=lambda n: "#B6392B" if n == FOCAL else "#222222",
    )
    for txt in ax.texts:
        txt.set_fontsize(7.5)
        if txt.get_text().strip().startswith("S."):
            txt.set_fontstyle("italic")

    ax.set_xlabel("substitutions per site", fontsize=7.5)
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)

    n_tips = len(tree.get_terminals())
    print(f"{args.treefile}: {n_tips} tips, rooted on {TIP_NAMES[ROOT_AT]}")
    print(f"  node labels ({args.support}): "
          + ", ".join(sorted(labels.values(), reverse=True)))

    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"Fig2b_species_tree_8genomes.{ext}",
                    dpi=300 if ext == "png" else None, bbox_inches="tight")
    print("written: Fig2b_species_tree_8genomes.svg / .png")


if __name__ == "__main__":
    main()
