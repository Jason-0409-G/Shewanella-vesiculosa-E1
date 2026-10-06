"""
Prepare PAML-compatible tree files from IQ-TREE output.

For each OG:
  1. read {og}_tree.treefile
  2. strip support values (codeml rejects internal-node labels other than #N marks);
     branch lengths are kept as starting values
  3. write {og}_unlabeled.nwk  (no marks; M0 model)
  4. write {og}_L5fg.nwk       (E1 leaf branch marked #1; branch and branch-site models)

Strain E1 appears under its internal key L5 in file names and sequence IDs.
PAML convention: "leaf_name #1" marks the branch leading to that leaf as foreground.

OG0001992_GH3_NagZ was part of the original run of this list; the gene lies outside
the island and was not carried into the manuscript, so it is omitted.
"""
from __future__ import annotations
import re
from pathlib import Path

PAML_DIR = Path("04_selection_pressure/paml/results")

OGS = [
    "OG0001451_GH1",
    "OG0001452_GH3_e227",
    "OG0002997_GH3_e108",
]


def strip_supports(newick: str) -> str:
    """Remove support values from internal nodes: ')96.7/92:0.02' or ')100:0.01' -> '):...'."""
    cleaned = re.sub(r'\)[\d.]+/[\d.]+:', '):', newick)
    cleaned = re.sub(r'\)[\d.]+:', '):', cleaned)
    return cleaned


def find_e1_leaf(newick: str) -> str:
    """Return the E1 leaf label, e.g. 'L5_PB002_01809'."""
    m = re.search(r'(L5_[A-Z0-9_]+):', newick)
    if not m:
        raise ValueError("No E1 (L5_) leaf found in tree")
    return m.group(1)


def mark_foreground(newick: str, label: str) -> str:
    """'L5_PB002_01809:0.0065' -> 'L5_PB002_01809 #1:0.0065'."""
    return newick.replace(f"{label}:", f"{label} #1:", 1)


def add_paml_header(newick: str, n_taxa: int) -> str:
    """codeml tree files start with '  n_taxa  n_trees'."""
    return f"  {n_taxa}  1\n{newick.strip()}\n"


def main():
    for og in OGS:
        treefile = PAML_DIR / og / f"{og}_tree.treefile"
        if not treefile.exists():
            print(f"{treefile} not found, skipping")
            continue

        print(f"\n=== {og} ===")
        nw_clean = strip_supports(treefile.read_text().strip())
        e1_label = find_e1_leaf(nw_clean)
        n_taxa = nw_clean.count(",") + 1   # leaves = commas + 1 (no internal labels remain)
        print(f"  E1 leaf: {e1_label}; taxa: {n_taxa}")

        unlabeled_out = PAML_DIR / og / f"{og}_unlabeled.nwk"
        unlabeled_out.write_text(add_paml_header(nw_clean, n_taxa))
        print(f"  [saved] {unlabeled_out}")

        fg_out = PAML_DIR / og / f"{og}_L5fg.nwk"
        fg_out.write_text(add_paml_header(mark_foreground(nw_clean, e1_label), n_taxa))
        print(f"  [saved] {fg_out}")

    print("\nDone. Next step: 06_run_codeml_four_models.py")


if __name__ == "__main__":
    main()
