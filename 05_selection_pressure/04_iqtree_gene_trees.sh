#!/bin/bash
# 04_iqtree_gene_trees.sh -- single-gene maximum-likelihood protein trees (IQ-TREE)
#
# The trees supply the topology for the codeml branch and branch-site models.
# Reference: Nguyen et al. 2015 MBE 32:268-274 (IQ-TREE);
#            Kalyaanamoorthy et al. 2017 Nat Methods 14:587-589 (ModelFinder)
#
# Run on the final four island genes. OG0001992_GH3_NagZ (PB002_02917), part of the
# earlier gene set, lies outside the island and is not included.

set -euo pipefail

PAML_DIR=04_selection_pressure/paml/results

OGS=(
    OG0002997_GH3_e108     # PB002_01800, 5' end of island
    OG0001451_GH1          # PB002_01809, bglB
    OG0000506_MFS          # PB002_01810, MFS transporter
    OG0001452_GH3_e227     # PB002_01811, 3' end of island
)

echo "=== IQ-TREE version ==="
iqtree --version 2>&1 | head -1

for og in "${OGS[@]}"; do
    indir="$PAML_DIR/$og"
    prot_aln="$indir/${og}_prot.aln"

    if [ ! -f "$prot_aln" ]; then
        echo "ERROR: $prot_aln not found, skipping"
        continue
    fi

    echo ""
    echo "=== Building tree for $og ==="

    cd "$indir"

    # -m TEST     ModelFinder selects the substitution model
    # -B 1000     ultrafast bootstrap, 1000 replicates
    # -alrt 1000  SH-aLRT, 1000 replicates
    # -T 4        4 threads
    # --seed 42   fixed seed (weakly supported nodes otherwise vary between runs)
    iqtree \
        -s "${og}_prot.aln" \
        -m TEST \
        -B 1000 \
        -alrt 1000 \
        --prefix "${og}_tree" \
        -T 4 \
        --seed 42 \
        --redo \
        > iqtree.log 2>&1

    # Outputs: {og}_tree.treefile (Newick with support values), {og}_tree.iqtree (report)
    if [ -f "${og}_tree.treefile" ]; then
        grep "Best-fit model" "${og}_tree.iqtree" | head -1 || true
    else
        echo "  Tree building failed for $og, see iqtree.log"
    fi

    cd - > /dev/null
done

echo ""
echo "Done. Next: 05_label_E1_foreground_branch.py"
