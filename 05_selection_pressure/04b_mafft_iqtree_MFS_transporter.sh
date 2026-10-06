#!/usr/bin/env bash
# MAFFT alignment and IQ-TREE gene tree for the MFS transporter (PB002_01810, OG0000506).
# Run from inside the OG0000506_MFS results directory containing OG0000506_MFS_prot.faa.
#   - MAFFT L-INS-i (Katoh & Standley 2013, MBE 30:772); MAFFT v7.526
#   - IQ-TREE 3.1.1, same options as 04_iqtree_gene_trees.sh
set -euo pipefail

cd "$(dirname "$0")"

PREFIX=OG0000506_MFS
IN_FAA=${PREFIX}_prot.faa
ALN=${PREFIX}_prot.aln
TREE_PREFIX=${PREFIX}_tree

mafft  --version 2>&1 | head -1
iqtree --version 2>&1 | head -1

# L-INS-i expands to: mafft --maxiterate 1000 --localpair
echo "[1] MAFFT L-INS-i"
mafft --maxiterate 1000 --localpair --thread 8 "$IN_FAA" > "$ALN" 2> mafft.log
echo "  -> $ALN ($(grep -c '^>' "$ALN") sequences)"

# -m TEST: ModelFinder; -B 1000: ultrafast bootstrap; -alrt 1000: SH-aLRT; -T 4: threads
echo "[2] IQ-TREE"
iqtree -s "$ALN" \
       -m TEST \
       -B 1000 \
       -alrt 1000 \
       -T 4 \
       --seed 42 \
       --prefix "$TREE_PREFIX" \
       --redo \
       2>&1 | tee iqtree.log

echo ""
echo "treefile: $(pwd)/${TREE_PREFIX}.treefile"
echo "contree:  $(pwd)/${TREE_PREFIX}.contree"
echo "report:   $(pwd)/${TREE_PREFIX}.iqtree"
