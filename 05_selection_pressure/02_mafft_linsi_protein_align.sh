#!/bin/bash
# 02_mafft_linsi_protein_align.sh -- protein alignment with MAFFT L-INS-i
# (--localpair --maxiterate 1000; Katoh & Standley 2013, MBE 30:772-780)
#
# Output: {og}_prot.aln (aligned FASTA)
#
# OG0001992_GH3_NagZ was part of the original run of this list; the gene lies
# outside the island and was not carried into the manuscript, so it is omitted.

set -euo pipefail

PAML_DIR=04_selection_pressure/paml/results

OGS=(
    OG0001451_GH1
    OG0001452_GH3_e227
    OG0002997_GH3_e108
)

echo "=== MAFFT version ==="
mafft --version 2>&1 | head -1

for og in "${OGS[@]}"; do
    indir="$PAML_DIR/$og"
    prot_in="$indir/${og}_prot.faa"
    prot_aln="$indir/${og}_prot.aln"

    if [ ! -f "$prot_in" ]; then
        echo "ERROR: $prot_in not found, skipping"
        continue
    fi

    n_seqs=$(grep -c "^>" "$prot_in")
    echo ""
    echo "=== Aligning $og ($n_seqs sequences) ==="

    mafft --localpair --maxiterate 1000 \
          --thread 8 \
          "$prot_in" > "$prot_aln" 2> "$indir/mafft.log"

    python3 -c "
seqs = []
cur = []
for line in open('$prot_aln'):
    if line.startswith('>'):
        if cur: seqs.append(''.join(cur))
        cur = []
    else:
        cur.append(line.strip())
if cur: seqs.append(''.join(cur))
lengths = set(len(s) for s in seqs)
print(f'  N seqs: {len(seqs)}, alignment length: {lengths}')
"
    echo "  [saved] $prot_aln"
done

echo ""
echo "Done. Next step: 03_pal2nal_codon_alignment.sh"
