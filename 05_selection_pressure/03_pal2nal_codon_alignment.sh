#!/bin/bash
# 03_pal2nal_codon_alignment.sh -- codon-aware back-translation with PAL2NAL
# (Suyama et al. 2006, Nucleic Acids Res 34:W609-W612)
#
# Input:  protein alignment (.aln) + matching CDS nucleotides (.ffn)
# Output: codon alignment in PAML format (.phy, input to codeml) and FASTA (.fa)
#
# Flags:
#   -nogap          remove all columns containing a gap (and in-frame stop codons),
#                   so the codon alignment is shorter than the protein alignment
#   -output paml    PAML-compatible PHYLIP
#   -codontable 11  bacterial genetic code
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

echo "=== PAL2NAL version ==="
pal2nal.pl --help 2>&1 | head -1

for og in "${OGS[@]}"; do
    indir="$PAML_DIR/$og"
    prot_aln="$indir/${og}_prot.aln"
    nucl_ffn="$indir/${og}_nucl.ffn"
    codon_phy="$indir/${og}_codon.phy"
    codon_fa="$indir/${og}_codon.fa"

    if [ ! -f "$prot_aln" ] || [ ! -f "$nucl_ffn" ]; then
        echo "ERROR: missing files for $og, skipping"
        continue
    fi

    n_prot=$(grep -c "^>" "$prot_aln")
    n_nucl=$(grep -c "^>" "$nucl_ffn")
    echo ""
    echo "=== $og  (prot: $n_prot, nucl: $n_nucl) ==="

    # PAML format, for codeml
    pal2nal.pl "$prot_aln" "$nucl_ffn" \
        -output paml \
        -codontable 11 \
        -nogap \
        > "$codon_phy" 2> "$indir/pal2nal.log" || {
        echo "  pal2nal failed for $og, see $indir/pal2nal.log"
        continue
    }

    # FASTA format (gaps retained)
    pal2nal.pl "$prot_aln" "$nucl_ffn" \
        -output fasta \
        -codontable 11 \
        > "$codon_fa" 2>> "$indir/pal2nal.log" || true

    if [ -s "$codon_phy" ]; then
        echo "  PHYLIP header: $(head -1 "$codon_phy")"
        echo "  [saved] $codon_phy"
        echo "  [saved] $codon_fa"
    else
        echo "  $codon_phy is empty"
    fi
done

echo ""
echo "Done. Next step: 04_iqtree_gene_trees.sh"
