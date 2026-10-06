#!/usr/bin/env bash
# Genus-wide maximum-likelihood tree of 141 Shewanellaceae and outgroup genomes (Fig. 2a).
#
# Commands follow the record of the original run (Colab notebook for CheckM and
# IQ-TREE; shewanella_tree.log of 2026-04-13).
#
# The marker alignment is not built from orthogroups: CheckM v1 extracts its own set
# of conserved single-copy marker proteins, aligns them with hmmalign and writes the
# concatenation. The concatenation is already aligned, so the conditional below takes
# the direct-copy branch and MAFFT is not invoked. trimAl is not used for this tree.
# (The eight-genome tree of Fig. 2b is built by build_species_tree_8genomes.sh.)
#
# Input   ${GENOMES}/*.fa   140 reference genomes (NCBI assemblies of Supplementary
#                           Data 2) plus S. vesiculosa E1 (PRJNA1478518), nucleotide FASTA
# Output  ${TREE_DIR}/shewanella_tree.treefile   141 tips
#
# Published run: IQ-TREE 3.1.1, automatically generated seed 163984, model Q.INSECT+R7
# (BIC), lnL -137069.838, alignment 141 sequences x 6,988 columns.
set -euo pipefail

GENOMES="${GENOMES:-./checkm_genomes}"
CHECKM_OUT="${CHECKM_OUT:-./checkm_tree_output}"
TREE_DIR="${TREE_DIR:-./iqtree_output}"
ALIGNED="${ALIGNED:-./phylo_aligned.fasta}"
mkdir -p "${CHECKM_OUT}" "${TREE_DIR}"

# 1. Marker extraction and concatenation (CheckM v1, database checkm_data_2015_01_16)
checkm tree "${GENOMES}" "${CHECKM_OUT}" -x fa -t 4
checkm tree_qa "${CHECKM_OUT}" \
    -f "${CHECKM_OUT}/marker_summary.tsv" --tab_table

# 2. Alignment
CONCAT="${CHECKM_OUT}/storage/tree/concatenated.fasta"
HAS_GAPS=$(head -100 "${CONCAT}" | grep -c '-' || true)
if [ "${HAS_GAPS}" -gt 10 ]; then
    cp "${CONCAT}" "${ALIGNED}"
else
    mafft --auto --thread 4 "${CONCAT}" > "${ALIGNED}"
fi

# 3. Maximum likelihood
iqtree -s "${ALIGNED}" \
       -m MFP \
       -bb 1000 \
       -nt AUTO \
       -pre "${TREE_DIR}/shewanella_tree"
