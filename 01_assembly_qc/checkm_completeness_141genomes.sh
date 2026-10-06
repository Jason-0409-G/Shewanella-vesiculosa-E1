#!/usr/bin/env bash
# Genome completeness of the 141-genome panel with CheckM v1 (lineage_wf --genes).
# The values populate the Completeness column of Supplementary Data 2 and the
# completeness track of Fig. 2a.
#
# Input   REF_FAA_DIR  one protein FASTA per reference genome (*.faa), 140 files,
#                      each named by its tree leaf label (Supplementary Data 2,
#                      "Tree leaf label"). The genomes are the NCBI assemblies listed
#                      in Supplementary Data 2; CDS were predicted with Prodigal.
#         E1_FAA       protein FASTA of S. vesiculosa E1 (Prokka output, see
#                      ../02_annotation/prokka_annotate_8genomes.sh; genome under
#                      BioProject PRJNA1478518)
# Output  ${OUT_DIR}/checkm_summary.tsv. These completeness values are published as the
#         Completeness column of Supplementary Data 2, which 09_figures/Fig2a_itol_annotations.py
#         reads to build the iTOL tracks for Fig. 2a.
#
# The E1 leaf label in the tree files is PB002_S_vesiculosa_L5 (internal strain name).
set -euo pipefail

REF_FAA_DIR="${REF_FAA_DIR:-./proteins}"
E1_FAA="${E1_FAA:-./8strains_prokka_v1156/L5_PB002_prokka/L5.faa}"
OUT_DIR="${OUT_DIR:-./checkm_141strains}"
THREADS="${THREADS:-8}"
E1_LEAF="PB002_S_vesiculosa_L5"

# CheckM expects one sub-directory per genome
BINS_DIR="${OUT_DIR}/bins"
mkdir -p "${BINS_DIR}"

for faa in "${REF_FAA_DIR}"/*.faa; do
    strain=$(basename "${faa}" .faa)
    mkdir -p "${BINS_DIR}/${strain}"
    cp "${faa}" "${BINS_DIR}/${strain}/${strain}.faa"
done

mkdir -p "${BINS_DIR}/${E1_LEAF}"
cp "${E1_FAA}" "${BINS_DIR}/${E1_LEAF}/${E1_LEAF}.faa"

CHECKM_OUT="${OUT_DIR}/lineage_wf"
mkdir -p "${CHECKM_OUT}"

checkm lineage_wf \
    --genes \
    -x faa \
    -t "${THREADS}" \
    --pplacer_threads 4 \
    "${BINS_DIR}" \
    "${CHECKM_OUT}"

checkm qa \
    "${CHECKM_OUT}/lineage.ms" \
    "${CHECKM_OUT}" \
    -o 2 \
    --tab_table \
    -f "${OUT_DIR}/checkm_summary.tsv"
