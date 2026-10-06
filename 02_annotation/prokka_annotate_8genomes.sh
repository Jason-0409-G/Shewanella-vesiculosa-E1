#!/usr/bin/env bash
# Uniform Prokka v1.15.6 annotation of the eight focal genomes.
#
# Input   GENOME_DIR  one nucleotide FASTA per genome, named as in the table below.
#                     E1: JBZUHQ000000000 (BioProject PRJNA1478518); the other seven are
#                     the RefSeq assemblies of Supplementary Data 1:
#                       M7  GCF_021560015.1   S. frigidimarina    GCF_000014705.1
#                       S. livingstonensis    GCF_003855395.1     S. polaris          GCF_006385555.1
#                       S. psychromarinicola  GCF_003855155.1     Shewanella sp002836315  GCF_010092445.1
#                       Shewanella sp014164505  GCF_014164505.1
# Output  ${OUTPUT_ROOT}/{label}_prokka/{prefix}.{faa,ffn,gff,gbk,txt,...}
#         The .faa files feed eggnog_mapper_8genomes.sh and dbcan_cazyme_8genomes.sh.
#
# L5 / PB002 are the internal names of strain E1 and are kept in file names and
# locus tags so that they match the intermediate files used by later steps.
set -euo pipefail

GENOME_DIR="${GENOME_DIR:-./genomes_fna}"
OUTPUT_ROOT="${OUTPUT_ROOT:-./8strains_prokka_v1156}"
CPUS="${CPUS:-8}"

mkdir -p "${OUTPUT_ROOT}"
LOG_FILE="${OUTPUT_ROOT}/_batch_log.txt"
prokka --version 2>&1 | tee "${LOG_FILE}"

#   label | species | prefix | locus tag | input file
STRAINS=(
    "L5_PB002|vesiculosa|L5|PB002|PB002_L5.fna"
    "M7_REF|vesiculosa|M7|KDH10|S_vesiculosa_M7.fna"
    "S_frigidimarina|frigidimarina|frig|FRIG|S_frigidimarina.fna"
    "S_livingstonensis|livingstonensis|livi|LIVI|S_livingstonensis.fna"
    "S_polaris|polaris|pola|POLA|S_polaris.fna"
    "S_psychromarinicola|psychromarinicola|psyc|PSYC|S_psychromarinicola.fna"
    "S_sp002836315|sp|sp02|SP02|S_sp002836315.fna"
    "S_sp014164505|sp|sp14|SP14|S_sp014164505.fna"
)

for entry in "${STRAINS[@]}"; do
    IFS="|" read -r label species prefix locustag fna <<< "${entry}"
    echo "${label}" | tee -a "${LOG_FILE}"

    prokka --outdir "${OUTPUT_ROOT}/${label}_prokka" \
           --prefix "${prefix}" \
           --locustag "${locustag}" \
           --kingdom Bacteria \
           --gcode 11 \
           --genus Shewanella \
           --species "${species}" \
           --strain "${prefix}" \
           --mincontiglen 200 \
           --evalue 1e-9 \
           --rfam \
           --addgenes \
           --cpus "${CPUS}" \
           --force \
           "${GENOME_DIR}/${fna}" >> "${LOG_FILE}" 2>&1
done
