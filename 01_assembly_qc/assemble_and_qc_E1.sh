#!/usr/bin/env bash
# Read filtering, assembly and quality assessment of the S. vesiculosa E1 chromosome.
#
# Commands follow the record of the original run (Colab notebook for Filtlong and
# CheckM2; Flye log of 2026-04-12, Flye 2.9.6-b1802, for the assembly).
# Reads were delivered as HiFi CCS by the sequencing provider; Filtlong removes only
# a few short or low-quality reads.
#
# Input   RAW_FQ   PacBio HiFi CCS reads, fastq.gz (SRA submission SUB16537385,
#                  BioProject PRJNA1478518; see Data availability in the manuscript)
# Output  ${WORK}/filt.fastq.gz                filtered reads
#         ${WORK}/flye_out/assembly.fasta      one closed circular chromosome, 4,858,980 bp
#         ${WORK}/checkm2_out/quality_report.tsv   100% complete, 0.51% contamination
#
# CheckM2 is run on the assembly split into one file per contig, hence the Bin Id
# "assembly.part_contig_1" in the report.
set -euo pipefail

WORK="${WORK:-./work}"
RAW_FQ="${RAW_FQ:?set RAW_FQ to the PacBio HiFi CCS fastq.gz}"
CHECKM2_DB="${CHECKM2_DB:-./checkm2_db}"
mkdir -p "${WORK}"

# 1. Length and quality filtering
filtlong \
    --min_length 1000 \
    --min_mean_q 90 \
    "${RAW_FQ}" | gzip > "${WORK}/filt.fastq.gz"

# 2. Assembly
flye --pacbio-hifi "${WORK}/filt.fastq.gz" \
     --out-dir "${WORK}/flye_out" \
     --threads 4

# 3. Completeness and contamination
mkdir -p "${WORK}/contigs_split"
seqkit split -i "${WORK}/flye_out/assembly.fasta" \
             -O "${WORK}/contigs_split" --by-id-prefix ""

checkm2 database --download --path "${CHECKM2_DB}"
checkm2 predict --threads 4 \
    --input "${WORK}/contigs_split/" \
    --output-directory "${WORK}/checkm2_out" \
    -x fasta --force
