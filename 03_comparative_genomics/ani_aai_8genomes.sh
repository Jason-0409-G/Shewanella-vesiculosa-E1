#!/usr/bin/env bash
# Pairwise ANI and AAI among the eight Shewanella genomes (Fig. 2c).
#
# FastANI parameters are those recorded in the log of the original run
# (k = 16, fragment length 3,000, 8 threads, all defaults for v1.34); the command-line
# form is a reconstruction consistent with the inputs and outputs that survive.
# EzAAI was run with its default identity and coverage thresholds (0.4 / 0.5), as
# recorded in the header of AAI_8x8.tsv; the thread count and search backend used
# at the time were not recorded.
#
# Input   genomes/{E1,M7,frig,livi,pola,psyc,sp02,sp14}.fna   for ANI
#         proteins/{...}.faa                                  Prokka proteins, for AAI
# Output  ANI_8x8.tsv(.matrix), AAI_8x8.tsv(.mtx), read by 09_figures/Fig2c_ani_aai_heatmap.R
#
# E1 (L5) vs M7: ANI 98.5 %, AAI 99.0 %. fastANI is not symmetric, so the two
# directions differ slightly; the heatmap script reads the long-format table and
# takes one direction per triangle.
set -euo pipefail

GENOMES="${GENOMES:-./genomes}"
PROTEINS="${PROTEINS:-./proteins}"
# Internal strain keys, as used in the stored results and read by
# 09_figures/Fig2c_ani_aai_heatmap.R. L5 is the internal name of E1.
LABELS=(L5 M7 frig livi pola psyc sp02 sp14)

# -- ANI ----------------------------------------------------------------------
ls "${GENOMES}"/*.fna > genome_list.txt
fastANI --ql genome_list.txt --rl genome_list.txt \
        -o ANI_8x8.tsv -t 8 --matrix

# -- AAI ----------------------------------------------------------------------
mkdir -p ezaai_db
for lab in "${LABELS[@]}"; do
    ezaai convert -i "${PROTEINS}/${lab}.faa" -s prot \
                  -o "ezaai_db/${lab}.db" -l "${lab}"
done
ezaai calculate -i ezaai_db -j ezaai_db -o AAI_8x8.tsv -mtx AAI_8x8.mtx
