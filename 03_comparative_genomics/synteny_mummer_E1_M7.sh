#!/usr/bin/env bash
# Whole-chromosome synteny between S. vesiculosa E1 and the reference M7 (Fig. 3a).
#
# Commands recovered from the .delta headers of the original run; re-execution
# reproduced the manuscript output files. Settings that determine the result:
#   * --maxmatch (nucmer default gives 365 alignments instead of 625)
#   * delta-filter thresholds (give the final 127 blocks)
#   * M7 is reverse-complemented first, to match the E1 orientation
#
# Input   E1.fna, M7.fna   chromosomes (Prokka input sequences)
# Output  *.coords   127 blocks = 120 forward + 7 inverted, read by
#                    09_figures/Fig3a_circular_synteny_E1_M7.py
#         *.report   91.78 % of E1 and 93.34 % of M7 aligned, 98.65 % mean identity
#
# dnadiff also reports a 1-to-1 block count (178); the figure and text use the
# delta-filter count of 127.
set -euo pipefail

E1="${E1:-E1.fna}"
M7="${M7:-M7.fna}"
M7RC="${M7RC:-M7_rc.fna}"

nucmer --version
seqkit seq -p -r "${M7}" > "${M7RC}"           # reverse complement

# -- alignment in E1 orientation, for the synteny figure ----------------------
nucmer --maxmatch -p E1_vs_M7rc "${E1}" "${M7RC}"
delta-filter -1 -i 95 -l 1000 E1_vs_M7rc.delta > E1_vs_M7rc.filtered.delta
show-coords -T -l -c E1_vs_M7rc.filtered.delta > E1_vs_M7rc.coords

awk 'NR>4' E1_vs_M7rc.coords | wc -l            # 127 blocks

# -- genome-wide coverage and identity ---------------------------------------
dnadiff -p E1_vs_M7_diff "${E1}" "${M7}"
grep -E "AlignedBases|AvgIdentity" E1_vs_M7_diff.report
