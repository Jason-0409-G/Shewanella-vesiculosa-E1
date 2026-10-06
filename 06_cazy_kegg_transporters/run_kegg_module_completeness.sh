#!/usr/bin/env bash
# run_kegg_module_completeness.sh
#
# Purpose
#   Run kegg-pathways-completeness (v1.4.3, command `give_completeness`) on the
#   KO annotations of the nine proteomes used in the manuscript (krill host plus
#   eight Shewanella genomes) and write one per-organism module table each.
#
# Inputs
#   ${KO_DIR}/{org}_kegg.ko.txt   two columns, `gene<TAB>KO` (one KO per line),
#                                 for org in host E1 M7 frig livi pola psyc sp002 sp014.
#                                 KO_DIR defaults to ./ko_annotations.
#
# Outputs
#   ${OUT_DIR}/{org}/{org}_pathways.tsv   columns: module_accession, completeness (%),
#                                         pathway_name, pathway_class, matching_ko, missing_ko
#   ${OUT_DIR}/{org}/log.txt              tool log
#   OUT_DIR defaults to ./data/modcomp/out (the path read by Fig4a_kegg_module_completeness.py).
#
# Relation to the paper
#   These per-organism tables are the input of build_module_completeness_table.py,
#   which assembles Supplementary Data 4 and the module matrix behind Fig. 4a.
#
# Notes
#   Built-in graphs.pkl and modules_table.tsv shipped with the package are used
#   (no -g / -t override); the KEGG module definitions are therefore fixed by the
#   installed package version. Other options are left at defaults.

set -euo pipefail

KO_DIR="${KO_DIR:-./ko_annotations}"
OUT_DIR="${OUT_DIR:-./data/modcomp/out}"
ORGS=(host E1 M7 frig livi pola psyc sp002 sp014)

command -v give_completeness >/dev/null || { echo "give_completeness not found" >&2; exit 1; }
echo "kegg-pathways-completeness: $(give_completeness --version)"

for org in "${ORGS[@]}"; do
  in="${KO_DIR}/${org}_kegg.ko.txt"
  [[ -s "$in" ]] || { echo "missing input: $in" >&2; exit 1; }
  mkdir -p "${OUT_DIR}/${org}"
  give_completeness -i "$in" -o "${OUT_DIR}/${org}" -r "${org}" > "${OUT_DIR}/${org}/log.txt" 2>&1
  echo "done: ${org}"
done
