#!/usr/bin/env bash
# CAZyme annotation of the eight focal proteomes with dbCAN3 (run_dbcan v5.2.8).
#
# Input   PROKKA_ROOT  Prokka output of prokka_annotate_8genomes.sh (uses {prefix}.faa)
#         DBCAN_DB     dbCAN database v5-2_9-13-2025, from `run_dbcan database --db_dir ...`
# Output  ${OUTPUT_ROOT}/{label}_dbcan/overview.tsv and the per-tool result files.
#         The two-of-three consensus (HMMER, DIAMOND, dbCAN_sub) is applied by
#         06_cazy_kegg_transporters/build_cazy_consensus_matrix.py.
set -euo pipefail

PROKKA_ROOT="${PROKKA_ROOT:-./8strains_prokka_v1156}"
OUTPUT_ROOT="${OUTPUT_ROOT:-./8strains_dbcan_v528}"
DBCAN_DB="${DBCAN_DB:-./dbcan_v528_db}"
CPUS="${CPUS:-8}"

mkdir -p "${OUTPUT_ROOT}"
LOG_FILE="${OUTPUT_ROOT}/_batch_log.txt"
run_dbcan version 2>&1 | tail -1 | tee "${LOG_FILE}"

#   label | prokka directory | file prefix
STRAINS=(
    "L5_PB002|L5_PB002_prokka|L5"
    "M7_REF|M7_REF_prokka|M7"
    "S_frigidimarina|S_frigidimarina_prokka|frig"
    "S_livingstonensis|S_livingstonensis_prokka|livi"
    "S_polaris|S_polaris_prokka|pola"
    "S_psychromarinicola|S_psychromarinicola_prokka|psyc"
    "S_sp002836315|S_sp002836315_prokka|sp02"
    "S_sp014164505|S_sp014164505_prokka|sp14"
)

for entry in "${STRAINS[@]}"; do
    IFS="|" read -r label prokka_dir prefix <<< "${entry}"
    echo "${label}" | tee -a "${LOG_FILE}"
    mkdir -p "${OUTPUT_ROOT}/${label}_dbcan"

    run_dbcan CAZyme_annotation \
        --input_raw_data "${PROKKA_ROOT}/${prokka_dir}/${prefix}.faa" \
        --output_dir "${OUTPUT_ROOT}/${label}_dbcan" \
        --db_dir "${DBCAN_DB}" \
        --mode protein \
        --threads "${CPUS}" \
        --log-level INFO >> "${LOG_FILE}" 2>&1
done
