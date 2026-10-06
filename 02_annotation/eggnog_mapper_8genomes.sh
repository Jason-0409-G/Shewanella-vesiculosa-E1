#!/usr/bin/env bash
# eggNOG-mapper functional annotation of the eight Shewanella proteomes.
#
# Command line as recorded in line 3 of every *.emapper.annotations header
# (emapper-2.1.13-8d129f4, run 2026-04-25). The eight runs differ only in -i, -o
# and --output_dir.
#
# Input   Prokka protein FASTA for each genome (see prokka_annotate_8genomes.sh)
# Output  {out_dir}/{label}.emapper.annotations, consumed by
#         03_comparative_genomics/scan_functional_modules_8genomes.py,
#         03_comparative_genomics/scan_vesiculosa_shared_genes.py,
#         06_cazy_kegg_transporters/count_transporter_families.py and mine_symbiosis.py.
#         The S. vesiculosa E1 subset of this output is Supplementary Data 3.
#
# The eggNOG database version is not recorded in the output headers; see
# SOFTWARE_VERSIONS.md.
set -euo pipefail

PROKKA_ROOT="${PROKKA_ROOT:-./8strains_prokka_v1156}"

#   prokka sub-directory   file prefix   output label   output dir
STRAINS=(
  "L5_PB002_prokka                     L5            L5_v1156       L5_PB002_eggnog"
  "M7_REF_prokka                       M7            M7_v1156       M7_REF_eggnog"
  "S_frigidimarina_prokka              frig          frig_v1156     S_frigidimarina_eggnog"
  "S_livingstonensis_prokka            livi          livi_v1156     S_livingstonensis_eggnog"
  "S_polaris_prokka                    pola          pola_v1156     S_polaris_eggnog"
  "S_psychromarinicola_prokka          psyc          psyc_v1156     S_psychromarinicola_eggnog"
  "S_sp002836315_prokka                sp02          sp02_v1156     S_sp002836315_eggnog"
  "S_sp014164505_prokka                sp14          sp14_v1156     S_sp014164505_eggnog"
)

emapper.py --version

for entry in "${STRAINS[@]}"; do
    read -r subdir prefix label outdir <<< "${entry}"
    faa="${PROKKA_ROOT}/${subdir}/${prefix}.faa"
    mkdir -p "${outdir}"
    emapper.py \
        -i "${faa}" \
        -o "${label}" \
        --output_dir "${outdir}" \
        --cpu 8 \
        --itype proteins \
        --sensmode default \
        --dmnd_iterate no \
        --override
done
