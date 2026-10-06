#!/usr/bin/env bash
# Assign the island transporter PB002_01808 to BglT rather than YicJ.
#
# Prokka annotates PB002_01808 as "yicJ" (inner membrane symporter YicJ) from the
# E. coli protein. Pairwise blastp against the two candidate assignments shows the
# E1 protein is far closer to the Shewanella baltica beta-glucoside transporter
# BglT than to E. coli YicJ, which is why the manuscript calls it bglT.
#
# Reference proteins (downloaded from UniProt):
#   A3D1N7  S. baltica OS155 BglT, locus Sbal_1130, GenBank ABN60650
#   P31435  E. coli K-12 YicJ
#
# Result of the original run, reproduced with BLAST+ default parameters. The manuscript
# quotes these rounded, as 88 % and 45 %:
#   PB002_01808 vs A3D1N7   88.2 % identity over 439 aa, 99 % query coverage
#   PB002_01808 vs P31435   44.7 % identity over 441 aa, 99 % query coverage
#
# Input   L5.faa   Prokka proteins of S. vesiculosa E1
# Output  bglT_assignment.tsv   one row per reference
set -euo pipefail

FAA="${FAA:-L5.faa}"
LOCUS="${LOCUS:-PB002_01808}"
WORK="${WORK:-./bglT_assignment}"
mkdir -p "${WORK}"

blastp -version | head -1

python3 - "${FAA}" "${LOCUS}" "${WORK}/query.faa" <<'PY'
import sys
from Bio import SeqIO
faa, locus, out = sys.argv[1:4]
for rec in SeqIO.parse(faa, "fasta"):
    if rec.id == locus:
        with open(out, "w") as fh:
            fh.write(f">{rec.id}\n{rec.seq}\n")
        print(f"{locus}: {len(rec.seq)} aa")
        break
else:
    sys.exit(f"{locus} not found in {faa}")
PY

for acc in A3D1N7 P31435; do
    curl -sS -f "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${WORK}/${acc}.faa"
done
cat "${WORK}"/A3D1N7.faa "${WORK}"/P31435.faa > "${WORK}/refs.faa"

printf 'query\tsubject\tpident\tlength\tqcovs\tevalue\tbitscore\n' > "${WORK}/bglT_assignment.tsv"
blastp -query "${WORK}/query.faa" -subject "${WORK}/refs.faa" \
       -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
       >> "${WORK}/bglT_assignment.tsv"

column -t -s$'\t' "${WORK}/bglT_assignment.tsv"
