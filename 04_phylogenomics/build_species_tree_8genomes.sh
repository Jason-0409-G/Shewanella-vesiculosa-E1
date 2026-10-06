#!/usr/bin/env bash
# Maximum-likelihood species tree of the eight focal genomes from a random sample of
# single-copy orthogroups (Fig. 2b; Supplementary Note 2).
#
# Steps: sample N orthogroups with a fixed seed -> rename headers from locus tag to
# strain key -> MAFFT --auto per orthogroup -> trimAl -automated1 -> concatenate ->
# IQ-TREE (ModelFinder Plus, 1000 UFBoot, 1000 SH-aLRT).
#
# Input   SCO_DIR  OrthoFinder Single_Copy_Orthologue_Sequences/*.fa (2,417 orthogroups)
# Output  ${OUT_DIR}/Sv8_supermatrix.treefile
# Usage   bash build_species_tree_8genomes.sh [N_OG]      (default 100)
#
# Strain keys: L5 is the internal name of S. vesiculosa E1 (locus-tag prefix PB002_),
# M7 = S. vesiculosa M7 (KDH10_); the others follow prokka_annotate_8genomes.sh.
set -euo pipefail

N_OG=${1:-100}
SEED=42
THREADS=8

SCO_DIR="${SCO_DIR:-./orthofinder_8strains/Single_Copy_Orthologue_Sequences}"
OUT_DIR="${OUT_DIR:-./species_tree}"
export SCO_DIR OUT_DIR N_OG SEED

ALN_DIR="${OUT_DIR}/alignments"
TRIM_DIR="${OUT_DIR}/trimmed"
RENAMED_DIR="${OUT_DIR}/renamed"
SAMPLE_LIST="${OUT_DIR}/sampled_OGs.txt"
mkdir -p "${ALN_DIR}" "${TRIM_DIR}" "${RENAMED_DIR}"

# 1. Sample N orthogroups
ls "${SCO_DIR}" | grep '\.fa$' | sort \
  | python3 -c "
import os, sys, random
items = sys.stdin.read().split()
random.seed(int(os.environ['SEED']))
random.shuffle(items)
for x in sorted(items[:int(os.environ['N_OG'])]):
    print(x)
" > "${SAMPLE_LIST}"

# 2. Rename FASTA headers (locus tag -> strain key)
python3 - <<'PYEOF'
import os
from pathlib import Path

PREFIX_MAP = {
    "PB002_": "L5",
    "KDH10_": "M7",
    "FRIG_":  "frig",
    "LIVI_":  "livi",
    "POLA_":  "pola",
    "PSYC_":  "psyc",
    "SP02_":  "sp02",
    "SP14_":  "sp14",
}

OUT_DIR     = Path(os.environ["OUT_DIR"])
SCO_DIR     = Path(os.environ["SCO_DIR"])
RENAMED_DIR = OUT_DIR / "renamed"

sampled = [l.strip() for l in open(OUT_DIR / "sampled_OGs.txt") if l.strip()]

for og in sampled:
    seen = set()
    with open(SCO_DIR / og) as fi, open(RENAMED_DIR / og, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                tag = line[1:].strip().split()[0]
                key = next(v for p, v in PREFIX_MAP.items() if tag.startswith(p))
                if key in seen:
                    raise RuntimeError(f"{og}: duplicate strain '{key}'")
                seen.add(key)
                fo.write(f">{key}\n")
            else:
                fo.write(line)
PYEOF

# 3. MAFFT alignment per orthogroup
while read -r og; do
    mafft --auto --quiet --thread "${THREADS}" "${RENAMED_DIR}/${og}" \
        > "${ALN_DIR}/${og%.fa}.aln.faa"
done < "${SAMPLE_LIST}"

# 4. trimAl column trimming
while read -r og; do
    trimal -in "${ALN_DIR}/${og%.fa}.aln.faa" \
           -out "${TRIM_DIR}/${og%.fa}.trim.faa" \
           -automated1 -fasta
done < "${SAMPLE_LIST}"

# 5. Concatenate into the supermatrix
python3 - <<'PYEOF'
import os
from pathlib import Path

OUT_DIR  = Path(os.environ["OUT_DIR"])
TRIM_DIR = OUT_DIR / "trimmed"
STRAINS  = ["L5", "M7", "frig", "sp02", "sp14", "livi", "psyc", "pola"]

def read_fa(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

sampled = [l.strip() for l in open(OUT_DIR / "sampled_OGs.txt") if l.strip()]

concat = {s: [] for s in STRAINS}
n_used = n_skip = length = 0
for og in sampled:
    p = TRIM_DIR / f"{og.replace('.fa', '')}.trim.faa"
    seqs = read_fa(p) if p.exists() else {}
    L = len(next(iter(seqs.values()), ""))
    if set(seqs) != set(STRAINS) or L == 0 or any(len(s) != L for s in seqs.values()):
        n_skip += 1
        continue
    for s in STRAINS:
        concat[s].append(seqs[s])
    length += L
    n_used += 1

with open(OUT_DIR / "supermatrix.faa", "w") as f:
    for s in STRAINS:
        f.write(f">{s}\n{''.join(concat[s])}\n")

print(f"used {n_used} orthogroups ({n_skip} skipped); supermatrix length {length:,} aa")
PYEOF

# 6. IQ-TREE
cd "${OUT_DIR}"
iqtree \
    -s supermatrix.faa \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -ntmax "${THREADS}" \
    --seed "${SEED}" \
    --redo \
    -pre Sv8_supermatrix
