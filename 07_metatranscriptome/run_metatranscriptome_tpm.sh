#!/usr/bin/env bash
# Metatranscriptome quantification: bwa mem against the E1 gene sequences (L5.ffn), CoverM TPM, TPM matrix.
#
# Input   clean metatranscriptome reads (not host-depleted), {sample}_clean_1.fq.gz / _clean_2.fq.gz
# Output  ${OUTDIR}/04_tpm_matrix/L5_TPM_matrix.tsv and L5_trimmed_mean_matrix.tsv
# Usage   BASEDIR=/path/to/workdir PROJECT=<id> QUEUE=<queue> bash run_metatranscriptome_tpm.sh
#         (jobs are submitted with qsub; adapt the qsub calls to your scheduler)

set -euo pipefail

BASEDIR="${BASEDIR:?set BASEDIR to your working directory}"
CLEAN_DIR="${BASEDIR}/metatrans"                # clean reads (not host-depleted)
SUFFIX_R1="_clean_1.fq.gz"
SUFFIX_R2="_clean_2.fq.gz"
L5_FFN="${BASEDIR}/bins/L5.ffn"
L5_ANNOT="${BASEDIR}/bins/L5.tsv"
OUTDIR="${BASEDIR}/L5_metatrans/metatrans_nohr"
SCRIPTDIR="${BASEDIR}/jobs/metatrans_nohr"
COVERM="coverm"
PYTHON="python3"
THREADS=10
PROJECT="${PROJECT:?set your cluster project id}"
QUEUE="${QUEUE:?set your cluster queue}"

mkdir -p "${OUTDIR}"/{01_index,02_align,03_coverm,04_tpm_matrix,logs}
mkdir -p "${SCRIPTDIR}"

IDX="${OUTDIR}/01_index/L5_contig"

# STEP 1: bwa index (skipped if it already exists)
INDEX_JID=""
if [ ! -f "${IDX}.bwt" ]; then
    cat > "${SCRIPTDIR}/bwa_index.sh" << EOF
#!/usr/bin/env bash
bwa index -p ${IDX} ${L5_FFN} 2> ${OUTDIR}/logs/bwa_index.log
echo "[done] bwa index L5.ffn"
EOF
    INDEX_JID=$(qsub -cwd -l vf=5g,num_proc=1 -q "${QUEUE}" -P "${PROJECT}" \
                     -N "nohr_index" \
                     -o "${OUTDIR}/logs/bwa_index.log" \
                     -e "${OUTDIR}/logs/bwa_index.err" \
                     "${SCRIPTDIR}/bwa_index.sh" | awk '{print $3}')
    echo "  [STEP1] index job: ${INDEX_JID}"
else
    echo "  [STEP1] index exists, skipping"
fi

# STEP 2: alignment and CoverM, one job per sample
SAMPLES=()
for f in "${CLEAN_DIR}"/*"${SUFFIX_R1}"; do
    SAMPLES+=("$(basename "${f}" "${SUFFIX_R1}")")
done
echo "  [STEP2] ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

ALIGN_JOB_IDS=()
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${CLEAN_DIR}/${SAMPLE}${SUFFIX_R1}"
    R2="${CLEAN_DIR}/${SAMPLE}${SUFFIX_R2}"
    BAM="${OUTDIR}/02_align/${SAMPLE}.sorted.bam"
    LOG="${OUTDIR}/logs/${SAMPLE}_bwa.log"
    FLAGSTAT="${OUTDIR}/logs/${SAMPLE}_flagstat.txt"
    COVERM_OUT="${OUTDIR}/03_coverm/${SAMPLE}_coverm.tsv"
    RG="@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA"

    cat > "${SCRIPTDIR}/align_${SAMPLE}.sh" << EOF
#!/usr/bin/env bash
set -euo pipefail

# alignment
if [ ! -f "${BAM}" ]; then
    bwa mem -t ${THREADS} -M -R "${RG}" \\
        ${IDX} ${R1} ${R2} \\
        2> ${LOG} \\
    | samtools sort -@ ${THREADS} -o ${BAM}
    samtools index ${BAM}
    samtools flagstat ${BAM} > ${FLAGSTAT}
    echo "[done align] ${SAMPLE}: \$(grep 'primary mapped' ${FLAGSTAT})"
else
    echo "[skip] ${SAMPLE} BAM exists"
fi

# CoverM (no paired-mode options)
${COVERM} contig \\
    -m trimmed_mean tpm \\
    -b ${BAM} \\
    -t ${THREADS} \\
    --min-covered-fraction 0 \\
    > ${COVERM_OUT}

n=\$(awk 'NR>1 && \$NF+0>0 {c++} END{print c+0}' "${COVERM_OUT}")
echo "[done coverm] ${SAMPLE}: \${n} genes TPM>0"
EOF

    HOLD_OPTS=""
    [ -n "${INDEX_JID}" ] && HOLD_OPTS="-hold_jid ${INDEX_JID}"

    JID=$(qsub -cwd -l vf=20g,num_proc=${THREADS} -q "${QUEUE}" -P "${PROJECT}" \
               -N "nohr_${SAMPLE}" \
               -o "${OUTDIR}/logs/${SAMPLE}_align.log" \
               -e "${OUTDIR}/logs/${SAMPLE}_align.err" \
               ${HOLD_OPTS} \
               "${SCRIPTDIR}/align_${SAMPLE}.sh" | awk '{print $3}')
    ALIGN_JOB_IDS+=("${JID}")
    echo "    ${SAMPLE} job: ${JID}"
done

# STEP 3: merge the per-sample TPM tables
cat > "${SCRIPTDIR}/merge_tpm.py" << 'PYEOF'
"""merge_tpm.py: collect the per-sample CoverM output into TPM and trimmed-mean matrices."""
import argparse
import glob
import os
import sys
from functools import reduce

import pandas as pd


def get_sample(col: str, suffix: str) -> str:
    """Sample name from a CoverM column header (strip the path and the .sorted suffix)."""
    bam_path = col[: col.rfind(f" {suffix}")]
    return os.path.basename(bam_path).replace(".sorted.bam", "").replace(".sorted", "")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverm_dir", required=True)
    parser.add_argument("--outdir",     required=True)
    parser.add_argument("--annot",      default="")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.coverm_dir, "*_coverm.tsv")))
    if not files:
        print(f"[ERROR] no *_coverm.tsv found in {args.coverm_dir}", file=sys.stderr)
        sys.exit(1)

    tpm_dfs = []
    tm_dfs  = []

    for f in files:
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            print(f"[WARN] skipping empty file: {f}", file=sys.stderr)
            continue
        contig_col = df.columns[0]
        # CoverM column headers are "{bam} Trimmed Mean" and "{bam} TPM"
        tm_col  = next(c for c in df.columns if c.endswith(" Trimmed Mean"))
        tpm_col = next(c for c in df.columns if c.endswith(" TPM"))
        sample  = get_sample(tpm_col, "TPM")

        tpm_dfs.append(df[[contig_col, tpm_col]].rename(columns={contig_col: "Gene", tpm_col: sample}))
        tm_dfs.append( df[[contig_col, tm_col ]].rename(columns={contig_col: "Gene", tm_col:  sample}))

    tpm_matrix = reduce(lambda a, b: pd.merge(a, b, on="Gene", how="outer"), tpm_dfs).fillna(0)
    tm_matrix  = reduce(lambda a, b: pd.merge(a, b, on="Gene", how="outer"), tm_dfs ).fillna(0)

    annot_cols = []
    if args.annot and os.path.exists(args.annot):
        annot = pd.read_csv(args.annot, sep="\t")
        annot = annot[annot["ftype"] != "gene"]   # drop empty "gene" rows; keep CDS/rRNA/tRNA
        key   = annot.columns[0]
        tpm_matrix = annot.merge(tpm_matrix, left_on=key, right_on="Gene", how="right").drop(columns=["Gene"])
        tm_matrix  = annot.merge(tm_matrix,  left_on=key, right_on="Gene", how="right").drop(columns=["Gene"])
        annot_cols = list(annot.columns)

    os.makedirs(args.outdir, exist_ok=True)
    tpm_out = os.path.join(args.outdir, "L5_TPM_matrix.tsv")
    tm_out  = os.path.join(args.outdir, "L5_trimmed_mean_matrix.tsv")
    tpm_matrix.to_csv(tpm_out, sep="\t", index=False)
    tm_matrix.to_csv( tm_out,  sep="\t", index=False)

    samples = [c for c in tpm_matrix.columns if c not in annot_cols and c != "Gene"]
    print(f"[done] {tpm_matrix.shape[0]} genes x {len(samples)} samples")
    print(f"  TPM          -> {tpm_out}")
    print(f"  trimmed mean -> {tm_out}")


if __name__ == "__main__":
    main()
PYEOF

MERGE_HOLD=$(IFS=,; echo "${ALIGN_JOB_IDS[*]}")

cat > "${SCRIPTDIR}/merge_tpm.sh" << EOF
#!/usr/bin/env bash
set -euo pipefail
echo "[start merge] \$(date)"
${PYTHON} ${SCRIPTDIR}/merge_tpm.py \\
    --coverm_dir ${OUTDIR}/03_coverm \\
    --outdir     ${OUTDIR}/04_tpm_matrix \\
    --annot      ${L5_ANNOT}
echo "[done merge] \$(date)"
EOF

MERGE_JID=$(qsub -cwd -l vf=8g,num_proc=2 -q "${QUEUE}" -P "${PROJECT}" \
                 -N "nohr_merge" \
                 -o "${OUTDIR}/logs/merge.log" \
                 -e "${OUTDIR}/logs/merge.err" \
                 -hold_jid "${MERGE_HOLD}" \
                 "${SCRIPTDIR}/merge_tpm.sh" | awk '{print $3}')

echo ""
echo "Jobs submitted (reads not host-depleted):"
echo "  index:     ${INDEX_JID:-existing}"
echo "  alignment: ${ALIGN_JOB_IDS[*]}"
echo "  merge TPM: ${MERGE_JID}"
echo "  output:    ${OUTDIR}/04_tpm_matrix/"
