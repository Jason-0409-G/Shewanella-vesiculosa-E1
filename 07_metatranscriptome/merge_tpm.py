#!/usr/bin/env python3
# Identical to the copy written out by run_metatranscriptome_tpm.sh; keep the two in sync.
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
