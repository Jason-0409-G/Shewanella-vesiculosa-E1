"""Classify OrthoFinder orthogroups across the eight Shewanella genomes.

Categories (E1 is labelled L5 in the OrthoFinder input):
  L5_strict_unique   E1 has >= 1 gene, the other seven genomes have none
  Sv_species_marker  E1 and M7 both have >= 1 gene, the other six have none
  Core_conserved     all eight genomes have >= 1 gene
  Single-copy core   exactly one gene in each of the eight genomes (reported only)
"""
from __future__ import annotations
import pandas as pd
from pathlib import Path

OF_ROOT = Path("orthofinder_8strains/input/OrthoFinder/Results_run_v1")
OUT_DIR = Path("orthofinder_analysis")
OUT_DIR.mkdir(exist_ok=True)

STRAINS = ["L5", "M7", "frig", "livi", "pola", "psyc", "sp02", "sp14"]
SV_STRAINS = ["L5", "M7"]


def categorize_og(df: pd.DataFrame) -> dict:
    cats = {"L5_strict_unique": [], "Sv_species_marker": [], "Core_conserved": []}
    for _, row in df.iterrows():
        og = row["Orthogroup"]
        l5_n, m7_n = row["L5"], row["M7"]

        if l5_n >= 1 and all(row[s] == 0 for s in STRAINS if s != "L5"):
            cats["L5_strict_unique"].append(og)

        if l5_n >= 1 and m7_n >= 1 and all(row[s] == 0 for s in STRAINS if s not in SV_STRAINS):
            cats["Sv_species_marker"].append(og)

        if all(row[s] >= 1 for s in STRAINS):
            cats["Core_conserved"].append(og)
    return cats


def gene_ids(og_genes: pd.DataFrame, ogs: list, strain: str, col: str) -> pd.DataFrame:
    """List the gene IDs of one strain within the given orthogroups."""
    sub = og_genes[og_genes["Orthogroup"].isin(ogs)]
    rows = []
    for _, row in sub.iterrows():
        cell = str(row[strain]) if pd.notna(row[strain]) else ""
        for g in cell.split(", "):
            g = g.strip()
            if g:
                rows.append({"OG": row["Orthogroup"], col: g})
    return pd.DataFrame(rows)


def main():
    counts_df = pd.read_csv(OF_ROOT / "Orthogroups" / "Orthogroups.GeneCount.tsv", sep="\t")
    og_genes = pd.read_csv(OF_ROOT / "Orthogroups" / "Orthogroups.tsv", sep="\t")
    print(f"Orthogroups: {len(counts_df)}")

    cats = categorize_og(counts_df)
    for cat, ogs in cats.items():
        print(f"  {cat:<20} {len(ogs):>5} OGs")
        (OUT_DIR / f"OGs_{cat}.txt").write_text("\n".join(ogs))

    n_single_copy = int((counts_df[STRAINS] == 1).all(axis=1).sum())
    print(f"  {'Single_copy_core':<20} {n_single_copy:>5} OGs")

    unique_df = gene_ids(og_genes, cats["L5_strict_unique"], "L5", "gene_id")
    unique_df.to_csv(OUT_DIR / "L5_strict_unique_genes.tsv", sep="\t", index=False)
    print(f"E1 strictly unique genes: {len(unique_df)}")

    sv_df = gene_ids(og_genes, cats["Sv_species_marker"], "L5", "L5_gene_id")
    sv_df.to_csv(OUT_DIR / "Sv_species_markers_L5_genes.tsv", sep="\t", index=False)
    print(f"E1 genes in S. vesiculosa-restricted OGs: {len(sv_df)}")


if __name__ == "__main__":
    main()
