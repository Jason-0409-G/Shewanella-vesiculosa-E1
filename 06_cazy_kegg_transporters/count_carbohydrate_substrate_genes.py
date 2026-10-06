#!/usr/bin/env python3
"""Tally carbohydrate-active genes by substrate group across the eight Shewanella genomes.

Produces the substrate x strain counts reported in Supplementary Data 5, part B.
Input is the dbCAN3 >=2-of-3-tool consensus family matrix written by
build_cazy_consensus_matrix.py; the same matrix underlies the Fig. 3c heatmap.
"""

import pandas as pd
import os

INPUT_FILE = "all_CAZy_family_matrix.tsv"   # written by build_cazy_consensus_matrix.py
OUTPUT_DIR = "carbohydrate_substrate_counts"

SUBSTRATE_GROUPS = {
    "Cellulose_BG":         ["GH1","GH3"],
    "Cellulase":            ["GH5","GH6","GH7","GH9","GH12","GH44","GH45","GH48"],
    "Chitin":               ["GH18","GH19","GH20","AA10"],
    "Starch_alpha_glucan":  ["GH13","GH15","GH31","GH57","GH77","GT5"],
    "Laminarin_beta_glucan":["GH16","GH17","GH55","GH64","GH81","GH128"],
    "Xylan":                ["GH8","GH10","GH11","GH30","GH43","GH51","GH62","GH67","GH115","GH120"],
    "Mannan":               ["GH2","GH26","GH36","GH76","GH92","GH125"],
    "Pectin":               ["PL1","PL2","PL3","PL4","PL10","PL11","GH28","GH78","GH88","GH105"],
    "Alginate":             ["PL5","PL7","PL14","PL15","PL17","PL18","PL36","PL39"],
    # build_cazy_consensus_matrix.py collapses subfamily suffixes, so a subfamily such as
# GH16_3 can never match; agarose is therefore carried by GH86, GH96 and PL6 alone.
    "Agarose":              ["GH86","GH96","PL6"],
    "Peptidoglycan":        ["GH73","GH102","GH103","GH104"],
    "Other_polysaccharide": ["GH23","GH24","GH25","GT2","GT4","GT26","GT51"],
}

STRAIN_MAP = {
    "L5_PB002":"E1","M7_REF":"M7","S_frigidimarina":"frig",
    "S_livingstonensis":"livi","S_polaris":"pola",
    "S_psychromarinicola":"psyc","S_sp002836315":"sp02","S_sp014164505":"sp14",
}

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(INPUT_FILE, sep="\t")
    print(f"Read {len(df)} CAZy families.")

    strain_cols = list(STRAIN_MAP.keys())
    short_names = list(STRAIN_MAP.values())
    df = df.set_index("family")

    substrate_rows, fam_rows = [], []
    for substrate, families in SUBSTRATE_GROUPS.items():
        found = [f for f in families if f in df.index]
        counts = df.loc[found, strain_cols].sum(axis=0) if found else pd.Series(0, index=strain_cols)

        row = {"substrate_group": substrate}
        for col, short in zip(strain_cols, short_names):
            row[short] = int(counts[col])
        present = [s for col, s in zip(strain_cols, short_names) if counts[col] > 0]
        row["present_in_strains"] = ";".join(present) or "none"
        substrate_rows.append(row)

        for fam in families:
            fr = {"substrate_group": substrate, "family": fam, "found_in_matrix": fam in df.index}
            if fam in df.index:
                for col, short in zip(strain_cols, short_names):
                    fr[short] = int(df.loc[fam, col])
                fr["total_count"] = sum(int(df.loc[fam, col]) for col in strain_cols)
            else:
                for short in short_names:
                    fr[short] = 0
                fr["total_count"] = 0
            fam_rows.append(fr)

    result_df = pd.DataFrame(substrate_rows, columns=["substrate_group"]+short_names+["present_in_strains"])
    fam_df    = pd.DataFrame(fam_rows,       columns=["substrate_group","family","found_in_matrix"]+short_names+["total_count"])

    result_df.to_csv(os.path.join(OUTPUT_DIR, "carbohydrate_substrate_matrix.tsv"), sep="\t", index=False)
    fam_df.to_csv(   os.path.join(OUTPUT_DIR, "families_per_group.tsv"),             sep="\t", index=False)

    print("\nSubstrate x strain matrix")
    print(result_df.to_string(index=False))

    print("\nTotals per substrate group")
    for _, row in result_df.iterrows():
        sub = row["substrate_group"]
        total = sum(row[s] for s in short_names)
        e1 = row["E1"]; m7 = row["M7"]
        print(f"  {sub:25s}: total={total:4.0f}  E1={e1}  M7={m7}  present={row['present_in_strains']}")

    print(f"\nDone. Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
