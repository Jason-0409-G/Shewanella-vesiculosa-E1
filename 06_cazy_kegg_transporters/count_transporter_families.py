#!/usr/bin/env python3
"""Count genes of twelve transporter/surface-protein families (Pfam matches in eggNOG annotations) in the eight genomes.

Produces the transporter-family counts of Supplementary Data 5 (Supplementary Note 4.5).
"""

import os
import pandas as pd

FAMILIES = {
    "MFS":        ["MFS_1", "MFS_2", "MFS_3", "Sugar_tr", "Lactate_perm"],
    "ABC":        ["ABC_tran", "ABC_membrane", "ABC_membrane_2", "BPD_transp_1", "BPD_transp_2"],
    "EamA_DMT":   ["EamA", "DMT"],
    "TonB_dep":   ["TonB_dep_Rec", "Plug"],
    "T2SS_E":     ["GspE", "T2SSE", "PilB"],
    "Flp_Tad":    ["Flp", "TadB", "CpaB", "CpaC", "TadC"],
    "MSHA":       ["MSHA_pilin", "MshA"],
    "Curli":      ["CsgA", "CsgB"],
    "Autotransp": ["Autotransporter", "Autotransporter_2", "POTRA"],
    "SecYEG_Tat": ["SecY", "TatA", "TatC"],
    "Trk_Kdp":    ["Trk_H_K", "Kdp_ATPase", "KdpA"],
    "OmpA_Porin": ["OmpA", "OmpA_C", "Porin_1", "Gram_neg_porin"],
}

# eggNOG-mapper output of each genome (see 02_annotation/eggnog_mapper_8genomes.sh)
STRAINS = {
    "L5":   "L5_PB002_eggnog/L5_v1156.emapper.annotations",
    "M7":   "M7_REF_eggnog/M7_v1156.emapper.annotations",
    "frig": "S_frigidimarina_eggnog/frig_v1156.emapper.annotations",
    "livi": "S_livingstonensis_eggnog/livi_v1156.emapper.annotations",
    "pola": "S_polaris_eggnog/pola_v1156.emapper.annotations",
    "psyc": "S_psychromarinicola_eggnog/psyc_v1156.emapper.annotations",
    "sp02": "S_sp002836315_eggnog/sp02_v1156.emapper.annotations",
    "sp14": "S_sp014164505_eggnog/sp14_v1156.emapper.annotations",
}
OUT_DIR = "transporter_family_count"

PFAM_COL = 20   # 0-based index of the PFAMs column in .emapper.annotations


def has_pfam_match(pfam_string, family_pfams):
    if not pfam_string or pfam_string.strip() in ("-", "nan"):
        return False
    tokens = [t.strip().lower() for t in str(pfam_string).split(",")]
    fam_lower = [f.lower() for f in family_pfams]
    for token in tokens:
        for fp in fam_lower:
            if fp in token or token in fp:
                return True
    return False


def count_transporters(filepath):
    df = pd.read_csv(filepath, sep="\t", comment="#", header=None,
                     low_memory=False, on_bad_lines="skip")
    pfam_col = df.iloc[:, PFAM_COL].astype(str)
    gene_ids = df.iloc[:, 0].astype(str)
    counts = {}
    for fam, pfams in FAMILIES.items():
        mask = pfam_col.apply(lambda x: has_pfam_match(x, pfams))
        counts[fam] = int(gene_ids[mask].nunique())
    return counts


def main():
    strain_order = list(STRAINS.keys())
    results = {}
    for strain in strain_order:
        results[strain] = count_transporters(STRAINS[strain])

    rows = []
    for fam in FAMILIES:
        row = {"family": fam}
        for s in strain_order:
            row[s] = results[s].get(fam, 0)
        rows.append(row)

    df = pd.DataFrame(rows)
    print()
    print(df.to_string(index=False))

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, "transporter_matrix.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
