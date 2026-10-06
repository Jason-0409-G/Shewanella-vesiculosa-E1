#!/usr/bin/env python3
"""Keyword mining of symbiosis-related gene categories in the eggNOG annotations of the eight genomes.

Produces the per-category gene counts of Supplementary Data 6 (Supplementary Note 4.6).
"""

import os
import pandas as pd
from collections import defaultdict

BASE = "."   # directory holding the *_eggnog output directories (02_annotation/eggnog_mapper_8genomes.sh)
OUT_DIR = "symbiosis_gene_mining"
os.makedirs(OUT_DIR, exist_ok=True)

STRAINS = ["L5", "M7", "frig", "livi", "pola", "psyc", "sp02", "sp14"]

EGGNOG_FILES = {
    "L5":   os.path.join(BASE, "L5_PB002_eggnog",          "L5_v1156.emapper.annotations"),
    "M7":   os.path.join(BASE, "M7_REF_eggnog",             "M7_v1156.emapper.annotations"),
    "frig": os.path.join(BASE, "S_frigidimarina_eggnog",    "frig_v1156.emapper.annotations"),
    "livi": os.path.join(BASE, "S_livingstonensis_eggnog",  "livi_v1156.emapper.annotations"),
    "pola": os.path.join(BASE, "S_polaris_eggnog",          "pola_v1156.emapper.annotations"),
    "psyc": os.path.join(BASE, "S_psychromarinicola_eggnog","psyc_v1156.emapper.annotations"),
    "sp02": os.path.join(BASE, "S_sp002836315_eggnog",      "sp02_v1156.emapper.annotations"),
    "sp14": os.path.join(BASE, "S_sp014164505_eggnog",      "sp14_v1156.emapper.annotations"),
}

# optional extra eggNOG table of E1 genes, merged into the L5 hits
L5_ONLY_FILE = os.path.join(BASE, "functional_mining_v2", "L5_only_eggnog_v2.tsv")

SYMBIOSIS_KEYWORDS = {
    "T4P_Flp_pilus":  ["cpaB","cpaC","tadB","tadC","PilB","PilQ","fimbriae","Flp","Tad"],
    "MSHA_pilus":     ["MSHA","MshA","MshB","MshC","MshD","MshE"],
    "Curli":          ["curli","csgA","csgB","CsgA","CsgB"],
    "T6SS":           ["T6SS","VgrG","Hcp","VipA","VipB","ImpA","ImpB","ClpV"],
    "Quorum_sensing": ["LuxR","LuxI","autoinducer","quorum"],
    "Iron_acq":       ["TonB","siderophore","enterobactin","aerobactin","IucA","IucC"],
    "Biofilm":        ["biofilm","LapA","LapG","LapD"],
    "Anaerobic_resp": ["fumarate reductase","fccA","TMAO reductase","torA","DMS reductase","dmsA"],
}

COLS = ["query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl",
        "COG_category","Description","Preferred_name","GOs","EC","KEGG_ko",
        "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE",
        "KEGG_TC","CAZy","BiGG_Reaction","PFAMs"]


def parse_eggnog(filepath):
    rows = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(COLS):
                parts += [""] * (len(COLS) - len(parts))
            rows.append(parts[:len(COLS)])
    return pd.DataFrame(rows, columns=COLS)


def match_keywords(df, keywords_dict):
    results = defaultdict(list)
    for _, row in df.iterrows():
        desc  = str(row.get("Description", ""))
        pref  = str(row.get("Preferred_name", ""))
        kegg  = str(row.get("KEGG_ko", ""))
        pfam  = str(row.get("PFAMs", ""))
        query = str(row.get("query", ""))
        text  = f"{desc} {pref} {kegg} {pfam}"
        for cat, kws in keywords_dict.items():
            already_ids = {g for g, _, _, _ in results[cat]}
            if query in already_ids:
                continue
            for kw in kws:
                if kw.lower() in text.lower():
                    results[cat].append((query, kw,
                                         pref if pref not in ("","nan") else "",
                                         desc if desc not in ("","nan") else ""))
                    break
    return results


def main():
    strain_data = {}
    for strain in STRAINS:
        fpath = EGGNOG_FILES[strain]
        strain_data[strain] = parse_eggnog(fpath)

    strain_hits = {s: match_keywords(strain_data[s], SYMBIOSIS_KEYWORDS) for s in STRAINS}

    # Merge hits from the optional extra E1 table into L5
    if os.path.exists(L5_ONLY_FILE):
        l5_only = pd.read_csv(L5_ONLY_FILE, sep="\t")
        for c in COLS:
            if c not in l5_only.columns:
                l5_only[c] = ""
        if "query" not in l5_only.columns:
            l5_only = l5_only.rename(columns={l5_only.columns[0]: "query"})
        extra_hits = match_keywords(l5_only, SYMBIOSIS_KEYWORDS)
        for cat, entries in extra_hits.items():
            existing = {g for g, _, _, _ in strain_hits["L5"][cat]}
            for e in entries:
                if e[0] not in existing:
                    strain_hits["L5"][cat].append(e)

    # Build count table
    rows_out = []
    for cat in SYMBIOSIS_KEYWORDS:
        row = {"category": cat}
        gene_lists = {}
        for s in STRAINS:
            hits = strain_hits[s].get(cat, [])
            row[s] = len(hits)
            gene_lists[s] = ";".join(sorted({pref if pref else q for q,_,pref,_ in hits})) or "-"
        row["E1_genes"] = gene_lists["L5"]
        row["M7_genes"] = gene_lists["M7"]
        rows_out.append(row)

    count_df = pd.DataFrame(rows_out, columns=["category"]+STRAINS+["E1_genes","M7_genes"])
    count_df.to_csv(os.path.join(OUT_DIR, "symbiosis_counts.tsv"), sep="\t", index=False)

    print(count_df.to_string(index=False))
    print(f"\nOutput: {OUT_DIR}")


if __name__ == "__main__":
    main()
