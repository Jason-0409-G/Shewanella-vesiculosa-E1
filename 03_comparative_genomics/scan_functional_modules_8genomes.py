"""Keyword scan of eggNOG-mapper annotations across the eight Shewanella genomes.

Modules scanned: Flp/Tad pilus, anti-phage toxin-antitoxin systems, chitin degradation.
Input   {label}.emapper.annotations for each genome, produced by
        02_annotation/eggnog_mapper_8genomes.sh (one directory per genome, as below)
Output  8strain_emapper_compare/8strain_hallmark_matrix.tsv   hits per module and genome
        8strain_emapper_compare/8strain_hallmark_details.tsv  the individual hits
"""
from __future__ import annotations
import pandas as pd
import re
from pathlib import Path

EMAPPER_FILES = {
    "E1":                  "L5_PB002_eggnog/L5_v1156.emapper.annotations",
    "M7":                  "M7_REF_eggnog/M7_v1156.emapper.annotations",
    "S_frigidimarina":     "S_frigidimarina_eggnog/frig_v1156.emapper.annotations",
    "S_livingstonensis":   "S_livingstonensis_eggnog/livi_v1156.emapper.annotations",
    "S_polaris":           "S_polaris_eggnog/pola_v1156.emapper.annotations",
    "S_psychromarinicola": "S_psychromarinicola_eggnog/psyc_v1156.emapper.annotations",
    "S_sp002836315":       "S_sp002836315_eggnog/sp02_v1156.emapper.annotations",
    "S_sp014164505":       "S_sp014164505_eggnog/sp14_v1156.emapper.annotations",
}
STRAIN_ORDER = list(EMAPPER_FILES.keys())

EMAPPER_COLS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
    "max_annot_lvl", "COG_category", "Description", "Preferred_name",
    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
    "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
    "BiGG_Reaction", "PFAMs",
]

OUTPUT_DIR = Path("8strain_emapper_compare")
OUTPUT_DIR.mkdir(exist_ok=True)

# Case-insensitive substring matches ("phrases") and whole-word matches ("tokens").
HALLMARKS = {
    "Flp_Tad_pilus": {
        "phrases": ["Flp pilus", "Flp pilus assembly", "tight adherence", "Tad pilus",
                    "GSP D family", "RcpC", "secretion system protein E"],
        "tokens": ["cpaB", "cpaC", "tadB", "tadC", "tadE", "tadG"],
    },
    "TA_phage_defense": {
        "phrases": ["abortive phage resistance", "AbiEii", "AbiGi",
                    "type II toxin-antitoxin", "RnlA"],
        "tokens": ["rnlA"],
    },
    "chitin_degradation": {
        "phrases": ["chitinase", "chitobiase", "N-acetyl-beta-D-glucosaminidase"],
        "tokens": ["chiA", "chiB"],
    },
}

SEARCH_COLS = ["Description", "Preferred_name", "PFAMs", "GOs", "KEGG_ko", "BRITE", "EC"]


def load_emapper(strain: str, path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", header=None,
                     names=EMAPPER_COLS, low_memory=False)
    df["strain"] = strain
    return df


def search_hallmark(row: pd.Series, hallmark: dict, columns: list) -> bool:
    text = " ".join(str(row[c]) for c in columns if pd.notna(row[c]))
    text_lower = text.lower()
    for phrase in hallmark["phrases"]:
        if phrase.lower() in text_lower:
            return True
    for token in hallmark["tokens"]:
        if re.search(r"\b" + re.escape(token) + r"\b", text):
            return True
    return False


combined = pd.concat([load_emapper(s, p) for s, p in EMAPPER_FILES.items()],
                     ignore_index=True)

counts, details = {}, {}
for module, hm in HALLMARKS.items():
    counts[module], details[module] = {}, {}
    for strain in STRAIN_ORDER:
        sub = combined[combined["strain"] == strain]
        hits = sub[sub.apply(lambda r: search_hallmark(r, hm, SEARCH_COLS), axis=1)]
        counts[module][strain] = len(hits)
        details[module][strain] = hits[["query", "Preferred_name", "Description"]]

matrix_df = pd.DataFrame([{"module": m, **c} for m, c in counts.items()])
print(matrix_df.to_string(index=False))
matrix_df.to_csv(OUTPUT_DIR / "8strain_hallmark_matrix.tsv", sep="\t", index=False)

detailed_rows = []
for module in HALLMARKS:
    for strain in STRAIN_ORDER:
        for _, h in details[module][strain].iterrows():
            detailed_rows.append({
                "module": module,
                "strain": strain,
                "query": h["query"],
                "Preferred_name": h["Preferred_name"],
                "Description": str(h["Description"])[:100] if pd.notna(h["Description"]) else "",
            })
detail_df = pd.DataFrame(detailed_rows)
detail_df.to_csv(OUTPUT_DIR / "8strain_hallmark_details.tsv", sep="\t", index=False)
print(f"{len(detail_df)} hits written to {OUTPUT_DIR}")
