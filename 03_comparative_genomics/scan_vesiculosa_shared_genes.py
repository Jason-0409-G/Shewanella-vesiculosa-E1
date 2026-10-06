"""Annotate the orthogroups restricted to S. vesiculosa E1 and M7.

Input   orthofinder_analysis/Sv_species_markers_L5_genes.tsv  from classify_orthogroups_8genomes.py
        L5_PB002_eggnog/L5_v1156.emapper.annotations          from 02_annotation/eggnog_mapper_8genomes.sh
Output  orthofinder_analysis/Sv_species_markers_with_emapper.tsv
The keyword scan below lists the chitin- and amino-sugar-related genes among them
(Supplementary Fig. 1).
"""
from __future__ import annotations
import re
import pandas as pd
from pathlib import Path

ANALYSIS_DIR = Path("orthofinder_analysis")
EMAPPER_PATH = "L5_PB002_eggnog/L5_v1156.emapper.annotations"
EMAPPER_COLS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
    "max_annot_lvl", "COG_category", "Description", "Preferred_name",
    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
    "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
    "BiGG_Reaction", "PFAMs",
]

sv_genes = pd.read_csv(ANALYSIS_DIR / "Sv_species_markers_L5_genes.tsv", sep="\t")
print(f"{len(sv_genes)} E1 genes from {sv_genes['OG'].nunique()} restricted orthogroups")

emap = pd.read_csv(EMAPPER_PATH, sep="\t", comment="#", header=None,
                   names=EMAPPER_COLS, low_memory=False)
sv_with_emap = sv_genes.merge(emap, left_on="L5_gene_id", right_on="query", how="left")
print(f"With eggNOG annotation: {sv_with_emap['query'].notna().sum()} / {len(sv_with_emap)}")

HALLMARKS = {
    "chitin_degradation": {
        "phrases": ["chitinase", "chitobiase", "chitin-binding", "GH18", "GH20"],
        "tokens": ["chiA", "chiB"],
    },
    "AA10_LPMO": {
        "phrases": ["lytic polysaccharide monooxygenase", "LPMO", "GbpA",
                    "AA10", "chitin-binding protein"],
        "tokens": [],
    },
}
SEARCH_COLS = ["Description", "Preferred_name", "PFAMs", "GOs", "KEGG_ko"]


def search_hallmark(row, hm, columns):
    text = " ".join(str(row[c]) for c in columns if pd.notna(row[c]))
    text_lower = text.lower()
    for p in hm["phrases"]:
        if p.lower() in text_lower:
            return True
    for t in hm["tokens"]:
        if re.search(r"\b" + re.escape(t) + r"\b", text):
            return True
    return False


for module, hm in HALLMARKS.items():
    mask = sv_with_emap.apply(lambda r: search_hallmark(r, hm, SEARCH_COLS), axis=1)
    hits = sv_with_emap[mask]
    print(f"\n{module}: {len(hits)} hits")
    for _, r in hits.iterrows():
        pref = str(r["Preferred_name"]) if pd.notna(r["Preferred_name"]) else "-"
        desc = str(r["Description"])[:80] if pd.notna(r["Description"]) else ""
        cog = str(r["COG_category"]) if pd.notna(r["COG_category"]) else "?"
        print(f"  {r['L5_gene_id']:<14} [{cog:<3}] {pref[:15]:<16} {desc}")

keep_cols = ["OG", "L5_gene_id", "COG_category", "Preferred_name",
             "Description", "PFAMs", "KEGG_ko", "EC"]
out = ANALYSIS_DIR / "Sv_species_markers_with_emapper.tsv"
sv_with_emap[keep_cols].to_csv(out, sep="\t", index=False)
print(f"\nSaved {out}")
