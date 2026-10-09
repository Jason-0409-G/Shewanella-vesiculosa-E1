#!/usr/bin/env python3
"""Build the CAZy family consensus matrix (families x 8 strains) from dbCAN3 output.

Consensus rule:
  - keep a protein only if at least 2 of the 3 tools support it (column "#ofTools")
  - split the "Recommend Results" column on "|" (multi-domain separator), so that a family carried
    as a non-first domain (e.g. CBM41|CBM48|GH13_13) is counted
  - strip the subfamily suffix (GH13_18 -> GH13, GH13_e144 -> GH13)
  - each protein contributes at most +1 to each base family it carries

Input   8strains_dbcan_v528/{strain}_dbcan/overview.tsv, written by
        02_annotation/dbcan_cazyme_8genomes.sh
Output  all_CAZy_family_matrix.tsv (header: family, L5_PB002, M7_REF, S_frigidimarina, ...),
        read by count_carbohydrate_substrate_genes.py and used for the Fig. 3c heatmap
"""
import re
from pathlib import Path

DBCAN = Path("8strains_dbcan_v528")
OUT_FILE = Path("all_CAZy_family_matrix.tsv")

# dbCAN directory -> strain column name (L5_PB002 = S. vesiculosa E1)
STRAINS = [
    ("L5_PB002_dbcan", "L5_PB002"),
    ("M7_REF_dbcan", "M7_REF"),
    ("S_frigidimarina_dbcan", "S_frigidimarina"),
    ("S_livingstonensis_dbcan", "S_livingstonensis"),
    ("S_polaris_dbcan", "S_polaris"),
    ("S_psychromarinicola_dbcan", "S_psychromarinicola"),
    ("S_sp002836315_dbcan", "S_sp002836315"),
    ("S_sp014164505_dbcan", "S_sp014164505"),
]

MIN_TOOLS = 2
FAM_RE = re.compile(r"^([A-Z]+\d+)")


def families_in_protein(recommend: str) -> set:
    fams = set()
    if recommend in ("-", "", None):
        return fams
    for domain in recommend.split("|"):
        m = FAM_RE.match(domain.strip())
        if m:
            fams.add(m.group(1))
    return fams


def count_strain(overview: Path) -> dict:
    counts = {}
    with open(overview) as fh:
        fh.readline()
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            try:
                if int(cols[5]) < MIN_TOOLS:
                    continue
            except ValueError:
                continue
            for f in families_in_protein(cols[6]):
                counts[f] = counts.get(f, 0) + 1
    return counts


def fam_sort_key(f):
    m = re.match(r"([A-Z]+)(\d+)", f)
    return (m.group(1), int(m.group(2))) if m else (f, 0)


def main():
    per_strain = {}
    all_fams = set()
    for d, label in STRAINS:
        c = count_strain(DBCAN / d / "overview.tsv")
        per_strain[label] = c
        all_fams.update(c.keys())

    fams = sorted(all_fams, key=fam_sort_key)
    labels = [lab for _, lab in STRAINS]

    lines = ["family\t" + "\t".join(labels)]
    for f in fams:
        lines.append(f + "\t" + "\t".join(str(per_strain[l].get(f, 0)) for l in labels))
    OUT_FILE.write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nWritten: {OUT_FILE}")


if __name__ == "__main__":
    main()
