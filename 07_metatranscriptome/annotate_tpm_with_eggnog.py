#!/usr/bin/env python3
"""Join the eggNOG-mapper annotation into the E1 TPM matrix (writes a new file; the input is not modified).

- The Prokka columns (gene, EC_number, COG, product) are kept.
- eggNOG columns (name, EC, KO, COG category, CAZy, Pfam, description) are inserted after "product".
- EC_flag marks Prokka EC numbers that conflict with the eggNOG EC (e.g. PB002_01811: Prokka 3.2.1.52,
  eggNOG 3.2.1.21); eggNOG is taken as authoritative.
- Rows without an eggNOG annotation (rRNA, tRNA, ...) get "-" in the eggNOG columns.

Usage: python annotate_tpm_with_eggnog.py [tpm.tsv] [eggnog.annotations] [out.tsv]
"""
import csv, sys

TPM   = sys.argv[1] if len(sys.argv) > 1 else "L5_TPM_matrix.tsv"
EGG   = sys.argv[2] if len(sys.argv) > 2 else "L5_PB002_eggnog/L5_v1156.emapper.annotations"
OUT   = sys.argv[3] if len(sys.argv) > 3 else "L5_TPM_matrix_eggNOG_annotated.tsv"

EMAPPER_COLS = ["query","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl",
                "COG_category","Description","Preferred_name","GOs","EC","KEGG_ko",
                "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE",
                "KEGG_TC","CAZy","BiGG_Reaction","PFAMs"]

def parse_eggnog(path):
    egg, header = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#query"):
                header = line[1:].split("\t"); continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            cols = header if header else EMAPPER_COLS
            row = dict(zip(cols, parts))
            q = row.get("query", parts[0])
            egg[q] = {
                "name":  row.get("Preferred_name", "-"),
                "ec":    row.get("EC", "-"),
                "ko":    row.get("KEGG_ko", "-"),
                "cog":   row.get("COG_category", "-"),
                "cazy":  row.get("CAZy", "-"),
                "pfams": row.get("PFAMs", "-"),
                "desc":  row.get("Description", "-"),
            }
    return egg

def ec_set(x):
    return {e.strip() for e in (x or "").replace(";", ",").split(",")
            if e.strip() and e.strip() != "-"}

egg = parse_eggnog(EGG)

with open(TPM) as f:
    rows = list(csv.reader(f, delimiter="\t"))

head = rows[0]
if any(c.startswith("eggNOG_") for c in head):
    sys.exit(f"Input already contains eggNOG columns: {TPM}\n"
             f"Pass the un-annotated TPM matrix as the first argument.")
i_locus = head.index("locus_tag")
i_ec    = head.index("EC_number")
i_prod  = head.index("product")
newcols = ["eggNOG_name","eggNOG_EC","eggNOG_KO","eggNOG_COGcat",
           "eggNOG_CAZy","eggNOG_PFAMs","eggNOG_desc","EC_flag"]

out = [head[:i_prod+1] + newcols + head[i_prod+1:]]
conflicts, annotated = [], 0
for row in rows[1:]:
    locus = row[i_locus]
    pec   = row[i_ec]
    e = egg.get(locus)
    if e:
        annotated += 1
        eec = e["ec"]
        pset, eset = ec_set(pec), ec_set(eec)
        if pset and eset:
            if pset <= eset:                       # Prokka EC contained in the eggNOG set
                flag = "ok" if pset == eset else "eggNOG-superset"
            else:                                  # Prokka EC not in the eggNOG set: conflict
                flag = "CONFLICT"
                conflicts.append((locus, pec, eec, e["name"], e["desc"][:55]))
        elif eset and not pset:
            flag = "eggNOG-only"
        else:
            flag = ""
        vals = [e["name"], e["ec"], e["ko"], e["cog"], e["cazy"], e["pfams"], e["desc"], flag]
    else:
        vals = ["-","-","-","-","-","-","-",""]
    out.append(row[:i_prod+1] + vals + row[i_prod+1:])

with open(OUT, "w", newline="") as f:
    csv.writer(f, delimiter="\t").writerows(out)

print(f"Written: {OUT}")
print(f"Rows: {len(rows)-1} | with eggNOG annotation: {annotated} | EC conflicts: {len(conflicts)}")
print("\nProkka vs eggNOG EC conflicts (eggNOG taken as authoritative):")
for locus, pec, eec, name, desc in conflicts:
    print(f"  {locus}: Prokka EC {pec:<10} -> eggNOG EC {eec:<10} ({name}) {desc}")
