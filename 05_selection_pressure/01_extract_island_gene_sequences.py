"""
Extract protein and nucleotide sequences for the PAML target orthogroups.

Each OG gets two FASTA files:
  {og}_prot.faa  protein sequences, input to MAFFT
  {og}_nucl.ffn  nucleotide CDS, input to PAL2NAL

Notes:
  - OG0001992_GH3_NagZ was included in the original run of this script. The gene lies
    outside the island and was not carried into the manuscript; it is omitted here.
  - The MFS transporter (OG0000506) is handled by 04b and 06b.
  - The GH3_e227 alignment behind the reported omega was rebuilt from reassembled
    full-length ORFs (see README.md).
"""
from __future__ import annotations
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(".")
PROKKA = ROOT / "8strains_prokka_v1156"
OF = ROOT / "orthofinder_8strains/input/OrthoFinder/Results_run_v1"
OUT = ROOT / "04_selection_pressure/paml/results"

STRAIN_INFO = {
    "L5":   {"prokka": "L5_PB002_prokka",          "prefix": "L5",   "tag_prefix": "PB002"},
    "M7":   {"prokka": "M7_REF_prokka",             "prefix": "M7",   "tag_prefix": "KDH10"},
    "frig": {"prokka": "S_frigidimarina_prokka",   "prefix": "frig", "tag_prefix": "FRIG"},
    "livi": {"prokka": "S_livingstonensis_prokka", "prefix": "livi", "tag_prefix": "LIVI"},
    "pola": {"prokka": "S_polaris_prokka",          "prefix": "pola", "tag_prefix": "POLA"},
    "psyc": {"prokka": "S_psychromarinicola_prokka","prefix": "psyc", "tag_prefix": "PSYC"},
    "sp02": {"prokka": "S_sp002836315_prokka",     "prefix": "sp02", "tag_prefix": "SP02"},
    "sp14": {"prokka": "S_sp014164505_prokka",     "prefix": "sp14", "tag_prefix": "SP14"},
}

TARGETS = {
    "OG0001451_GH1": {          # bglB
        "og": "OG0001451",
        # sp14 copy is truncated by the contig end (JACJFI010000009.1, 39.5 kb;
        # gene ends 99 bp from the boundary, 89 aa shorter): assembly artefac
        "exclude": ["sp14"],
    },
    "OG0001452_GH3_e227": {
        "og": "OG0001452",
        "exclude": ["M7"],      # M7 copy split across two ORFs (01695 + 01696)
    },
    "OG0002997_GH3_e108": {
        "og": "OG0002997",
        "exclude": ["psyc"],    # no S. psychromarinicola member in this OG
    },
}

def load_orthogroups() -> dict[str, dict[str, list[str]]]:
    """Parse Orthogroups.tsv into nested dict: og -> strain -> [gene_ids]."""
    df = pd.read_csv(OF / "Orthogroups/Orthogroups.tsv", sep="\t")
    result = {}
    for _, row in df.iterrows():
        og = row["Orthogroup"]
        result[og] = {}
        for col in df.columns[1:]:
            cell = row[col]
            if pd.notna(cell) and str(cell).strip():
                result[og][col] = [g.strip() for g in str(cell).split(",")]
            else:
                result[og][col] = []
    return resul


def load_seqs(strain: str, kind: str) -> dict[str, str]:
    """Load all .faa or .ffn for one strain."""
    info = STRAIN_INFO[strain]
    if kind == "prot":
        path = PROKKA / info["prokka"] / f"{info['prefix']}.faa"
    else:
        path = PROKKA / info["prokka"] / f"{info['prefix']}.ffn"
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}


def main():
    og2members = load_orthogroups()

    print("Loading sequences for 8 strains ...")
    prot_cache = {s: load_seqs(s, "prot") for s in STRAIN_INFO}
    nucl_cache = {s: load_seqs(s, "nucl") for s in STRAIN_INFO}
    print("  Loaded.\n")

    for tgt_name, tgt in TARGETS.items():
        og = tgt["og"]
        exclude = set(tgt["exclude"])
        outdir = OUT / tgt_name
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"=== {tgt_name} ({og}) ===")
        print(f"  Excluding strains: {sorted(exclude) if exclude else 'None'}")

        prot_records = []
        nucl_records = []
        members = og2members.get(og, {})

        for strain in STRAIN_INFO:
            if strain in exclude:
                continue
            genes = members.get(strain, [])
            if not genes:
                print(f"  WARNING: {strain}: no gene in {og}")
                continue
            if len(genes) > 1:
                print(f"  WARNING: {strain}: multi-copy ({len(genes)}), taking first")
            gene = genes[0]

            prot_seq = prot_cache[strain].get(gene)
            nucl_seq = nucl_cache[strain].get(gene)

            if prot_seq is None or nucl_seq is None:
                print(f"  WARNING: {strain} {gene}: protein={'OK' if prot_seq else 'MISSING'}, "
                      f"nucl={'OK' if nucl_seq else 'MISSING'}")
                continue

            # CDS length should be 3 x protein length, with or without the stop codon
            expected_3p = 3 * len(prot_seq) + 3
            expected_3p_nostop = 3 * len(prot_seq)
            if len(nucl_seq) not in (expected_3p, expected_3p_nostop):
                print(f"  WARNING: {strain} {gene}: prot={len(prot_seq)}aa, "
                      f"nucl={len(nucl_seq)}bp (expected {expected_3p_nostop} or {expected_3p})")

            # strain prefix in the ID gives readable tree labels
            seq_id = f"{strain}_{gene}"

            prot_records.append(SeqRecord(Seq(prot_seq), id=seq_id, description=""))
            nucl_records.append(SeqRecord(Seq(nucl_seq), id=seq_id, description=""))
            print(f"  {strain:<5} {gene:<14} prot={len(prot_seq)}aa, nucl={len(nucl_seq)}bp")

        prot_out = outdir / f"{tgt_name}_prot.faa"
        nucl_out = outdir / f"{tgt_name}_nucl.ffn"
        SeqIO.write(prot_records, prot_out, "fasta")
        SeqIO.write(nucl_records, nucl_out, "fasta")
        print(f"  [saved] {prot_out} ({len(prot_records)} seqs)")
        print(f"  [saved] {nucl_out} ({len(nucl_records)} seqs)\n")

    print("Done. Next step: 02_mafft_linsi_protein_align.sh")


if __name__ == "__main__":
    main()
