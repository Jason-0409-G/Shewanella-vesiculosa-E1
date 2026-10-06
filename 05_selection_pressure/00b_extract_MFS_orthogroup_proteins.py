#!/usr/bin/env python3
"""Write the protein FASTA for the MFS transporter orthogroup (OG0000506).

The MFS transporter was added to the selection-pressure analysis after the gene
list of 01_extract_island_gene_sequences.py was fixed, so its protein FASTA was
originally produced outside that script. This reproduces it: take the members
OrthoFinder assigned to OG0000506 and pull their sequences from the same input
proteomes OrthoFinder was given. Headers carry the gene ID only, with no strain
prefix, which is how 04b and 06b expect the file.

Input   Orthogroups.tsv            OrthoFinder output
        <input_dir>/{strain}.faa   the proteomes OrthoFinder was run on
Output  OG0000506_MFS_prot.faa

With --check <file> the result is compared against an existing copy and the
script exits non-zero if they differ.
"""
import argparse
import sys
from pathlib import Path

from Bio import SeqIO

ORTHOGROUP = "OG0000506"
OUT_NAME = "OG0000506_MFS_prot.faa"


def members(tsv: Path, og: str) -> list[str]:
    for line in tsv.read_text().splitlines():
        if line.startswith(og + "\t"):
            cells = line.rstrip("\n").split("\t")[1:]
            return [g.strip() for cell in cells for g in cell.split(",") if g.strip()]
    sys.exit(f"{og} not found in {tsv}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthogroups", type=Path, default=Path("Orthogroups.tsv"))
    ap.add_argument("--input-dir", type=Path, default=Path("."),
                    help="directory holding the proteomes OrthoFinder was run on")
    ap.add_argument("--out", type=Path, default=Path(OUT_NAME))
    ap.add_argument("--check", type=Path, help="compare against an existing copy")
    args = ap.parse_args()

    wanted = members(args.orthogroups, ORTHOGROUP)
    print(f"{ORTHOGROUP}: {len(wanted)} members")

    pool = {}
    for faa in sorted(args.input_dir.glob("*.faa")):
        for rec in SeqIO.parse(faa, "fasta"):
            pool[rec.id] = str(rec.seq)

    missing = [g for g in wanted if g not in pool]
    if missing:
        sys.exit(f"not found in {args.input_dir}: {', '.join(missing)}")

    text = "".join(f">{g}\n{pool[g]}\n" for g in wanted)
    args.out.write_text(text)
    for g in wanted:
        print(f"  {g}  {len(pool[g])} aa")
    print(f"written: {args.out}")

    if args.check:
        ref = {r.id: str(r.seq) for r in SeqIO.parse(args.check, "fasta")}
        new = {r.id: str(r.seq) for r in SeqIO.parse(args.out, "fasta")}
        if ref == new:
            print(f"check: identical to {args.check}")
        else:
            sys.exit(f"check FAILED against {args.check}")


if __name__ == "__main__":
    main()
