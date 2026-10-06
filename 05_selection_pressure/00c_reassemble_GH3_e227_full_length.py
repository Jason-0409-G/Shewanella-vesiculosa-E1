#!/usr/bin/env python3
"""Reassemble the full-length GH3_e227 coding sequence in S. vesiculosa E1 and M7.

Both genomes carry one spurious extra base inside a homopolymer run, a long-read
assembly artifact, which shifts the reading frame of GH3_e227. In E1 the annotated
ORF is therefore truncated to 838 aa, and in M7 the gene is split across two ORFs.
Deleting the extra base restores the frame in both, giving the 873-aa (E1) and
878-aa (M7) proteins that the eight-taxon codon alignment and the reported
dN/dS of 0.0759 are based on.

Corrections applied, as recorded in the provenance note of the 2026-07-16 run:

  E1  contig_1, plus strand.  AAA homopolymer at 2,017,760; delete one A.
      Coding sequence then runs 2,622 nt from 2,015,263, the annotated start of
      PB002_01811, and reads through what Prokka had called PB002_01812.

  M7  CP073588.1, minus strand.  TTTTTTTT run at 1,828,099-1,828,106 on the plus
      strand (an A-run on the coding strand) at the junction between KDH10_01696
      and KDH10_01695; delete the base at 1,828,106. Coding sequence then runs
      2,637 nt over 1,826,433-1,829,070, starting 18 bp upstream of the annotated
      start of KDH10_01696.

This was originally done by hand. The script makes it reproducible and asserts the
result against the stored sequences, so a failed assertion means an input changed.

Input   L5.fna, M7.fna     the assembled chromosomes
Output  E1_e227_reassembled_CDS.fna, M7_e227_reassembled_CDS.fna
        and the translated proteins, for the eight-taxon alignment

With --check <dir> the two CDS are compared against existing copies.
"""
import argparse
import sys
from pathlib import Path

from Bio.Seq import Seq

# genome, 1-based position of the base to delete, 1-based CDS span on the plus
# strand after that deletion, strand, expected protein length
CORRECTIONS = {
    "E1": dict(fasta="L5.fna", delete_at=2_017_760, start=2_015_263,
               length=2_622, strand="+", protein_aa=873, locus="PB002_01811"),
    "M7": dict(fasta="M7.fna", delete_at=1_828_106, start=1_826_433,
               length=2_637, strand="-", protein_aa=878, locus="KDH10_01696+01695"),
}


def read_chromosome(path: Path) -> str:
    return "".join(l.strip() for l in path.read_text().splitlines()
                   if not l.startswith(">")).upper()


def reassemble(genome: str, spec: dict) -> tuple[str, str]:
    pos = spec["delete_at"]
    removed = genome[pos - 1]
    corrected = genome[:pos - 1] + genome[pos:]
    # the deletion shifts every base after `pos` one to the left
    start0 = spec["start"] - 1
    cds = corrected[start0:start0 + spec["length"]]
    if spec["strand"] == "-":
        cds = str(Seq(cds).reverse_complement())
    protein = str(Seq(cds).translate(table=11)).rstrip("*")
    if len(protein) != spec["protein_aa"]:
        sys.exit(f"expected {spec['protein_aa']} aa, got {len(protein)}")
    if "*" in protein:
        sys.exit(f"internal stop codon in the reassembled protein")
    return cds, protein, removed


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-dir", type=Path, default=Path("."))
    ap.add_argument("--out-dir", type=Path, default=Path("."))
    ap.add_argument("--check", type=Path, help="directory holding existing *_reassembled_CDS.fna")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for name, spec in CORRECTIONS.items():
        genome = read_chromosome(args.genome_dir / spec["fasta"])
        cds, protein, removed = reassemble(genome, spec)
        out = args.out_dir / f"{name}_e227_reassembled_CDS.fna"
        out.write_text(f">{name}_GH3_e227_reassembled {spec['locus']}\n{cds}\n")
        print(f"{name}: deleted {removed} at {spec['delete_at']:,} -> "
              f"{len(cds)} nt, {len(protein)} aa ({spec['locus']})")

        if args.check:
            ref_path = args.check / f"{name}_e227_reassembled_CDS.fna"
            ref = "".join(l.strip() for l in ref_path.read_text().splitlines()
                          if not l.startswith(">")).upper()
            print(f"  check: {'identical to' if cds == ref else 'DIFFERS from'} {ref_path.name}")
            if cds != ref:
                sys.exit(1)


if __name__ == "__main__":
    main()
