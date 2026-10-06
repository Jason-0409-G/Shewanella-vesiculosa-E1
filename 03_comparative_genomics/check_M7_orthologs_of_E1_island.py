#!/usr/bin/env python3
"""Check whether M7 carries intact orthologs of the E1 island genes.

For each E1 island protein: blastp against the M7 proteome (best hit and length
ratio) and tblastn against the M7 genome (HSP count, query coverage, weighted
identity), which also detects pseudogenized or unannotated loci.
"""
import os, subprocess

# Prokka protein FASTA files (see 02_annotation/prokka_annotate_8genomes.sh) and the M7 genome.
L5_FAA = os.environ.get("L5_FAA", "L5_PB002_prokka/L5.faa")
M7_FAA = os.environ.get("M7_FAA", "M7_REF_prokka/M7.faa")
M7_GENOME = os.environ.get("M7_GENOME", "M7.fna")
OUT = os.environ.get("OUT_DIR", "M7_island_check")
os.makedirs(OUT, exist_ok=True)

ISLAND = {  # E1 locus tag -> gene label
    "PB002_01800": "GH3_e108",
    "PB002_01809": "bglB",
    "PB002_01810": "MFS transporter",
    "PB002_01811": "GH3_e227 (frameshifted in M7)",
}

def read_fasta(path):
    s, n, b = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if n: s[n] = "".join(b)
            n = line[1:].split()[0]; b = []
        else: b.append(line.strip())
    if n: s[n] = "".join(b)
    return s

l5 = read_fasta(L5_FAA)

# E1 island proteins as query
qfa = f"{OUT}/L5_island_proteins.faa"
with open(qfa, "w") as o:
    for loc in ISLAND:
        o.write(f">{loc}\n{l5[loc]}\n")

m7pdb = f"{OUT}/M7_prot"
m7gdb = f"{OUT}/M7_genome"
subprocess.run(["makeblastdb", "-in", M7_FAA, "-dbtype", "prot", "-out", m7pdb], check=True, capture_output=True)
subprocess.run(["makeblastdb", "-in", M7_GENOME, "-dbtype", "nucl", "-out", m7gdb], check=True, capture_output=True)

def blastp(q):
    r = subprocess.run(["blastp", "-query", q, "-db", m7pdb, "-outfmt",
        "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore", "-max_target_seqs", "3"],
        capture_output=True, text=True)
    return [ln.split("\t") for ln in r.stdout.strip().split("\n") if ln]

def tblastn(q):
    r = subprocess.run(["tblastn", "-query", q, "-db", m7gdb, "-outfmt",
        "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore", "-max_target_seqs", "5"],
        capture_output=True, text=True)
    return [ln.split("\t") for ln in r.stdout.strip().split("\n") if ln]

print("M7 island-component integrity check (vs E1 island proteins)")
for loc, label in ISLAND.items():
    qlen = len(l5[loc])
    print(f"\n### {loc}  {label}   (E1 = {qlen} aa)")
    bp = [h for h in blastp(qfa) if h[0] == loc]
    if bp:
        h = bp[0]
        sid, pid, aln, ql, sl, qcov = h[1], float(h[2]), int(h[3]), int(h[4]), int(h[5]), h[6]
        ratio = sl / ql * 100
        print(f"  blastp->M7 proteome: best {sid} | id {pid:.0f}% | qcov {qcov}% | "
              f"M7prot {sl} aa = {ratio:.0f}% of E1 length")
    else:
        print("  blastp->M7 proteome: NO annotated hit")
    # tblastn vs M7 genome
    tb = [h for h in tblastn(qfa) if h[0] == loc]
    if tb:
        best = tb[0]
        sid = best[1]
        same = [h for h in tb if h[1] == sid]
        tot_aln = sum(int(h[3]) for h in same)
        qcov_dna = 100 * tot_aln / qlen
        n_hsp = len(same)
        wpid = sum(float(h[2]) * int(h[3]) for h in same) / tot_aln
        print(f"  tblastn->M7 genome:  {n_hsp} HSP on {sid} | covers {qcov_dna:.0f}% of E1 protein | wgt_id {wpid:.0f}%")
        verdict = ("INTACT ortholog" if (qcov_dna >= 85 and n_hsp <= 2 and wpid >= 60)
                   else "TRUNCATED / split / pseudogenized")
        print(f"  --> {verdict}")
    else:
        print("  tblastn->M7 genome:  NO hit  --> ABSENT")
