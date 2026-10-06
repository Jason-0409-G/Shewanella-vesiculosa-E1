"""Genome features of the eight Shewanella genomes (Supplementary Data 1).

Computed from the Prokka v1.15.6 GenBank and CDS FASTA files; all contigs of draft
assemblies are included.
"""
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from pathlib import Path
import pandas as pd
from collections import Counter
import numpy as np

PROKKA_ROOT = Path("8strains_prokka_v1156")
OUTPUT_DIR = Path("genome_features")
OUTPUT_DIR.mkdir(exist_ok=True)

STRAINS = [
    ("E1",                  "L5_PB002_prokka",            "L5"),
    ("M7_REF",              "M7_REF_prokka",              "M7"),
    ("S_frigidimarina",     "S_frigidimarina_prokka",     "frig"),
    ("S_livingstonensis",   "S_livingstonensis_prokka",   "livi"),
    ("S_polaris",           "S_polaris_prokka",           "pola"),
    ("S_psychromarinicola", "S_psychromarinicola_prokka", "psyc"),
    ("S_sp002836315",       "S_sp002836315_prokka",       "sp02"),
    ("S_sp014164505",       "S_sp014164505_prokka",       "sp14"),
]


def gc_skew(seq):
    """Genome-wide GC skew, (G - C) / (G + C)."""
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    return (g - c) / (g + c) if (g + c) > 0 else 0.0


def gc3_content(cds_records):
    """GC content at codon third positions (%)."""
    third_pos = []
    for rec in cds_records:
        seq = str(rec.seq).upper()
        for i in range(2, len(seq), 3):
            third_pos.append(seq[i])
    if not third_pos:
        return 0.0
    n = Counter(third_pos)
    gc = n["G"] + n["C"]
    total = sum(n.values())
    return gc / total * 100 if total > 0 else 0.0


def features(strain_label, prokka_dir, prefix):
    """Compute the feature row of one genome."""
    gbk = PROKKA_ROOT / prokka_dir / f"{prefix}.gbk"
    ffn = PROKKA_ROOT / prokka_dir / f"{prefix}.ffn"

    records = list(SeqIO.parse(gbk, "genbank"))
    n_contigs = len(records)
    seq = "".join(str(r.seq).upper() for r in records)

    cds_features = [f for r in records for f in r.features if f.type == "CDS"]
    cds_lens = [len(f) for f in cds_features]
    coding_bp = sum(cds_lens)

    length = len(seq)
    if length == 0:
        return None

    gc_count = seq.count("G") + seq.count("C")
    gc = gc_count / length * 100

    n_cds = len(cds_features)
    n_trna = sum(1 for r in records for f in r.features if f.type == "tRNA")
    n_rrna = sum(1 for r in records for f in r.features if f.type == "rRNA")

    skew = gc_skew(seq) * 100
    intergenic_bp = length - coding_bp
    intergenic_pct = intergenic_bp / length * 100

    cds_records = list(SeqIO.parse(ffn, "fasta")) if ffn.exists() else []
    gc3 = gc3_content(cds_records) if cds_records else None

    cds_mean = np.mean(cds_lens) if cds_lens else 0
    cds_median = np.median(cds_lens) if cds_lens else 0
    cds_std = np.std(cds_lens) if cds_lens else 0

    return {
        "strain": strain_label,
        "n_contigs": n_contigs,
        "length_bp": length,
        "GC_pct": round(gc, 2),
        "GC_skew_pct": round(skew, 3),
        "GC3_pct": round(gc3, 2) if gc3 else None,
        "CDS": n_cds,
        "CDS_per_Mb": round(n_cds / (length / 1e6), 1),
        "CDS_mean_len_bp": round(cds_mean, 1),
        "CDS_median_len_bp": round(cds_median, 1),
        "CDS_stdev_len_bp": round(cds_std, 1),
        "coding_density_pct": round(coding_bp / length * 100, 2),
        "intergenic_pct": round(intergenic_pct, 2),
        "tRNA": n_trna,
        "rRNA": n_rrna,
    }


if __name__ == "__main__":
    rows = [features(*s) for s in STRAINS]
    rows = [r for r in rows if r is not None]
    df = pd.DataFrame(rows)

    print(df.to_string(index=False))

    out = OUTPUT_DIR / "8strains_genome_features.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\n[saved] {out}")
