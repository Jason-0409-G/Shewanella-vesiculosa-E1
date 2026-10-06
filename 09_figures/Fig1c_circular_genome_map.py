#!/usr/bin/env python3
"""Fig. 1c - circular map of the closed S. vesiculosa E1 chromosome.

Tracks, outermost first: coordinate scale, CDS on the plus strand, CDS on the
minus strand, GC content and GC skew, both deviations from the genome mean.

Inputs
------
L5.fna  closed chromosome, one contig (Prokka input/output)
L5.gff  Prokka annotation of the same assembly

Output
------
Fig1c_circular_genome_map.{svg,png}

The published panel was relabelled and recoloured in a vector editor; this script
draws the data content (track order, feature positions, window statistics), not
that final styling.
"""
import sys
from pathlib import Path

from pycirclize import Circos
from pycirclize.parser import Gff

FNA = Path(sys.argv[1] if len(sys.argv) > 1 else "L5.fna")
GFF = Path(sys.argv[2] if len(sys.argv) > 2 else "L5.gff")
OUT = Path("Fig1c_circular_genome_map")

WINDOW, STEP = 5000, 2500          # GC content / skew window and step, in bp
COL_PLUS, COL_MINUS = "#C04E2C", "#BEA06B"
COL_GC, COL_SKEW_HI, COL_SKEW_LO = "#7F1C1F", "#C59B2B", "#535454"


def read_sequence(path: Path) -> str:
    return "".join(l.strip() for l in path.read_text().splitlines()
                   if not l.startswith(">")).upper()


def sliding(seq: str, window: int, step: int):
    """Yield (window midpoint, GC fraction, GC skew) for each window."""
    n = len(seq)
    for start in range(0, n, step):
        chunk = seq[start:start + window]
        if len(chunk) < window // 2:
            break
        g, c = chunk.count("G"), chunk.count("C")
        gc = (g + c) / len(chunk)
        skew = (g - c) / (g + c) if (g + c) else 0.0
        yield start + len(chunk) / 2, gc, skew


def main() -> None:
    seq = read_sequence(FNA)
    size = len(seq)
    gff = Gff(str(GFF))
    print(f"{FNA.name}: {size:,} bp")

    circos = Circos(sectors={gff.name: size})
    sector = circos.sectors[0]

    # -- scale -------------------------------------------------------------
    scale = sector.add_track((97, 100))
    scale.axis(fc="#F5EEE8")
    scale.xticks_by_interval(
        1_000_000, outer=False, show_bottom_line=True,
        label_formatter=lambda v: f"{v / 1_000_000:.0f} Mb",
        label_orientation="vertical", line_kws=dict(ec="#535454"),
    )
    scale.xticks_by_interval(100_000, outer=False, tick_length=1,
                             show_label=False, line_kws=dict(ec="#888888"))

    # -- CDS, by strand ----------------------------------------------------
    for track_range, strand, colour in (((86, 94), 1, COL_PLUS),
                                        ((77, 85), -1, COL_MINUS)):
        track = sector.add_track(track_range)
        feats = gff.extract_features("CDS", target_strand=strand)
        track.genomic_features(feats, plotstyle="arrow", fc=colour, lw=0)
        print(f"  CDS strand {strand:+d}: {len(feats)}")

    # -- GC content, as deviation from the genome mean ---------------------
    stats = list(sliding(seq, WINDOW, STEP))
    x = [p for p, _, _ in stats]
    gc_mean = (seq.count("G") + seq.count("C")) / size
    gc_dev = [gc - gc_mean for _, gc, _ in stats]
    print(f"  genome GC {gc_mean * 100:.2f} %   windows {len(stats)}")

    gc_track = sector.add_track((62, 75))
    lim = max(abs(v) for v in gc_dev)
    gc_track.fill_between(x, gc_dev, 0, vmin=-lim, vmax=lim, fc=COL_GC, lw=0)

    # -- GC skew -----------------------------------------------------------
    skew = [s for _, _, s in stats]
    sk_track = sector.add_track((47, 60))
    lim = max(abs(v) for v in skew)
    sk_track.fill_between(x, [max(v, 0) for v in skew], 0,
                          vmin=-lim, vmax=lim, fc=COL_SKEW_HI, lw=0)
    sk_track.fill_between(x, [min(v, 0) for v in skew], 0,
                          vmin=-lim, vmax=lim, fc=COL_SKEW_LO, lw=0)

    # text must be registered before plotfig(), which renders all elements
    circos.text(f"Shewanella vesiculosa E1\n{size / 1e6:.2f} Mb",
                size=13, r=0, color="#222222")
    fig = circos.plotfig()
    circos.ax.legend(
        handles=[
            circos.ax.plot([], [], marker="s", ls="", color=COL_PLUS, label="CDS (+)")[0],
            circos.ax.plot([], [], marker="s", ls="", color=COL_MINUS, label="CDS (-)")[0],
            circos.ax.plot([], [], marker="s", ls="", color=COL_GC, label="GC content")[0],
            circos.ax.plot([], [], marker="s", ls="", color=COL_SKEW_HI, label="GC skew")[0],
        ],
        loc="lower right", bbox_to_anchor=(1.1, 0.0), fontsize=8, frameon=False,
    )
    for ext in ("svg", "png"):
        fig.savefig(f"{OUT}.{ext}", dpi=300 if ext == "png" else None,
                    bbox_inches="tight")
    print(f"written: {OUT}.svg / .png")


if __name__ == "__main__":
    main()
