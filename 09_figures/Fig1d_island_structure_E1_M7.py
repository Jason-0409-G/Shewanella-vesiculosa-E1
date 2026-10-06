#!/usr/bin/env python3
"""Fig. 1d - the beta-glucoside-processing island in S. vesiculosa E1 and M7.

Draws the six island genes of each genome to scale on two rows, with the
intervening region between GH3_e108 and galR collapsed and labelled by its length.
Coordinates, strands and lengths are read from the Prokka GFF files, so the
drawn sizes are the annotated ones and are not hard-coded.

Gene symbols follow the manuscript, not the Prokka product names: Prokka calls
PB002_01808 "yicJ", PB002_01810 "fucP_1" and PB002_01811 "nagZ_1"; the manuscript
uses bglT (assigned by blastp against S. baltica BglT), "MFS transporter" and
GH3_e227 respectively.

In M7 the GH3_e227 ortholog is split across two ORFs (KDH10_01696 and
KDH10_01695) by a frameshift; only KDH10_01696 is drawn, labelled by locus tag,
but the end coordinate printed under the M7 row also covers KDH10_01695. The M7
island is in the opposite orientation, so its axis runs from high to low
coordinates.

Inputs
------
L5.gff, M7.gff   Prokka annotations

Output
------
Fig1d_island_structure_E1_M7.{svg,png}
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

matplotlib.rcParams.update({"svg.fonttype": "none", "font.family": "Arial",
                            "font.size": 8})

L5_GFF = Path(sys.argv[1] if len(sys.argv) > 1 else "L5.gff")
M7_GFF = Path(sys.argv[2] if len(sys.argv) > 2 else "M7.gff")
OUT = Path("Fig1d_island_structure_E1_M7")

# published symbol -> locus tag, in the order drawn (left to right)
E1_GENES = [("GH3_e108", "PB002_01800"), ("galR", "PB002_01807"),
            ("bglT", "PB002_01808"), ("bglB", "PB002_01809"),
            ("MFS", "PB002_01810"), ("GH3_e227", "PB002_01811")]
M7_GENES = [("GH3_e108", "KDH10_01707"), ("galR", "KDH10_01700"),
            ("bglT", "KDH10_01699"), ("bglB", "KDH10_01698"),
            ("MFS", "KDH10_01697"), ("01696", "KDH10_01696")]

M7_SPAN_ALSO = ["KDH10_01695"]   # second ORF of M7 GH3_e227: counted in the span, not drawn
# Gene colours here follow the Fig. 4b palette, not the brown-gold scheme of the
# published Fig. 1d, which was recoloured in the vector editor.
COLOURS = {"GH3_e108": "#266C9A", "galR": "#8B5E3C", "bglT": "#D78C3C",
           "bglB": "#1C8458", "MFS": "#763A8E",
           "GH3_e227": "#B6392B", "01696": "#B6392B"}
GAP_FRAC = 0.12          # share of the axis given to the collapsed gap


def read_gff(path: Path) -> dict:
    """locus tag -> (start, end, strand); 1-based inclusive, as in the GFF."""
    out = {}
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("#") or "\tCDS\t" not in line:
            continue
        f = line.split("\t")
        m = re.search(r"ID=([A-Za-z0-9]+_\d+)", f[8])
        if m:
            out[m.group(1)] = (int(f[3]), int(f[4]), f[6])
    return out


def layout(genes, ann, reverse):
    """Map genomic coordinates onto a 0-1 axis with the long gap collapsed."""
    coords = [ann[tag] for _, tag in genes]
    gap_from, gap_to = coords[0][1], coords[1][0]           # GH3_e108 -> galR
    if reverse:
        gap_from, gap_to = coords[1][1], coords[0][0]
    gap_len = abs(gap_to - gap_from)

    cluster = coords[1:]
    c_lo = min(a for a, b, s in cluster)
    c_hi = max(b for a, b, s in cluster)
    drawn = (coords[0][1] - coords[0][0]) + (c_hi - c_lo)   # bp actually drawn
    left_frac = (coords[0][1] - coords[0][0]) / drawn * (1 - GAP_FRAC)

    boxes = []
    for (label, tag), (a, b, strand) in zip(genes, coords):
        if (label, tag) == genes[0]:
            x0, w = 0.0, left_frac
        else:
            span = c_hi - c_lo
            off = (c_hi - b) if reverse else (a - c_lo)
            x0 = left_frac + GAP_FRAC + off / span * (1 - GAP_FRAC - left_frac)
            w = (b - a) / span * (1 - GAP_FRAC - left_frac)
        pointing_right = (strand == "+") != reverse
        boxes.append((label, tag, x0, w, b - a + 1, pointing_right))
    return boxes, gap_len, (c_lo, c_hi), coords[0]


def draw_row(ax, y, name, boxes, gap_len, span_lo, span_hi, reverse):
    ax.plot([0, 1], [y, y], color="#999999", lw=0.8, zorder=1)
    ax.text(-0.035, y, name, ha="right", va="center", fontsize=9, fontweight="bold")

    for label, tag, x0, w, length, right in boxes:
        head = min(w * 0.35, 0.018)
        ax.add_patch(FancyArrow(
            x0 if right else x0 + w, y, (w - head) if right else -(w - head), 0,
            width=0.30, head_width=0.30, head_length=head,
            length_includes_head=True, lw=0,
            fc=COLOURS[label], zorder=3))
        ax.text(x0 + w / 2, y - 0.42, label, ha="center", va="top",
                fontsize=7.5, style="italic" if label.islower() else "normal")
        ax.text(x0 + w / 2, y + 0.42, f"{length / 1000:.1f} kb",
                ha="center", va="bottom", fontsize=6.5, color="#555555")

    gx = boxes[0][2] + boxes[0][3] + GAP_FRAC / 2
    ax.plot([boxes[0][2] + boxes[0][3], boxes[1][2]], [y, y],
            color="#999999", lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.text(gx, y + 0.42, f"{gap_len / 1000:.1f} kb", ha="center", va="bottom",
            fontsize=6.5, color="#555555")

    lo, hi = (span_hi, span_lo) if reverse else (span_lo, span_hi)
    ax.text(0.0, y - 0.95, f"{lo / 1e6:.3f} Mb", ha="left", va="top", fontsize=6.5)
    ax.text(1.0, y - 0.95, f"{hi / 1e6:.3f} Mb" + (" (rev.)" if reverse else ""),
            ha="right", va="top", fontsize=6.5)


def main() -> None:
    e1, m7 = read_gff(L5_GFF), read_gff(M7_GFF)
    fig, ax = plt.subplots(figsize=(7.2, 2.9))

    for y, name, genes, ann, rev in ((2.3, "E1", E1_GENES, e1, False),
                                     (0.0, "M7", M7_GENES, m7, True)):
        boxes, gap, (c_lo, c_hi), first = layout(genes, ann, rev)
        lo = min(first[0], c_lo)
        hi = max(first[1], c_hi)
        if name == "M7":
            lo = min([lo] + [ann[t][0] for t in M7_SPAN_ALSO])
            hi = max([hi] + [ann[t][1] for t in M7_SPAN_ALSO])
        draw_row(ax, y, name, boxes, gap, lo, hi, rev)
        print(f"{name}: {lo:,}-{hi:,} = {hi - lo + 1:,} bp "
              f"({(hi - lo + 1) / 1000:.1f} kb); GH3_e108-galR gap "
              f"{gap:,} bp ({gap / 1000:.1f} kb)")

    ax.text(0.5, 2.9, "β-glucoside-processing island",
            ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xlim(-0.08, 1.04)
    ax.set_ylim(-1.3, 3.9)
    ax.axis("off")
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(f"{OUT}.{ext}", dpi=300 if ext == "png" else None,
                    bbox_inches="tight")
    print(f"written: {OUT}.svg / .png")


if __name__ == "__main__":
    main()
