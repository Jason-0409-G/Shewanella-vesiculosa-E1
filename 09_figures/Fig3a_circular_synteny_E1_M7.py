#!/usr/bin/env python3
"""Fig. 3a - circular whole-genome synteny between S. vesiculosa E1 and M7.

Each MUMmer alignment block of the filtered E1-vs-M7 (reverse-complemented)
comparison is drawn as one thin Bezier link between the block midpoints: blue
for forward blocks, red for inversions. The two chromosomes form a 344-degree
circle with a 16-degree gap at the top.

Input
-----
E1_vs_M7rc.coords   show-coords -T -l -c output (synteny_mummer_E1_M7.sh);
                    127 blocks = 120 forward + 7 inverted

Output
------
Fig3a_circular_synteny_E1_M7.{png,pdf,svg}; the SVG is split into three
Inkscape layers (genome arcs, synteny lines, text and legend).

Requires pycirclize, pandas, numpy, matplotlib and lxml. The published panel
was relabelled in a vector editor.
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from lxml import etree
from pycirclize import Circos

COORDS = Path(sys.argv[1] if len(sys.argv) > 1 else "E1_vs_M7rc.coords")
OUT = Path("Fig3a_circular_synteny_E1_M7")

E1_LEN = 4_858_980
M7_LEN = 4_782_877

M7_EDGE = "#2B4C7E"
M7_FILL = "#EAF0F8"
E1_EDGE = "#5C3317"
E1_FILL = "#F7EDE4"
FWD_COL = "#4D8FBD"   # forward blocks
INV_COL = "#B03A2E"   # inverted blocks
TEXT_MID = "#555555"


def parse_coords(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path, sep="\t", skiprows=4, header=None,
        names=["S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY",
               "LEN_R", "LEN_Q", "COV_R", "COV_Q", "TAG_R", "TAG_Q"],
    )
    df["orientation"] = np.where(df["S2"] < df["E2"], "Forward", "Inverted")
    df["e1_mid"] = ((df[["S1", "E1"]].min(axis=1) + df[["S1", "E1"]].max(axis=1)) / 2).astype(int)
    df["m7_mid"] = ((df[["S2", "E2"]].min(axis=1) + df[["S2", "E2"]].max(axis=1)) / 2).astype(int)
    print(f"{len(df)} blocks: "
          f"{(df.orientation == 'Forward').sum()} forward, "
          f"{(df.orientation == 'Inverted').sum()} inverted")
    return df


def draw(coords: pd.DataFrame) -> None:
    sectors = {"M7": M7_LEN, "E1": E1_LEN}
    circos = Circos(sectors, space=12, start=-172, end=172)   # 344-degree arc, gap at the top

    for sector in circos.sectors:
        is_e1 = sector.name == "E1"
        track = sector.add_track((92, 100))
        edge = E1_EDGE if is_e1 else M7_EDGE
        track.axis(fc=E1_FILL if is_e1 else M7_FILL, ec=edge, lw=1.8)
        track.xticks_by_interval(
            1_000_000,
            label_formatter=lambda v: f"{v / 1e6:.0f} Mb",
            label_size=7,
            label_orientation="vertical",
            line_kws=dict(color=edge, lw=1.0),
            text_kws=dict(color=TEXT_MID),
        )
        track.xticks_by_interval(
            500_000,
            label_formatter=lambda v: "",
            label_size=0,
            line_kws=dict(color=edge, lw=0.45),
        )
        sector.text(r"$\it{S.\ vesiculosa}$ " + sector.name, r=116, size=10,
                    fontweight="bold", color=edge)

    n_fwd = int((coords.orientation == "Forward").sum())
    n_inv = int((coords.orientation == "Inverted").sum())
    for row in coords.itertuples(index=False):
        fwd = row.orientation == "Forward"
        circos.link_line(
            ("M7", row.m7_mid), ("E1", row.e1_mid),
            color=FWD_COL if fwd else INV_COL,
            alpha=0.40 if fwd else 0.85,
            lw=0.7 if fwd else 1.1,
        )

    fig = circos.plotfig(figsize=(9, 8.5))
    fig.axes[0].legend(
        handles=[
            mpatches.Patch(fc=FWD_COL, ec="none", alpha=0.75,
                           label=f"Forward synteny  (n={n_fwd})"),
            mpatches.Patch(fc=INV_COL, ec="none", alpha=0.9,
                           label=f"Inversion  (n={n_inv})"),
        ],
        loc="lower center", bbox_to_anchor=(0.5, -0.04),
        ncol=2, fontsize=8, frameon=False,
    )

    for suffix, kw in (("png", dict(dpi=300)), ("pdf", {}), ("svg", {})):
        fig.savefig(f"{OUT}.{suffix}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig)
    add_inkscape_layers(Path(f"{OUT}.svg"))


def add_inkscape_layers(svg_path: Path) -> None:
    """Split the contents of axes_1 into three Inkscape layers.

    The layers are created inside axes_1, so that clipPath references stay valid.
    Synteny links are recognised by their stroke width (0.7 or 1.1).
    """
    SVG_NS = "http://www.w3.org/2000/svg"
    INK_NS = "http://www.inkscape.org/namespaces/inkscape"
    LINK_WIDTHS = ("0.7", "1.1")

    tree = etree.parse(str(svg_path), etree.XMLParser(remove_blank_text=True))
    axes1 = next(g for g in tree.getroot().iter(f"{{{SVG_NS}}}g") if g.get("id") == "axes_1")

    def stroke_width(el):
        for tok in el.get("style", "").split(";"):
            if "stroke-width" in tok:
                return tok.split(":")[-1].strip()

    arcs, lines, texts = [], [], []
    for child in list(axes1):
        tag = child.tag.replace(f"{{{SVG_NS}}}", "")
        if tag == "text":
            texts.append(child)
        elif tag == "g":
            if any(stroke_width(p) in LINK_WIDTHS for p in child.iter(f"{{{SVG_NS}}}path")):
                lines.append(child)
            elif any(d.tag == f"{{{SVG_NS}}}text" for d in child.iter()):
                texts.append(child)
            else:
                arcs.append(child)
        elif stroke_width(child) in LINK_WIDTHS:
            lines.append(child)
        else:
            arcs.append(child)
        axes1.remove(child)
    print(f"layers: arcs={len(arcs)} links={len(lines)} text={len(texts)}")

    for label, nodes in (("Genome arcs", arcs), ("Synteny lines", lines),
                         ("Text & legend", texts)):
        layer = etree.SubElement(axes1, f"{{{SVG_NS}}}g")
        layer.set(f"{{{INK_NS}}}groupmode", "layer")
        layer.set(f"{{{INK_NS}}}label", label)
        for n in nodes:
            layer.append(n)

    tree.write(str(svg_path), xml_declaration=True, encoding="utf-8", pretty_print=True)


if __name__ == "__main__":
    draw(parse_coords(COORDS))
    print(f"written: {OUT}.png / .pdf / .svg")
