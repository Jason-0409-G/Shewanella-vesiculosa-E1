#!/usr/bin/env python3
"""iTOL annotation files for the genus-wide tree of Fig. 2a.

Fig. 2a was rendered on the iTOL web service from the 141-tip tree plus the
annotation files written here. Every value is taken from Supplementary Data 2,
so the panel can be rebuilt from the published table alone.

Tracks, in the order they appear around the tree:
  iTOL_01_focal_highlight  the focal strain in red, bold italic
  iTOL_02_genus_strip   clade colours, by the genus prefix of the leaf label
  iTOL_04_labels        leaf labels replaced by the published species names
  iTOL_05_completeness  CheckM v1.2 completeness, as a gradient
  iTOL_06_GC_gradient   GC content, as a gradient
  iTOL_07_genome_size   genome size in Mb, as a bar
Bootstrap symbols and the tree scale are drawn by iTOL itself from the newick
file and are not configured here.

The tree keeps the internal leaf name for the focal strain (PB002_S_vesiculosa_L5);
the labels file is what renames it to S. vesiculosa E1 for display.

Input   Supplementary Data.xlsx, worksheet "Supplementary Data 2"
Output  iTOL_*.txt in the current directory, to be dragged onto the iTOL tree
"""
import sys
import zipfile
from pathlib import Path

from lxml import etree

XLSX = Path(sys.argv[1] if len(sys.argv) > 1 else "Supplementary Data.xlsx")
SHEET = "xl/worksheets/sheet2.xml"
OUT = Path(".")

FOCAL_LEAF = "PB002_S_vesiculosa_L5"
FOCAL_NAME = "S. vesiculosa E1"
FOCAL_COLOUR = "#B6392B"

GENUS_COLOURS = {
    "Shewanella":     "#C04E2C",
    "Parashewanella": "#BEA06B",
    "Ferrimonas":     "#266C9A",
    "Paraferrimonas": "#8B5E3C",
    "SM1919":         "#1C8458",
    "UBA2521":        "#763A8E",
}


def read_supp_data_2(path: Path) -> list[dict]:
    z = zipfile.ZipFile(path)
    shared = ["".join(t.text or "" for t in si.iter() if t.tag.endswith("}t"))
              for si in etree.fromstring(z.read("xl/sharedStrings.xml"))]
    root = etree.fromstring(z.read(SHEET))
    ns = root.tag.split("}")[0] + "}"
    rows = []
    for row in root.iter(ns + "row"):
        cells = []
        for c in row.iter(ns + "c"):
            v = c.find(ns + "v")
            cells.append(shared[int(v.text)] if (c.get("t") == "s" and v is not None)
                         else (v.text if v is not None else ""))
        rows.append(cells)
    header = next(r for r in rows if "Tree leaf label" in r)
    idx = {name: i for i, name in enumerate(header)}
    out = []
    for r in rows[rows.index(header) + 1:]:
        if len(r) <= idx["Tree leaf label"] or not r[idx["Tree leaf label"]]:
            continue
        out.append({k: r[i] for k, i in idx.items() if i < len(r)})
    return out


def leaf_of(row: dict) -> str:
    """Leaf name as it appears in the newick file."""
    label = row["Tree leaf label"]
    return FOCAL_LEAF if label == "Shewanella_vesiculosa_E1" else label


def genus_of(leaf: str) -> str:
    return "Shewanella" if leaf == FOCAL_LEAF else leaf.split("_", 1)[0]


def write(name: str, header: list[str], data: list[str]) -> None:
    (OUT / name).write_text("\n".join(header + ["DATA"] + data) + "\n")
    print(f"  {name}: {len(data)} rows")


def main() -> None:
    rows = read_supp_data_2(XLSX)
    print(f"{XLSX.name}: {len(rows)} genomes")

    genera = sorted({genus_of(leaf_of(r)) for r in rows})
    print(f"  genera: {', '.join(genera)}")

    write("iTOL_02_genus_strip.txt",
          ["DATASET_COLORSTRIP", "SEPARATOR TAB", "DATASET_LABEL\tGenus",
           "COLOR\t#888888", "STRIP_WIDTH\t28", "MARGIN\t4", "BORDER_WIDTH\t0",
           "LEGEND_TITLE\tGenus",
           "LEGEND_SHAPES\t" + "\t".join("1" for _ in genera),
           "LEGEND_COLORS\t" + "\t".join(GENUS_COLOURS.get(g, "#CCCCCC") for g in genera),
           "LEGEND_LABELS\t" + "\t".join(genera)],
          [f"{leaf_of(r)}\t{GENUS_COLOURS.get(genus_of(leaf_of(r)), '#CCCCCC')}"
           for r in rows])

    write("iTOL_04_labels.txt",
          ["LABELS", "SEPARATOR TAB"],
          [f"{leaf_of(r)}\t{FOCAL_NAME if leaf_of(r) == FOCAL_LEAF else r['Species']}"
           for r in rows])

    for fname, label, column, colour in (
            ("iTOL_05_completeness.txt", "CheckM completeness (%)", "Completeness (%)", "#1C8458"),
            ("iTOL_06_GC_gradient.txt",  "GC (%)",                  "GC (%)",           "#266C9A")):
        vals = [(leaf_of(r), r[column]) for r in rows if r.get(column)]
        lo = min(float(v) for _, v in vals)
        write(fname,
              ["DATASET_GRADIENT", "SEPARATOR TAB", f"DATASET_LABEL\t{label}",
               f"COLOR\t{colour}", "STRIP_WIDTH\t28", "MARGIN\t4", "BORDER_WIDTH\t0",
               "COLOR_MIN\t#F5EEE8", f"COLOR_MAX\t{colour}",
               f"LEGEND_TITLE\t{label}", "LEGEND_SHAPES\t1\t1",
               f"LEGEND_COLORS\t#F5EEE8\t{colour}",
               f"LEGEND_LABELS\t{lo:.1f}\t100"],
              [f"{leaf}\t{float(v):.2f}" for leaf, v in vals])

    sizes = [(leaf_of(r), float(r["Genome size (bp)"]) / 1e6)
             for r in rows if r.get("Genome size (bp)")]
    write("iTOL_07_genome_size.txt",
          ["DATASET_SIMPLEBAR", "SEPARATOR TAB", "DATASET_LABEL\tGenome size (Mb)",
           "COLOR\t#BEA06B", "WIDTH\t120", "MARGIN\t6", "BORDER_WIDTH\t0",
           "DATASET_SCALE\t2\t4\t6"],
          [f"{leaf}\t{mb:.3f}" for leaf, mb in sizes])

    write("iTOL_01_focal_highlight.txt",
          ["TREE_COLORS", "SEPARATOR TAB"],
          [f"{FOCAL_LEAF}\tlabel\t{FOCAL_COLOUR}\tbold-italic",
           f"{FOCAL_LEAF}\tbranch\t{FOCAL_COLOUR}\tnormal\t3"])

    print(f"\nFocal strain {FOCAL_LEAF} -> {FOCAL_NAME}, highlighted in {FOCAL_COLOUR}")


if __name__ == "__main__":
    main()
