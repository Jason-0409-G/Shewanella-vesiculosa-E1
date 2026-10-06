#!/usr/bin/env python3
"""Assemble Supplementary Data 4 (KEGG module completeness, 9 organisms).

Inputs
------
``<out_dir>/{org}/{org}_pathways.tsv`` for org in host, E1, M7, frig, livi, pola,
psyc, sp002, sp014, as produced by ``run_kegg_module_completeness.sh``
(kegg-pathways-completeness v1.4.3, ``give_completeness``).  Each file lists only
the modules with completeness > 0 in that organism, with columns
``module_accession, completeness, pathway_name, pathway_class, matching_ko, missing_ko``.

Outputs
-------
1. ``<res_dir>/TableS3_module_completeness_allmodules.tsv``: Supplementary Data 4
   (sheet 4 of the submitted workbook).  Columns: Module ID, Module name, KEGG class
   (leading "Pathway modules; " removed), then completeness (%) for Krill (= host),
   S. vesiculosa E1, S. vesiculosa M7, S. frigidimarina, S. livingstonensis,
   S. polaris, S. psychromarinicola, S. sp002836315, S. sp014164505.  One row per module
   present (> 0 %) in at least one organism (331 modules); modules absent from an
   organism are reported as 0.0.  Rows are sorted by KEGG class, then module ID.
   Values are rounded to one decimal via ``float(f"{x:.1f}")`` (the tool prints two decimals).
2. ``<res_dir>/module_matrix_9org.tsv``: the same data in wide form with short organism
   keys, the full KEGG class string and the module name (the matrix used by Fig. 4a).

Usage
-----
    python build_module_completeness_table.py [out_dir] [res_dir]
    # defaults: ./data/modcomp/out  and  ./results
"""
import sys
from pathlib import Path

import pandas as pd

ORGS = ["host", "E1", "M7", "frig", "livi", "pola", "psyc", "sp002", "sp014"]
HEADERS = ["Krill %", "S. vesiculosa E1 (%)", "S. vesiculosa M7 (%)",
           "S. frigidimarina (%)", "S. livingstonensis (%)", "S. polaris (%)",
           "S. psychromarinicola (%)", "S. sp002836315 (%)", "S. sp014164505 (%)"]
PREFIX = "Pathway modules; "


def main(out_dir: Path, res_dir: Path) -> None:
    comp, meta = {}, {}
    for org in ORGS:
        df = pd.read_csv(out_dir / org / f"{org}_pathways.tsv", sep="\t")
        comp[org] = df.set_index("module_accession")["completeness"]
        for r in df.itertuples():
            meta.setdefault(r.module_accession, (r.pathway_name, r.pathway_class))

    mat = pd.DataFrame(comp).reindex(columns=ORGS).fillna(0.0)
    mat = mat.apply(lambda c: c.map(lambda x: float(f"{x:.1f}")))
    mat["class"] = [meta[m][1] for m in mat.index]
    mat["name"] = [meta[m][0] for m in mat.index]
    mat.index.name = "module"
    res_dir.mkdir(parents=True, exist_ok=True)
    mat.sort_index().reset_index().to_csv(res_dir / "module_matrix_9org.tsv",
                                          sep="\t", index=False)

    tab = pd.DataFrame({
        "Module ID": mat.index,
        "Module name": mat["name"].values,
        "KEGG class": mat["class"].str.replace(PREFIX, "", regex=False).values,
    })
    for h, o in zip(HEADERS, ORGS):
        tab[h] = mat[o].values
    tab = tab.sort_values(["KEGG class", "Module ID"], kind="stable")
    tab.to_csv(res_dir / "TableS3_module_completeness_allmodules.tsv",
               sep="\t", index=False)
    print(f"{len(tab)} modules x {len(ORGS)} organisms written to {res_dir}")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "data/modcomp/out"),
         Path(sys.argv[2] if len(sys.argv) > 2 else "results"))
