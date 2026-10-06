"""
Run PAML codeml: four models for each OG in OGS.

Models:
  M0        one-ratio model (model = 0, NSsites = 0)
  branch    two-ratio branch model (model = 2, NSsites = 0); E1 branch marked #1
  bsA       branch-site model A (model = 2, NSsites = 2); E1 branch marked #1
  bsA_null  branch-site model A with omega2 fixed at 1

Likelihood-ratio tests are computed in 07_likelihood_ratio_tests.py.
All runs use cleandata = 0.

OG0001992_GH3_NagZ was part of the original run of this list; the gene lies outside
the island and was not carried into the manuscript, so it is omitted. The reported
GH3_e227 omega comes from a separate eight-taxon re-run with cleandata = 1 (README.md).

References: Yang 2007 MBE 24:1586 (PAML 4); Zhang, Nielsen & Yang 2005 MBE 22:2472.
"""
from __future__ import annotations
import re
import shutil
import subprocess
from pathlib import Path

PAML_DIR = Path("04_selection_pressure/paml/results")

OGS = [
    "OG0001451_GH1",
    "OG0001452_GH3_e227",
    "OG0002997_GH3_e108",
]


# codeml control file templates
CTL_TEMPLATE_M0 = """
      seqfile = {seq}
     treefile = {tree}
      outfile = {out}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 0
      NSsites = 0
        icode = 0
        Mgene = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 1
        getSE = 0
 RateAncestor = 0
   Small_Diff = 5e-7
    cleandata = 0
       method = 0
"""

CTL_TEMPLATE_BRANCH = """
      seqfile = {seq}
     treefile = {tree}
      outfile = {out}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 2
      NSsites = 0
        icode = 0
        Mgene = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 1
        getSE = 0
 RateAncestor = 0
   Small_Diff = 5e-7
    cleandata = 0
       method = 0
"""

CTL_TEMPLATE_BSA = """
      seqfile = {seq}
     treefile = {tree}
      outfile = {out}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 2
      NSsites = 2
        icode = 0
        Mgene = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 1.5
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 4
        getSE = 0
 RateAncestor = 0
   Small_Diff = 5e-7
    cleandata = 0
       method = 0
"""

CTL_TEMPLATE_BSA_NULL = """
      seqfile = {seq}
     treefile = {tree}
      outfile = {out}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 2
      NSsites = 2
        icode = 0
        Mgene = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 1
        omega = 1.0
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 4
        getSE = 0
 RateAncestor = 0
   Small_Diff = 5e-7
    cleandata = 0
       method = 0
"""

MODELS = {
    "M0":       (CTL_TEMPLATE_M0,       "unlabeled"),
    "branch":   (CTL_TEMPLATE_BRANCH,   "L5fg"),
    "bsA":      (CTL_TEMPLATE_BSA,      "L5fg"),
    "bsA_null": (CTL_TEMPLATE_BSA_NULL, "L5fg"),
}


def extract_lnL(out_file: Path) -> str:
    if not out_file.exists():
        return "?"
    m = re.search(r"lnL\(.*?\):\s+(-?\d+\.\d+)", out_file.read_text())
    return m.group(1) if m else "?"


def run_codeml_for_og(og: str) -> None:
    ogdir = PAML_DIR / og

    if not ogdir.exists():
        print(f"  Directory {ogdir} does not exist, skipping")
        return

    for model_name, (template, tree_label) in MODELS.items():
        print(f"\n--- {og} :: {model_name} ---")

        model_dir = ogdir / f"paml_{model_name}"
        model_dir.mkdir(exist_ok=True)

        # paths relative to model_dir, where codeml is run
        seq_path = f"../{og}_codon.phy"
        tree_path = f"../{og}_{tree_label}.nwk"
        out_path = f"{og}_{model_name}.out"
        ctl_path = model_dir / "codeml.ctl"

        seq_abs = ogdir / f"{og}_codon.phy"
        tree_abs = ogdir / f"{og}_{tree_label}.nwk"
        if not seq_abs.exists():
            print(f"  Missing seqfile: {seq_abs}")
            continue
        if not tree_abs.exists():
            print(f"  Missing treefile: {tree_abs}")
            continue

        ctl_content = template.format(
            seq=seq_path,
            tree=tree_path,
            out=out_path,
        ).strip()
        ctl_path.write_text(ctl_content + "\n")

        try:
            result = subprocess.run(
                ["codeml", "codeml.ctl"],
                cwd=str(model_dir),
                capture_output=True,
                text=True,
                timeout=3600,
            )
            if result.returncode != 0:
                print(f"  codeml exit {result.returncode}")
                if result.stderr:
                    print(f"  STDERR (last 500): {result.stderr[-500:]}")
            elif "Time used" in result.stdout:
                print(f"  Done. lnL = {extract_lnL(model_dir / out_path)}")
            else:
                print("  codeml may not have finished cleanly")
                print(f"  STDOUT (last 300): {result.stdout[-300:]}")
        except subprocess.TimeoutExpired:
            print(f"  codeml timed out after 1 hour for {og} {model_name}")


def main():
    codeml_path = shutil.which("codeml")
    if not codeml_path:
        raise SystemExit("codeml not found in PATH")
    print(f"codeml: {codeml_path}")
    print(f"{len(OGS)} OGs x {len(MODELS)} models ({', '.join(MODELS)}) = {len(OGS) * len(MODELS)} runs")

    for og in OGS:
        print(f"\nOG: {og}")
        run_codeml_for_og(og)

    print("\nAll codeml runs done. Next: 07_likelihood_ratio_tests.py")


if __name__ == "__main__":
    main()
