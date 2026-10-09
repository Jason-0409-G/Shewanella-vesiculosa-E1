"""
Likelihood-ratio tests for the codeml branch and branch-site models.

Two LRTs per OG:
  LRT 1: M0 vs branch model
    H0: one omega across the tree; H1: the E1 branch has its own omega
    2 x (lnL_branch - lnL_M0) ~ chi2(1). These are the P values shown in
    Supplementary Fig. 4.
  LRT 2: branch-site model A vs A-null (omega2 fixed at 1)
    H0: no positive selection on the E1 branch; H1: omega2 > 1 at some sites.
    The statistic is referred to chi2(1), the conservative option; the manuscript
    uses the standard 50:50 mixture of chi2(0) and chi2(1). The choice makes no
    difference here. For bglB, GH3_e108 and the MFS transporter the statistic is
    exactly zero, the maximum-likelihood estimate of omega2 having settled on the
    boundary at 1; for the seven-taxon GH3_e227 run read by this script it is 0.19,
    and for the eight-taxon re-run behind the reported omega it is zero. The mixture
    rejects at 2.71, so none of them is close. These P values are not reported in
    the manuscript.

This script reads the seven-taxon GH3_e227 results, so the omega it prints for that
gene is 0.0717, not the 0.0759 reported in the manuscript, which comes from the
eight-taxon re-run described in README.md. The other three genes match the paper.

OG0001992_GH3_NagZ was part of the original run of this list; the gene lies outside
the island and was not carried into the manuscript, so it is omitted.

References: Yang & Reis 2011 MBE 28:1217; Zhang, Nielsen & Yang 2005 MBE 22:2472.
"""
from __future__ import annotations
import re
from pathlib import Path
from scipy.stats import chi2

PAML_DIR = Path("04_selection_pressure/paml/results")

OGS = [
    "OG0001451_GH1",
    "OG0001452_GH3_e227",
    "OG0002997_GH3_e108",
]

MODELS = ["M0", "branch", "bsA", "bsA_null"]


def extract_lnL(out_file: Path) -> float | None:
    if not out_file.exists():
        return None
    text = out_file.read_text()
    m = re.search(r"lnL\(.*?\):\s+(-?\d+\.\d+)", text)
    return float(m.group(1)) if m else None


def extract_omega_M0(out_file: Path) -> float | None:
    """Extract single ω from M0 output."""
    if not out_file.exists():
        return None
    text = out_file.read_text()
    m = re.search(r"omega \(dN/dS\)\s*=\s+(-?[\d.]+)", text)
    return float(m.group(1)) if m else None


def extract_omegas_branch(out_file: Path) -> tuple[float | None, float | None]:
    """
    Background and foreground (#1) omega from branch-model output, i.e. the line
      w (dN/dS) for branches:  0.07289 0.34521
    """
    if not out_file.exists():
        return (None, None)
    text = out_file.read_text()
    m = re.search(r"w \(dN/dS\) for branches:\s+(.+?)\n", text)
    if not m:
        return (None, None)
    values = [float(x) for x in m.group(1).strip().split()]
    if len(values) >= 2:
        # first value: background (label 0); second: foreground (label 1)
        return (values[0], values[1])
    return (None, None)


def extract_bsA_omegas(out_file: Path) -> dict | None:
    """
    Site classes from branch-site model A output:
      MLEs of dN/dS (w) for site classes:
      site class             0        1       2a       2b
      proportion       0.65329  0.21304  0.10473  0.02894
      background w     0.07289  1.00000  0.07289  1.00000
      foreground w     0.07289  1.00000  9.99999  9.99999
    """
    if not out_file.exists():
        return None
    text = out_file.read_text()
    m = re.search(
        r"site class\s+0\s+1\s+2a\s+2b\s*\n"
        r"\s*proportion\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\n"
        r"\s*background w\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\n"
        r"\s*foreground w\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
        text,
    )
    if not m:
        return None
    g = m.groups()
    return {
        "proportions": [float(g[0]), float(g[1]), float(g[2]), float(g[3])],
        "background_w": [float(g[4]), float(g[5]), float(g[6]), float(g[7])],
        "foreground_w": [float(g[8]), float(g[9]), float(g[10]), float(g[11])],
    }


def lrt_pvalue(lnL_alt: float, lnL_null: float, df: int = 1) -> tuple[float, float]:
    """LRT: 2 x (lnL_alt - lnL_null) ~ chi2(df). Returns (statistic, p-value)."""
    if lnL_alt is None or lnL_null is None:
        return (None, None)
    stat = 2 * (lnL_alt - lnL_null)
    if stat < 0:
        stat = 0  # null fitted as well or better
    pval = 1 - chi2.cdf(stat, df=df)
    return (stat, pval)


def lrt_branchsite(lnL_alt: float, lnL_null: float) -> tuple[float, float]:
    """Branch-site test against chi2(1); conservative relative to the 50:50 mixture."""
    return lrt_pvalue(lnL_alt, lnL_null, df=1)


def main():
    print("PAML LRT Summary")

    summary = []

    for og in OGS:
        ogdir = PAML_DIR / og
        out_files = {
            m: ogdir / f"paml_{m}" / f"{og}_{m}.out"
            for m in MODELS
        }
        lnLs = {m: extract_lnL(f) for m, f in out_files.items()}

        print(f"\n{og}")
        print("  Log-likelihoods:")
        for m in MODELS:
            lnL = lnLs[m]
            print(f"    {m:<10}  lnL = {lnL:.4f}" if lnL is not None else f"    {m:<10}  missing")

        # ----------------- M0 ω -----------------
        omega_M0 = extract_omega_M0(out_files["M0"])
        if omega_M0 is not None:
            print(f"\n  M0 single ω = {omega_M0:.4f}")

        # ----------------- Branch model ω -----------------
        omega_bg, omega_fg = extract_omegas_branch(out_files["branch"])
        if omega_bg is not None:
            print("  Branch model:")
            print(f"    Background omega                = {omega_bg:.4f}")
            print(f"    Foreground omega (E1 lineage)  = {omega_fg:.4f}")
            print(f"    E1/background ratio             = {omega_fg/omega_bg:.2f}x" if omega_bg > 0 else "")

        # ----------------- bsA site classes -----------------
        bsA_omegas = extract_bsA_omegas(out_files["bsA"])
        if bsA_omegas:
            print("\n  Branch-site Model A (E1 #1 foreground):")
            print("    Site class      proportion    bg ω      fg ω")
            classes = ["0 (purifying)", "1 (neutral)", "2a (pos sel BG=pur)", "2b (pos sel BG=neut)"]
            for i, cls in enumerate(classes):
                p = bsA_omegas["proportions"][i]
                bw = bsA_omegas["background_w"][i]
                fw = bsA_omegas["foreground_w"][i]
                print(f"    {cls:<22} {p:<10.4f}  {bw:<8.4f}  {fw:<8.4f}")

        # ----------------- LRT 1: M0 vs branch -----------------
        stat1, pval1 = lrt_pvalue(lnLs["branch"], lnLs["M0"], df=1)
        # ----------------- LRT 2: bsA vs bsA_null -----------------
        stat2, pval2 = lrt_branchsite(lnLs["bsA"], lnLs["bsA_null"])

        def fmt(stat, pval):
            if stat is None or pval is None:
                return "not computed (a codeml output is missing)"
            return f"2\u0394lnL = {stat:>7.3f}   p = {pval:.4f}"

        print("\n  LRT tests")
        print(f"  LRT 1 (M0 vs branch):     {fmt(stat1, pval1)}")
        print(f"  LRT 2 (bsA vs bsA null):  {fmt(stat2, pval2)}")

        print("\n  Interpretation:")
        if pval1 is not None and pval1 < 0.05 and None not in (omega_fg, omega_bg):
            if omega_fg > omega_bg:
                print(f"    E1 omega ({omega_fg:.4f}) significantly higher than background ({omega_bg:.4f})")
                if omega_fg > 1:
                    print("    E1 omega > 1: positive selection at lineage level")
                else:
                    print("    E1 omega < 1: purifying selection, relaxed relative to background")
            else:
                print("    E1 omega significantly lower than background (stronger purifying selection)")
        else:
            print("    No significant omega difference between E1 and background lineages")

        if pval2 is not None and pval2 < 0.05:
            print("    Branch-site test: positive selection at specific sites on the E1 lineage")
        else:
            print("    Branch-site test: no site-specific positive selection on the E1 lineage")

        summary.append({
            "OG": og,
            "lnL_M0": lnLs["M0"],
            "lnL_branch": lnLs["branch"],
            "lnL_bsA": lnLs["bsA"],
            "lnL_bsA_null": lnLs["bsA_null"],
            "omega_M0": omega_M0,
            "omega_bg_branch": omega_bg,
            "omega_fg_branch": omega_fg,
            "LRT1_stat": stat1,
            "LRT1_pval": pval1,
            "LRT2_stat": stat2,
            "LRT2_pval": pval2,
        })

    import csv
    out_tsv = PAML_DIR.parent / "paml_LRT_summary.tsv"
    with open(out_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    print(f"\nSummary table saved: {out_tsv}")


if __name__ == "__main__":
    main()
