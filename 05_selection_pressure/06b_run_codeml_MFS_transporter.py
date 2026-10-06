"""
PAML selection-pressure analysis for the MFS transporter (OG0000506, PB002_01810).

Performs PAL2NAL back-translation, PAML tree preparation and the four codeml models,
then adds the result row to paml_LRT_summary.tsv. The protein file and the IQ-TREE gene
tree must already exist (04b_mafft_iqtree_MFS_transporter.sh).

Usage (from the repository root):
    python 05_selection_pressure/06b_run_codeml_MFS_transporter.py

Requires: pal2nal.pl, codeml. All runs use cleandata = 0.
"""
from __future__ import annotations
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT     = Path(".")
PROKKA8  = ROOT / "02_annotation/annotation/8strains_prokka_v1156"
PAML_RES = ROOT / "04_selection_pressure/paml/results"
OG       = "OG0000506_MFS"
OGDIR    = PAML_RES / OG

# OG0000506 members (Orthogroups.txt):
#   FRIG_01384 KDH10_01697 LIVI_03053 PB002_01810
#   POLA_02649 PSYC_04061 SP02_02810 SP14_03070
# M7 = KDH10; its copy is intact and kept.
MEMBERS = {
    "L5_PB002_01810":   (PROKKA8 / "L5_PB002_prokka/L5.ffn",            "PB002_01810"),
    "M7_KDH10_01697":   (PROKKA8 / "M7_REF_prokka/M7.ffn",               "KDH10_01697"),
    "frig_FRIG_01384":  (PROKKA8 / "S_frigidimarina_prokka/frig.ffn",    "FRIG_01384"),
    "livi_LIVI_03053":  (PROKKA8 / "S_livingstonensis_prokka/livi.ffn",  "LIVI_03053"),
    "pola_POLA_02649":  (PROKKA8 / "S_polaris_prokka/pola.ffn",          "POLA_02649"),
    "psyc_PSYC_04061":  (PROKKA8 / "S_psychromarinicola_prokka/psyc.ffn","PSYC_04061"),
    "sp02_SP02_02810":  (PROKKA8 / "S_sp002836315_prokka/sp02.ffn",      "SP02_02810"),
    "sp14_SP14_03070":  (PROKKA8 / "S_sp014164505_prokka/sp14.ffn",      "SP14_03070"),
}

# codeml control parameters
_CTL_BASE = dict(
    noisy=3, verbose=1, runmode=0, seqtype=1, CodonFreq=2,
    clock=0, aaDist=0, icode=0, Mgene=0, fix_kappa=0, kappa=2,
    fix_alpha=1, alpha=0, Malpha=0, ncatG=1, getSE=0,
    RateAncestor=0, Small_Diff="5e-7", cleandata=0, method=0,
)

MODELS = {
    "M0": {**_CTL_BASE, "model": 0, "NSsites": 0, "fix_omega": 0, "omega": 0.4,
           "_tree": "unlabeled"},
    "branch": {**_CTL_BASE, "model": 2, "NSsites": 0, "fix_omega": 0, "omega": 0.045,
               "_tree": "L5fg"},
    "bsA": {**_CTL_BASE, "model": 2, "NSsites": 2, "fix_omega": 0, "omega": 1.5,
            "ncatG": 4, "_tree": "L5fg"},
    "bsA_null": {**_CTL_BASE, "model": 2, "NSsites": 2, "fix_omega": 1, "omega": 1.0,
                 "ncatG": 4, "_tree": "L5fg"},
}


def step3_pal2nal() -> bool:
    """Collect the member CDSs into _nucl.ffn and back-translate with PAL2NAL."""
    print("\n[3] CDS extraction and PAL2NAL")

    nucl_ffn = OGDIR / f"{OG}_nucl.ffn"
    codon_phy = OGDIR / f"{OG}_codon.phy"

    lines = []
    for label, (ffn_path, locus_tag) in MEMBERS.items():
        if not ffn_path.exists():
            print(f"  {ffn_path} not found, aborting")
            return False
        found = False
        with open(ffn_path) as fh:
            seq_lines: list[str] = []
            capture = False
            for line in fh:
                if line.startswith(">"):
                    if capture:
                        lines.append(f">{label}\n" + "".join(seq_lines))
                        found = True
                    capture = locus_tag in line
                    seq_lines = []
                elif capture:
                    seq_lines.append(line)
            if capture and seq_lines:
                lines.append(f">{label}\n" + "".join(seq_lines))
                found = True
        if not found:
            print(f"  {locus_tag} not found in {ffn_path}, aborting")
            return False
        print(f"  {label} ({locus_tag})")

    nucl_ffn.write_text("\n".join(lines) + "\n")
    print(f"  [saved] {nucl_ffn}  ({len(MEMBERS)} sequences)")

    prot_aln = OGDIR / f"{OG}_prot.aln"
    if not prot_aln.exists():
        print(f"  protein alignment not found: {prot_aln}")
        return False

    result = subprocess.run(
        ["pal2nal.pl", str(prot_aln), str(nucl_ffn),
         "-output", "paml", "-codontable", "11", "-nogap"],
        capture_output=True, text=True
    )
    if not result.stdout.strip():
        print(f"  PAL2NAL output empty\n  STDERR: {result.stderr[:400]}")
        return False
    codon_phy.write_text(result.stdout)
    print(f"  PHYLIP header: {result.stdout.split(chr(10))[0].strip()}")
    print(f"  [saved] {codon_phy}")
    return True


def step5_prep_trees() -> bool:
    """Write PAML trees (unlabeled and E1-foreground) from the IQ-TREE treefile."""
    print("\n[5] PAML tree files")

    treefile = OGDIR / f"{OG}_tree.treefile"
    if not treefile.exists():
        print(f"  {treefile} not found")
        return False

    # strip support values (')96.7/92:' or ')100:')
    nw = re.sub(r'\)[\d.]+/[\d.]+:', '):', treefile.read_text().strip())
    nw = re.sub(r'\)[\d.]+:', '):', nw)

    # leaves carry the raw locus tag; E1 = PB002_01810
    l5_locus = MEMBERS["L5_PB002_01810"][1]
    m = re.search(rf'({re.escape(l5_locus)}):', nw)
    if not m:
        print(f"  E1 leaf {l5_locus} not found in tree")
        return False
    l5_label = m.group(1)

    n_taxa = nw.count(",") + 1
    print(f"  E1 leaf: {l5_label}; taxa: {n_taxa}")

    def paml_header(tree: str) -> str:
        return f"  {n_taxa}  1\n{tree.strip()}\n"

    unlabeled_out = OGDIR / f"{OG}_unlabeled.nwk"
    unlabeled_out.write_text(paml_header(nw))
    print(f"  [saved] {unlabeled_out}")

    nw_l5fg = nw.replace(f"{l5_label}:", f"{l5_label} #1:", 1)
    l5fg_out = OGDIR / f"{OG}_L5fg.nwk"
    l5fg_out.write_text(paml_header(nw_l5fg))
    print(f"  [saved] {l5fg_out}")
    return True


def _write_ctl(ctl_path: Path, params: dict, seq: str, tree: str, out: str) -> None:
    lines = [f"      seqfile = {seq}", f"     treefile = {tree}",
             f"      outfile = {out}", ""]
    for k, v in params.items():
        if k.startswith("_"):
            continue
        lines.append(f"  {k:>12} = {v}")
    ctl_path.write_text("\n".join(lines) + "\n")


def _extract_lnL(out_file: Path) -> str:
    if not out_file.exists():
        return "?"
    m = re.search(r"lnL\(.*?\):\s+(-?\d+\.\d+)", out_file.read_text())
    return m.group(1) if m else "?"


def step6_codeml() -> bool:
    """Run the four codeml models."""
    print("\n[6] codeml")

    codon_phy = OGDIR / f"{OG}_codon.phy"
    if not codon_phy.exists():
        print(f"  {codon_phy} missing; run step 3 first")
        return False

    for model_name, params in MODELS.items():
        tree_label = params["_tree"]
        tree_nwk = OGDIR / f"{OG}_{tree_label}.nwk"
        if not tree_nwk.exists():
            print(f"  tree file {tree_nwk} missing, skipping {model_name}")
            continue

        model_dir = OGDIR / f"paml_{model_name}"
        model_dir.mkdir(exist_ok=True)

        seq_rel  = f"../{OG}_codon.phy"
        tree_rel = f"../{OG}_{tree_label}.nwk"
        out_name = f"{OG}_{model_name}.out"

        ctl_path = model_dir / "codeml.ctl"
        _write_ctl(ctl_path, params, seq_rel, tree_rel, out_name)

        print(f"  {model_name} ...", end=" ", flush=True)
        result = subprocess.run(
            ["codeml", "codeml.ctl"],
            cwd=str(model_dir),
            capture_output=True, text=True, timeout=3600,
        )
        out_file = model_dir / out_name
        if "Time used" in result.stdout and out_file.exists():
            print(f"lnL = {_extract_lnL(out_file)}")
        else:
            print("may not have finished")
            if result.stderr:
                print(f"     STDERR: {result.stderr[:300]}")
    return True


def step7_update_summary() -> None:
    """Add the MFS row to paml_LRT_summary.tsv (columns as in 07_likelihood_ratio_tests.py)."""
    print("\n[7] LRT summary")
    import math

    summary_file = ROOT / "04_selection_pressure/paml/paml_LRT_summary.tsv"

    def lnL(model: str) -> float:
        out = OGDIR / f"paml_{model}" / f"{OG}_{model}.out"
        v = _extract_lnL(out)
        return float(v) if v != "?" else float("nan")

    lnL_M0     = lnL("M0")
    lnL_branch = lnL("branch")
    lnL_bsA    = lnL("bsA")
    lnL_bsA_n  = lnL("bsA_null")

    def get_omega_M0() -> float:
        out_file = OGDIR / "paml_M0" / f"{OG}_M0.out"
        if not out_file.exists():
            return float("nan")
        m = re.search(r"omega \(dN/dS\)\s*=\s*([\d.]+)", out_file.read_text())
        return float(m.group(1)) if m else float("nan")

    omega_M0 = get_omega_M0()

    # codeml prints "w (dN/dS) for branches:  <background> <foreground (#1)>"
    def get_branch_omegas() -> tuple[float, float]:
        out_file = OGDIR / "paml_branch" / f"{OG}_branch.out"
        if not out_file.exists():
            return float("nan"), float("nan")
        m = re.search(r"w \(dN/dS\) for branches:\s+([\d.]+)\s+([\d.]+)",
                      out_file.read_text())
        if not m:
            return float("nan"), float("nan")
        return float(m.group(1)), float(m.group(2))

    omega_bg, omega_fg = get_branch_omegas()

    # LRT 1: M0 vs branch (df = 1); LRT 2: bsA vs bsA_null (df = 1)
    lrt1 = 2 * (lnL_branch - lnL_M0) if not math.isnan(lnL_branch) else float("nan")
    lrt2 = 2 * (lnL_bsA - lnL_bsA_n) if not math.isnan(lnL_bsA) else float("nan")

    def chi2_p(stat: float, df: int = 1) -> float:
        if math.isnan(stat) or stat < 0:
            return 1.0
        from scipy.stats import chi2
        return chi2.sf(stat, df)

    p1 = chi2_p(lrt1)
    p2 = chi2_p(lrt2)

    row = (f"OG0000506_MFS\t{lnL_M0:.6f}\t{lnL_branch:.6f}\t"
           f"{lnL_bsA:.6f}\t{lnL_bsA_n:.6f}\t"
           f"{omega_M0:.5f}\t{omega_bg:.5g}\t{omega_fg:.5g}\t"
           f"{lrt1:.6f}\t{p1:.10f}\t{lrt2:.6f}\t{p2:.10f}")

    print(f"  {row}")

    # replace an existing MFS row, otherwise append
    lines = summary_file.read_text().splitlines() if summary_file.exists() else []
    replaced = False
    for i, ln in enumerate(lines):
        if ln.startswith("OG0000506_MFS\t") or ln.strip() == "OG0000506_MFS":
            lines[i] = row
            replaced = True
            break
    if not replaced:
        lines.append(row)
    summary_file.write_text("\n".join(lines) + "\n")
    print(f"  [{'replaced' if replaced else 'appended'}] {summary_file}")


def main() -> None:
    for tool in ["pal2nal.pl", "codeml"]:
        if not shutil.which(tool):
            print(f"{tool} not found in PATH")
            sys.exit(1)

    print(f"OG: {OG}\ndirectory: {OGDIR}\nmembers: {len(MEMBERS)} strains (including M7/KDH10)")

    if not step3_pal2nal():
        sys.exit(1)
    if not step5_prep_trees():
        sys.exit(1)

    step6_codeml()
    step7_update_summary()
    print("\nDone.")


if __name__ == "__main__":
    main()
