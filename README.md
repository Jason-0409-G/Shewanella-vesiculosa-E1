# Krill gut *Shewanella vesiculosa* E1 — analysis code

Code accompanying the manuscript:

> **Krill gut *Shewanella vesiculosa* supports cellobiose utilization with implications for Southern Ocean carbon cycling**
> Jian Gao, HongZhi Tang, Sheng Du, PengFei Zheng, QunJian Yin, MengYu Liu, FuJian Peng, ChangWei Shao\*, Liang Meng\*, ZhanFei Wei\*
> Submitted to *Communications Earth & Environment*

---

## Data availability

| Record | Accession |
|---|---|
| BioProject | [PRJNA1478518](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1478518) |
| BioSample | SAMN60895750 |
| Genome assembly (WGS) | JBZUHQ000000000 — single closed circular chromosome, 4,858,980 bp, locus-tag prefix `AC4LHR` |
| Raw PacBio HiFi reads | SRA submission SUB16537385 |
| Krill-gut metatranscriptomes | PRJNA1210771 (previously published) |

Genome and read records are held until publication.

---

## Repository layout

Folders follow the order of the Methods section.

| Folder | Contents |
|---|---|
| `01_assembly_qc/` | Read filtering, assembly and CheckM2 assessment of the E1 chromosome, and the CheckM completeness screen of the 141-genome panel |
| `02_annotation/` | Prokka, dbCAN3 and eggNOG-mapper runners for the eight focal genomes, and the blastp assignment of `bglT` |
| `03_comparative_genomics/` | Orthogroup classification, genome-feature statistics, functional-module scans, the *S. vesiculosa* shared-gene scan, verification of island orthologs in M7, whole-genome synteny, ANI and AAI |
| `04_phylogenomics/` | Genus-wide tree of 141 genomes (Fig. 2a) and species tree from 100 random single-copy orthogroups (Fig. 2b) |
| `05_selection_pressure/` | PAML pipeline (**see the folder's own README — the scripts are not a single consistent pass**), numbered by execution order: sequence extraction → MAFFT L-INS-i → PAL2NAL → IQ-TREE gene trees → foreground labelling → codeml (four models) → likelihood-ratio tests. `04b`/`06b` add the MFS transporter, `00b`/`00c` prepare two inputs that the numbered steps assume already exist |
| `06_cazy_kegg_transporters/` | CAZyme consensus matrix, carbohydrate-substrate gene tallies, transporter-family counts, symbiosis-gene mining, and the KEGG module-completeness run and table behind Fig. 4a and Supplementary Data 4 |
| `07_metatranscriptome/` | BWA-MEM alignment to the E1 coding sequences, CoverM TPM quantification, eggNOG annotation merge |
| `08_carbon_budget/` | Deterministic first-order carbon budget (v3.2) and E1-share sensitivity analysis |
| `09_figures/` | Figure generators, named by their figure number in the manuscript. `Fig2a_itol_annotations.py` writes the annotation files that Fig. 2a was rendered with on iTOL |

### Carbon budget

`08_carbon_budget/carbon_load_v32_final.py` implements the deterministic budget of Supplementary Note 3 (`v32` denotes model version 3.2):

```
C_ing = (Q/PP) × NPP
F_hyd = C_ing × φ × ε          where ε = η × h
G     = K₂ × F_hyd             R = (1 − K₂) × F_hyd,  K₂ = GGE/AE
p_E1  = s_E1·m / (s_E1·m + f)  F_E1 = p_E1 × F_hyd
```

Standard library only; re-execution returns identical values. Earlier Monte Carlo formulations were superseded and are not included.

---

## Software

| Step | Tool |
|---|---|
| Read filtering | Filtlong v0.3.1 |
| Assembly | Flye v2.9.6-b1802 (`--pacbio-hifi`) |
| Assembly QC | CheckM2 v1.1.0, CheckM v1.2 |
| Annotation | Prokka v1.15.6, eggNOG-mapper v2.1.13 (eggNOG DB v5.0.2), dbCAN3 v5.2.8 |
| Orthology | OrthoFinder v3.1.2 |
| Alignment | MAFFT v7.526, PAL2NAL v14, trimAl v1.4.1 |
| Phylogeny | IQ-TREE v3.1.1 |
| Selection | PAML codeml v4.10.10 |
| ANI / AAI | FastANI v1.34, EzAAI v1.2.4 |
| Synteny | MUMmer4 v4.0.1 |
| Read QC / mapping | fastp v0.22.0, BWA-MEM v0.7.17, CoverM v0.6.1 |
| Pathway completeness | kegg-pathways-completeness v1.4.3 |
| Homology search | BLAST+ v2.17.0 |

Full version provenance is in `SOFTWARE_VERSIONS.md`. Python scripts use numpy, pandas, scipy, matplotlib, Pillow, biopython, pycirclize and lxml; R scripts use ComplexHeatmap, circlize and svglite. SOFTWARE_VERSIONS.md lists these together with the external tools.

---

## Scope and limitations

Please read this section before attempting to re-run anything.

**These scripts are a record of how the reported results were produced, not a turnkey pipeline.** They were written against the original project directory layout and have been collected here for publication.

**Steps without code in this repository.** Some analyses in the Methods were run interactively or on a cluster and no script was retained. For these, the Methods section gives the tool, version and parameters, but this repository provides no runnable code:

- the OrthoFinder run itself (several scripts here consume its output). The command was
  `orthofinder -f <input> -t 8 -a 8 -S diamond -n run_v1`
- the rendering of Fig. 2a, which was done on the iTOL web service; the annotation files it was
  given are produced by `09_figures/Fig2a_itol_annotations.py`, which builds them from Supplementary Data 2
- read quality filtering (fastp); the metatranscriptome script starts from clean reads
- the eight-taxon re-run of GH3_e227 that the reported omega of 0.0759 comes from.
  `05_selection_pressure/00c_reassemble_GH3_e227_full_length.py` rebuilds the two corrected
  coding sequences, but the alignment, tree and codeml steps that follow were run interactively
  and no script was kept; that run also used `cleandata = 1`, unlike the other genes
- extraction of the 16S rRNA gene from the assembly
- signal-peptide and localization prediction (SignalP 6.0, PSORTb 3.0), both run on their web servers
- Fig. 3b, a two-ellipse Venn drawn directly in the vector editor. Its three counts
  (264 E1-only / 3,904 shared / 443 M7-only) come from a pairwise OrthoFinder comparison of the
  two S. vesiculosa genomes. That two-genome run is not reproduced here and its counts do not
  fall out of the eight-genome run in `03_comparative_genomics/`, which uses a different
  orthogroup partition; the per-gene result it produced, the 264 E1 locus tags, was retained

**Several panels need inputs that are not redistributed here.** `Fig1b_TEM_crop.py` needs the raw
micrograph; `Fig4a_kegg_module_completeness.py` and `Fig4b_island_tpm_rank.py` need the TPM table
that `07_metatranscriptome/` rebuilds from the metatranscriptome accessions; `Fig2a_itol_annotations.py`
and the tree scripts need the Newick files that `04_phylogenomics/` and `05_selection_pressure/`
produce; `Fig5b_mechanism_panel.py` needs the icon set. Each script names the input it expects.

**Three figure scripts were rewritten after the fact.** The scripts that produced the published
Fig. 1c, Fig. 1d and Fig. 2b were not kept. `Fig1c_circular_genome_map.py`,
`Fig1d_island_structure_E1_M7.py` and `Fig2b_species_tree_8genomes.py` were written afterwards from
the same inputs and reproduce the data content of those panels (feature positions, window
statistics, gene order and lengths, topology, branch lengths and support values), not their
published styling.

**Figure scripts produce the base artwork, not the final panels.** Every panel was afterwards
arranged, relabelled and recoloured in a vector editor. A re-run therefore reproduces the data
content — values, positions, topologies, gene order and lengths — but not the published layout,
labels or palette. The supplementary figures are generated as editable SVG through
`09_figures/_svglib.py`. Two scripts exist for
Supplementary Fig. 3: `SuppFig3_island_gene_trees.py` reads the IQ-TREE `.treefile` outputs, while
`SuppFig3_island_gene_trees_layout.py` inlines the same newick strings and reproduces the published
layout. `Fig5b_mechanism_panel.py` additionally needs the icon set in `data/fig5img/`, which is not
redistributed here.

**Gene identifiers differ from the public assembly.** Scripts use the internal locus-tag prefix `PB002_` (for example `PB002_01811` for GH3_e227). The deposited assembly JBZUHQ000000000 uses the prefix `AC4LHR`. Gene identifiers therefore need to be mapped before these scripts can be run against the public record.

**Strain naming.** The strain appears as `L5` (laboratory), `PB002` (sequencing) and `E1` (publication). The manuscript uses `E1` throughout. Display labels in the figure scripts have been updated to `E1`; internal dictionary keys, data-file names and locus tags retain the earlier identifiers so that they still match the intermediate files they were written against.

**Paths.** Scripts reference data through paths relative to the original project layout, or through `.`, `data/` and environment variables where absolute paths were removed. Adjust them to your own locations. Large inputs — genomes, annotation tables, alignments, read data — are not included here; retrieve them from the accessions listed above.

---

## Citation

Please cite the paper once published. Until then, cite this repository and the BioProject accession.

## License

MIT — see `LICENSE`.

## Contact

Corresponding authors: ChangWei Shao (shaocw@ysfri.ac.cn), Liang Meng (mengliang1@genomics.cn), ZhanFei Wei (weizf@ysfri.ac.cn)
