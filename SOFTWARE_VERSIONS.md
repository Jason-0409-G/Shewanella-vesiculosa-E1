# Software and database versions

Versions of every tool used for the analyses reported in the paper, and of the libraries the
scripts in this repository depend on. Each version was read from a run log, a tool header or a
`--version` call rather than recalled, and the column "verified from" says which.

Paths such as `02_annotation/annotation/L5_prokka/L5.log` refer to the original project
directory, not to this repository; the data files themselves are not redistributed here and
should be retrieved from the accessions in the Data availability statement.

---

## 1. Assembly and annotation

| Tool | Version | Database | Verified from |
|---|---|---|---|
| Filtlong | 0.3.1 | — | `filtlong --version` in the assembly notebook. Run as `--min_length 1000 --min_mean_q 90`; it removed two reads |
| Flye | 2.9.6-b1802 | — | the Flye log of the E1 run, `INFO: Starting Flye 2.9.6-b1802`. Its `assembly_info.txt` records one closed circular contig, 4,858,980 bp at 133x. That log is not in the current project tree; the assembly summary is |
| CheckM2 | 1.1.0 | CheckM2 database, Zenodo record 14897628 | the CheckM2 run log, `INFO: Running CheckM2 version 1.1.0`. The surviving `quality_report.tsv` gives E1 as 100 % complete with 0.51 % contamination |
| CheckM | 1.2 | `checkm_data_2015_01_16` | used for the 141-genome panel (`lineage_wf --genes -x faa`, `qa -o 2 --tab_table`). These are the values in the Completeness column of Supplementary Data 2 and on Fig. 2a. CheckM and CheckM2 use different marker sets and their values are not comparable |
| Prokka | 1.15.6 | bundled UniProt and Pfam HMMs | `L5.log` line 1 |
| dbCAN3 | 5.2.8 | dbCAN-HMM, CAZy.dmnd, dbCAN-sub, all at the same release | batch log. DIAMOND step `--evalue 1e-102 --max-target-seqs 1`; HMM step `evalue 1e-15 cov 0.35` |
| eggNOG-mapper | 2.1.13-8d129f4 | eggNOG 5.0.2 | `emapper.py --version`, and the header of every `.emapper.annotations` file |
| DIAMOND | 2.1.24 | — | reported by `emapper.py --version` as the eggNOG-mapper back end |
| BLAST+ | 2.17.0+ | — | Prokka dependency log. Used for the blastp assignment of `bglT`; that comparison was re-run and reproduced with 2.16.0+, which gives the same percentage identities |

## 2. Phylogenetics

| Tool | Version | Verified from |
|---|---|---|
| IQ-TREE | 3.1.1, built 8 Apr 2026 | line 1 of every `*.log` |
| ModelFinder | bundled with IQ-TREE | `Best-fit model according to BIC` in each `.iqtree` report. `-m MFP` searched 1,232 protein models for the supermatrix trees; the single-gene trees used `-m TEST`, a smaller set |
| UFBoot2 | bundled with IQ-TREE | `-B 1000` / `-bb 1000` in the logs |
| SH-aLRT | bundled with IQ-TREE | `-alrt 1000` in the logs. Note that the 141-genome tree was run without it |
| MAFFT | 7.526 (2024-04-26) | `mafft --version` |
| trimAl | 1.4.1 | tool header. Used for the eight-strain supermatrix only, not for the 141-genome tree |

## 3. Orthology and comparative genomics

| Tool | Version | Verified from |
|---|---|---|
| OrthoFinder | 3.1.2 | `Log.txt` line 1 of the eight-genome run |
| FastANI | 1.34 | `fastANI --version`. Run with k = 16 and 3,000-bp fragments, the defaults for this version |
| EzAAI | 1.2.4 (Jul 2025) | `ezaai_calculate.log`. Run with the default identity and coverage thresholds, 0.4 and 0.5 |
| MUMmer4 | 4.0.1 | `nucmer --version`. The `dnadiff` wrapper reports 1.3 internally |

## 4. Selection pressure

| Tool | Version | Verified from |
|---|---|---|
| PAML codeml | 4.10.10 (29 Jan 2026) | the `CODONML (in paml version 4.10.10, 29 Jan 2026)` banner in each `*.out` file |
| PAL2NAL | 14 | `pal2nal.pl` header |

## 5. Metatranscriptome and metabolic reconstruction

| Tool | Version | Verified from |
|---|---|---|
| fastp | 0.22.0 | the quality-filtering step of the source study; the reads re-analysed here were obtained already filtered |
| BWA-MEM | 0.7.17 | pipeline script. Run as `bwa mem -M -R <RG>` |
| samtools | 1.15 | pipeline script, for sort, index and flagstat |
| CoverM | 0.6.1 | `coverm --version`. Run as `-m trimmed_mean tpm --min-covered-fraction 0` |
| kegg-pathways-completeness | 1.4.3 | `pip show kegg-pathways-completeness`; entry point `give_completeness`, using the module definitions bundled with that release |

## 6. Protein feature prediction

| Tool | Version | Note |
|---|---|---|
| SignalP | 6.0 | signal-peptide prediction for the island proteins; run on the web server, no script in this repository |
| PSORTb | 3.0 | subcellular localization; likewise run on the web server |

## 7. Libraries the scripts in this repository need

| Language | Packages |
|---|---|
| Python 3 | biopython, lxml, matplotlib, numpy, pandas, Pillow, pycirclize, scipy |
| R | ComplexHeatmap, circlize, svglite |

`09_figures/Fig5b_mechanism_panel.py` additionally needs the icon set under `data/fig5img/`,
which is not redistributed here.

---

## 8. Command lines behind the reported trees

```bash
# Genus-wide tree of 141 genomes (Fig. 2a). No -alrt, no fixed seed; the seed
# recorded in the log is 163984. Best-fit model Q.INSECT+R7.
iqtree -s phylo_aligned.fasta -m MFP -bb 1000 -nt AUTO

# Eight-strain supermatrix species tree (Fig. 2b). Best-fit model Q.PLANT+F+I+G4,
# rooted on S. livingstonensis afterwards.
iqtree -s supermatrix.faa -m MFP -bb 1000 -alrt 1000 -nt AUTO -ntmax 8 --seed 42 \
       --redo -pre Sv8_supermatrix

# Protein gene trees for the four island genes (Supplementary Fig. 3).
iqtree -s {OG}_prot.aln -m TEST -B 1000 -alrt 1000 --prefix {OG}_tree -T 4 --seed 42 --redo
#   GH3_e108 -> Q.PFAM+G4;  bglB -> WAG+G4;  MFS -> Q.PLANT+G4
#   GH3_e227 -> Q.PLANT+G4 for the seven-taxon tree; the model of the eight-taxon
#   re-run behind the reported omega was not recorded
```

Alignment settings differ between the two uses of MAFFT: `--auto` for the supermatrix
orthogroups, and L-INS-i (`--maxiterate 1000 --localpair`) for the single-gene island trees.
`trimal -automated1` was applied to the supermatrix only.

---

## 9. How to cite the tools

- IQ-TREE 3: Wong, T. K. F. *et al.* 2026. *Mol. Biol. Evol.* **43**, msag117
- ModelFinder: Kalyaanamoorthy, S. *et al.* 2017. *Nat. Methods* **14**, 587–589
- UFBoot2: Hoang, D. T. *et al.* 2018. *Mol. Biol. Evol.* **35**, 518–522
- SH-aLRT: Guindon, S. *et al.* 2010. *Syst. Biol.* **59**, 307–321
- MAFFT 7: Katoh, K. & Standley, D. M. 2013. *Mol. Biol. Evol.* **30**, 772–780
- trimAl: Capella-Gutiérrez, S. *et al.* 2009. *Bioinformatics* **25**, 1972–1973
- Flye: Kolmogorov, M. *et al.* 2019. *Nat. Biotechnol.* **37**, 540–546
- CheckM2: Chklovski, A. *et al.* 2023. *Nat. Methods* **20**, 1203–1212
- Prokka: Seemann, T. 2014. *Bioinformatics* **30**, 2068–2069
- OrthoFinder: Emms, D. M. & Kelly, S. 2019. *Genome Biol.* **20**, 238
- DIAMOND: Buchfink, B. *et al.* 2021. *Nat. Methods* **18**, 366–368
- eggNOG-mapper v2: Cantalapiedra, C. P. *et al.* 2021. *Mol. Biol. Evol.* **38**, 5825–5829
- dbCAN3: Zheng, J. *et al.* 2023. *Nucleic Acids Res.* **51**, W115–W121
- CAZy: Drula, E. *et al.* 2022. *Nucleic Acids Res.* **50**, D571–D577
- PAML: Yang, Z. 2007. *Mol. Biol. Evol.* **24**, 1586–1591
- PAL2NAL: Suyama, M. *et al.* 2006. *Nucleic Acids Res.* **34**, W609–W612
- FastANI: Jain, C. *et al.* 2018. *Nat. Commun.* **9**, 5114
- EzAAI: Kim, D. *et al.* 2021. *J. Microbiol.* **59**, 476–480
- MUMmer4: Marçais, G. *et al.* 2018. *PLoS Comput. Biol.* **14**, e1005944
- fastp: Chen, S. *et al.* 2018. *Bioinformatics* **34**, i884–i890
- BWA-MEM: Li, H. 2013. arXiv:1303.3997
- CoverM: Aroney, S. T. N. *et al.* 2025. *Bioinformatics* **41**, btaf147
- SignalP 6.0: Teufel, F. *et al.* 2022. *Nat. Biotechnol.* **40**, 1023–1025
- PSORTb 3.0: Yu, N. Y. *et al.* 2010. *Bioinformatics* **26**, 1608–1615
