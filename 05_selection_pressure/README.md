# Selection-pressure analysis — how these scripts relate to the published numbers

The four ω values reported in the manuscript and in Supplementary Fig. 4 were **not**
produced by a single pass through `01`–`07`. The gene set changed during the analysis
and the scripts were never retro-fitted, so they are kept here as the record of what
was actually run. Read this before trying to re-execute them.

## Execution history

| Gene | Orthogroup | ω (M0) | Produced by |
|---|---|---|---|
| bglB (GH1) | OG0001451 | 0.0425 | `01`–`03`, `04`, `05`–`07` |
| GH3_e108 | OG0002997 | 0.0634 | `01`–`03`, `04`, `05`–`07` |
| MFS transporter | OG0000506 | 0.0457 | `04b` (MAFFT) + `04` (tree) + `06b` (PAL2NAL → labelling → codeml) |
| GH3_e227 | OG0001452 | 0.0759 | eight-taxon re-run, see below |

## Three things the scripts do not show

**1. The `OGS` lists in `01`-`03` and `05`-`07` are the pre-May-2026 gene set.**
They do not contain the MFS transporter. They originally also contained
`OG0001992_GH3_NagZ`, which lies outside the island and was dropped from the manuscript;
it has been removed from the lists and each script header records that it was part of the
original run. Only `04_iqtree_gene_trees.sh` was updated (2026-05-28) to the final
four-gene set. The NagZ results on disk are superseded and are not reported anywhere in
the paper.

**2. The results path in `01`-`03` and `05`-`07` has been corrected.**
These scripts originally pointed to `paml_v2/results`, a directory name that no longer
exists; the real results live in `04_selection_pressure/paml/results`, which is now used
by all scripts (as it already was by `04`, `04b` and `06b`).

**3. GH3_e227 was re-run on eight taxa.**
Prokka called a truncated 838-aa ORF for the E1 copy (`PB002_01811`) because of a
homopolymer that shifts the reading frame, and the M7 copy was split across two ORFs
(`KDH10_01696`, `KDH10_01695`). Both were reassembled by removing one base from the homopolymer, giving 873 aa for E1
and 878 aa for M7; `00c_reassemble_GH3_e227_full_length.py` reproduces both and asserts the
result against the stored sequences. The six congeners were taken unchanged.
The resulting eight-taxon codon alignment was re-run through codeml on 2026-07-16 and is
the source of ω = 0.0759. That run also used `cleandata = 1`, whereas every other gene
here used `cleandata = 0`.

Two earlier GH3_e227 runs exist on disk and are **not** the published values: a
seven-taxon run on the Prokka ORF (ω = 0.0717) and a seven-taxon run on a trimmed
alignment (ω = 0.0723).

## Taxon counts

`bglB` and `GH3_e108` have seven taxa, not eight. For `GH3_e108` OrthoFinder recovered no
ortholog in *S. psychromarinicola*. For `bglB` the *S.* sp014164505 copy was excluded by hand:
it sits 99 bp from the end of its contig and is 89 aa short, which `01_extract_island_gene_sequences.py`
records as an assembly artefact rather than a real truncation. The MFS transporter and the
eight-taxon GH3_e227 alignment have all eight.

## Likelihood-ratio tests

`07_likelihood_ratio_tests.py` refers the branch-site statistic to χ²(1). The manuscript
uses the standard 50:50 χ²(0):χ²(1) mixture, which is the less conservative of the two;
every branch-site statistic is zero here, because the maximum-likelihood estimate of omega2
sits on the boundary at 1, so both referents give P = 1 and the two agree. The P values printed beside the bars in Supplementary
Fig. 4 are the **branch** LRT (two-ratio vs M0), for which χ²(1) is the correct null.

## One ordering constraint

`07_likelihood_ratio_tests.py` writes `paml_LRT_summary.tsv` from scratch, while
`06b_run_codeml_MFS_transporter.py` adds the MFS row to the same file. Run `07` first and
`06b` after it, which is the reverse of the numbering; running them the other way round
drops the MFS row.
