# R_Work — Bioinformatics Visualization & Analysis Scripts

A collection of R scripts for bioinformatics data analysis and publication-quality figure generation, developed for genomic and transcriptomic research. Each folder is a self-contained analysis module with its own input files, R script, and README.

> **Author**: Mohammed Swalih K T
> **Language**: R (98.3%) | CSS (1.7%)

---

## Repository Structure

```
R_Work/
├── DESeq2-GO-KEGG-CopperStress/       # Differential expression + functional enrichment
├── Pie_chart/                          # ACMG/variant classification pie charts
├── bar_plot/                           # Stacked horizontal bar plot (ACMG per gene)
├── lollpop_skt/                        # Lollipop plot (swalih-kt variant)
├── maf/                                # Population MAF comparison (gnomAD vs IndiGen)
├── motif_breakR/                       # TF binding motif disruption analysis
├── trackViewer-lollipop-exon-plot/     # Protein domain lollipop plot (HLCS)
└── zygosity/                           # Genotype count heatmaps
```

---

## Modules

### 🧬 DESeq2-GO-KEGG-CopperStress
Differential expression analysis of copper stress RNA-seq data followed by functional enrichment using Gene Ontology (GO) and KEGG pathway analysis.

| | |
|---|---|
| **Tool** | DESeq2, clusterProfiler |
| **Input** | Raw count matrix, sample metadata |
| **Output** | DE gene list, GO/KEGG enrichment plots |

---

### 🥧 Pie_chart
Generates pie charts summarizing variant classifications by ACMG category or database status across a gene/cohort.

| | |
|---|---|
| **Tool** | ggplot2 |
| **Input** | Variant table with ACMG classification |
| **Output** | Pie chart PNG |

---

### 📊 bar_plot
Produces a publication-quality stacked horizontal bar plot showing the number of variants per gene grouped by ACMG pathogenicity classification. Supports novel vs. known variant stratification.

| | |
|---|---|
| **Tool** | ggplot2, dplyr |
| **Input** | TSV with Gene, VCF, ACMG, Database_Status columns |
| **Output** | `Figure3.png` — stacked bar chart |

---

### 🍭 lollpop_skt
Lollipop plot for visualizing variant positions and frequencies along a gene or protein, customized for cohort-level variant data.

| | |
|---|---|
| **Tool** | trackViewer / ggplot2 |
| **Input** | Variant position and frequency table |
| **Output** | Lollipop plot PNG |

---

### 🌍 maf
Compares minor allele frequencies of variants of interest between the **IndiGen cohort** and global/sub-continental **gnomAD** populations. Performs Fisher's exact test with Bonferroni correction and generates a bubble plot annotated with ACMG classifications.

| | |
|---|---|
| **Tool** | ggplot2, dplyr, patchwork |
| **Input** | `all_variant.txt` (gnomAD), `var_indigen.txt` (IndiGen + ACMG) |
| **Output** | `Figure7.png` — population enrichment bubble plot |
| **Statistics** | Fisher's exact test, Bonferroni correction |

---

### 🔬 motif_breakR
Predicts the impact of regulatory SNPs on transcription factor (TF) binding motifs using **motifbreakR**. Scores variants against HOCOMOCO v11 and JASPAR 2024 databases, computes empirical p-values, and produces filtered hit tables and per-SNP/per-TF summaries.

| | |
|---|---|
| **Tool** | motifbreakR, MotifDb, BiocParallel |
| **Input** | `clean_fixed.bed` — BED file of regulatory variants |
| **Output** | Motif disruption TSVs, strong hit tables, SNP/TF summaries |
| **Databases** | HOCOMOCOv11-core-A/B/C, JASPAR 2024 |

---

### 🧩 trackViewer-lollipop-exon-plot
Generates a protein-level lollipop plot for **HLCS** gene mutations using `trackViewer`. Variants are color-coded by mutation type (RANK), sized by curation frequency, and overlaid on annotated protein domains (BPL_LPL Catalytic Domain, BPL C-Terminal Domain).

| | |
|---|---|
| **Tool** | trackViewer, grid |
| **Input** | `HLCS_lollipop_Final.tsv` — curated HLCS variants with ACMG, amino acid position, mutation type |
| **Output** | HLCS lollipop plot PNG |
| **Protein** | NM_001352514.2 (873 aa) |

**Mutation type color legend:**

| Mutation Type | Color |
|---|---|
| Frameshift deletion | ![#ff4000](https://placehold.co/15x15/ff4000/ff4000.png) `#ff4000` |
| Frameshift insertion | ![#990099](https://placehold.co/15x15/990099/990099.png) `#990099` |
| Frameshift substitution | ![#ff00d6](https://placehold.co/15x15/ff00d6/ff00d6.png) `#ff00d6` |
| Nonsynonymous | ![#118ab2](https://placehold.co/15x15/118ab2/118ab2.png) `#118ab2` |
| Stopgain | ![#66c2a5](https://placehold.co/15x15/66c2a5/66c2a5.png) `#66c2a5` |

---

### 🔥 zygosity
Generates a three-panel genotype count heatmap (Figure 6) showing the distribution of heterozygous and homozygous variants across individuals at both variant level and gene level. Color intensity represents the number of individuals carrying each genotype.

| | |
|---|---|
| **Tool** | ggplot2, ggpubr, dplyr |
| **Input** | `f6_var.txt`, `f6_gene_hethom.txt`, `f6_gene_uhet_hom.txt` |
| **Output** | `papper.png` — three-panel heatmap (6000 × 4000 px, 300 DPI) |

---

## Requirements

All scripts require **R v4.0 or later**. Package requirements vary by module — see each folder's README for details. Commonly used packages:

- `ggplot2`, `dplyr`, `tidyr`, `ggpubr`, `patchwork`
- `trackViewer`, `motifbreakR`, `MotifDb`
- `DESeq2`, `clusterProfiler`
- `BSgenome.Hsapiens.UCSC.hg38`, `SNPlocs.Hsapiens.dbSNP155.GRCh38`

---

## Notes

- Each folder contains its own `README.md` with detailed step-by-step documentation, input/output descriptions, and usage notes.
- All figures are exported as high-resolution PNG files (300 DPI) suitable for publication.
- Working directory paths in scripts are set to local machine paths — update `setwd()` before running.

---

## References

- [trackViewer — Bioconductor](https://bioconductor.org/packages/release/bioc/html/trackViewer.html)
- [motifbreakR — Bioconductor](https://bioconductor.org/packages/release/bioc/html/motifbreakR.html)
- [DESeq2 — Bioconductor](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [gnomAD Browser](https://gnomad.broadinstitute.org/)
- [JASPAR 2024](https://jaspar.elixir.no/)
- [HOCOMOCO v11](https://hocomoco11.autosome.org/)
