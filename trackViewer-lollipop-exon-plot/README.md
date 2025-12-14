# Lollipop Plot for Exonic Mutation Distribution using trackViewer

This repository contains an R script to generate a **lollipop plot** visualizing the distribution and frequency of genetic variants across **exons** of a gene using the **trackViewer** package.

The plot highlights:
- Mutation frequency per exon
- Mutation types using color coding
- Most frequent exon annotation
- Custom legend and axis formatting

---

## 📦 Requirements

Install the required R packages before running the script:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("trackViewer")
BiocManager::install("GenomicRanges")

install.packages("grid")
