# DESeq2 Differential Expression and GO/KEGG Enrichment Analysis

This repository contains an **RNA-seq differential expression and functional enrichment pipeline**
using **DESeq2** and **clusterProfiler**.

The analysis compares **Copper-treated control (Cu_Cnt)** vs **Copper-treated WD (Cu_WD)** samples
and performs:
- Differential gene expression analysis
- PCA visualization
- GO Biological Process enrichment
- KEGG pathway enrichment
- Publication-ready barplots

---

## 📦 Requirements

### R (≥ 4.1 recommended)

Install required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))

install.packages("ggplot2")
