# DESeq2 Differential Expression and GO/KEGG Enrichment Pipeline

This repository provides a complete **RNA-seq differential expression and functional enrichment analysis workflow**
using **DESeq2** and **clusterProfiler**.

The pipeline compares **Copper-treated Control (Cu_Cnt)** and **Copper-treated WD (Cu_WD)** samples and performs:

- Differential gene expression analysis
- PCA for sample clustering
- GO Biological Process enrichment
- KEGG pathway enrichment
- Publication-ready visualizations

---

## 📌 Analysis Overview

**Comparison:**  
- Cu_Cnt (n = 6) vs Cu_WD (n = 6)

**Organism:**  
- Human (*Homo sapiens*)

**Gene IDs:**  
- ENSEMBL (converted to ENTREZ for enrichment)

**Statistical thresholds:**  
- |log2FoldChange| ≥ 2  
- Adjusted p-value (BH) ≤ 0.05  

---

## 📦 Requirements

### R (version ≥ 4.1 recommended)

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
