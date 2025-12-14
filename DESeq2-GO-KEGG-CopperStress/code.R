#setwd("/home/swalih/Desktop/SS_CU/input_count_files")

# ----------------------------
# Load DESeq2
# ----------------------------
library(DESeq2)

# Read counts file
counts <- read.table("copper_work/ss_cu/cnt+cu&wd+cu/cont+cu&wd+cu.tsv", header=TRUE, sep="\t", row.names=1)

# Inspect
head(counts)

# ----------------------------
# Build sample information
# ----------------------------
# Columns: Cu_Cnt1a ... Cu_Cnt2c (6 samples) → group Cu_Cnt
#          Cu_WD2a ... Cu_WD1c (6 samples)   → group Cu_WD

sample_names <- colnames(counts)
condition <- c(rep("Cu_Cnt", 6), rep("Cu_WD", 6))   # 6 vs 6

colData <- data.frame(
  row.names = sample_names,
  condition = factor(condition)
)

# ----------------------------
# Create DESeq dataset
# ----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = colData,
  design    = ~ condition
)

# Pre-filter low count genes
dds <- dds[rowSums(counts(dds)) > 1, ]
# Run DESeq2
dds <- DESeq(dds)

# ----------------------------
# Get results
# ----------------------------
res <- results(dds, contrast = c("condition", "Cu_WD", "Cu_Cnt"))

# Order by adjusted p-value
res <- res[order(res$padj), ]

# Filter DEGs (|log2FC| >= 1, padj <= 0.05)
degs <- subset(res, abs(log2FoldChange) >= 2 & padj <= 0.05)

# Save results
write.csv(as.data.frame(res),  file="copper_work/ss_cu/DESeq2_all_cnt+cu&wd+cu.csv")
write.csv(as.data.frame(degs), file="DESeq2_filter__cnt+cu&wd+cu.csv")

# ----------------------------
# PCA Plot
# ----------------------------
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")



# ----------------------------
# Load libraries
# ----------------------------
library(clusterProfiler)
library(org.Hs.eg.db)   # Human annotation
library(enrichplot)
library(ggplot2)

# ----------------------------
# ----------------------------
# Split DEGs
# ----------------------------
up_genes   <- subset(degs, log2FoldChange > 2)
down_genes <- subset(degs, log2FoldChange < 2)

up_symbols   <- rownames(up_genes)
down_symbols <- rownames(down_genes)

# Background (all tested genes)
bg_symbols <- rownames(res)

# ----------------------------
# Convert ENSEMBL → ENTREZ
# ----------------------------
up_df <- bitr(
  up_symbols,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

down_df <- bitr(
  down_symbols,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

bg_df <- bitr(
  bg_symbols,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# ----------------------------
# GO Enrichment
# ----------------------------
ego_up <- enrichGO(
  gene     = up_df$ENTREZID,
  universe = bg_df$ENTREZID,
  OrgDb    = org.Hs.eg.db,
  keyType  = "ENTREZID",
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

ego_down <- enrichGO(
  gene     = down_df$ENTREZID,
  universe = bg_df$ENTREZID,
  OrgDb    = org.Hs.eg.db,
  keyType  = "ENTREZID",
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# ----------------------------
# KEGG Enrichment
# ----------------------------
ekegg_up <- enrichKEGG(
  gene     = up_df$ENTREZID,
  universe = bg_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

ekegg_down <- enrichKEGG(
  gene     = down_df$ENTREZID,
  universe = bg_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# ----------------------------
# Helper: capitalize terms
# ----------------------------
capitalize_terms <- function(x) {
  sapply(x, function(s) {
    s <- tolower(s)
    s <- gsub("(^|[[:space:]])([a-z])", "\\1\\U\\2", s, perl=TRUE)
    return(s)
  })
}

# Capitalize descriptions
ego_up@result$Description   <- capitalize_terms(ego_up@result$Description)
ego_down@result$Description <- capitalize_terms(ego_down@result$Description)
ekegg_up@result$Description <- capitalize_terms(ekegg_up@result$Description)
ekegg_down@result$Description <- capitalize_terms(ekegg_down@result$Description)

# ----------------------------
# Barplots: Upregulated
# ----------------------------
barplot(ego_up, showCategory=20, title="GO BP Upregulated", x="GeneRatio") +
  theme(axis.text.y=element_text(size=12, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13, face="bold"))

barplot(ekegg_up, showCategory=20, title="KEGG Upregulated", x="GeneRatio") +
  theme(axis.text.y=element_text(size=12, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13, face="bold"))

# ----------------------------
# Barplots: Downregulated
# ----------------------------
barplot(ego_down, showCategory=20, title="GO BP Downregulated", x="GeneRatio") +
  theme(axis.text.y=element_text(size=12, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13, face="bold"))

barplot(ekegg_down, showCategory=20, title="KEGG Downregulated", x="GeneRatio") +
  theme(axis.text.y=element_text(size=12, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        plot.title=element_text(size=16, face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13, face="bold"))

options(enrichplot.colours = c("#59a89c","#f0c571"))  #up

options(enrichplot.colours = c("#f0c571","#fbe9c4")) #down
options(enrichplot.colours = c("#75E0B0","#006072")) #kegg

# =======================================================
# Save enrichment results as TSV
# =======================================================
write.table(as.data.frame(ego_up), file="copper_work/ss_cu/cnt+cu&wd+cu/GO_UP_ORA.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(ego_down), file="copper_work/ss_cu/cnt+cu&wd+cu/GO_DOWN_ORA.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(ekegg_up), file="copper_work/ss_cu/cnt+cu&wd+cu/KEGG_UP_ORA.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(as.data.frame(ekegg_down), file="copper_work/ss_cu/cnt+cu&wd+cu/KEGG_DOWN_ORA.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
