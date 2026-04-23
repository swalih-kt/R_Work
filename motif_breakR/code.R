# ==========================================================
# motifbreakR Analysis Pipeline
# ==========================================================
#
# Purpose:
#   Predict the impact of regulatory variants (SNPs) on
#   transcription factor binding motifs using motifbreakR.
#   Identifies SNPs that disrupt or create TF binding sites
#   using HOCOMOCO and JASPAR 2024 motif databases.
#
# Input:
#   clean_fixed.bed   - BED file of variants to analyze
#
# Output:
#   snps_output.tsv                    - Loaded SNP details
#   motifbreakR_results_no_pval.tsv    - All motif disruption results (no p-value)
#   motifbreakR_results_with_pval.tsv  - Results with p-values
#   motifbreakR_final_hits_strong.tsv  - Strong effect hits only
#   SNP_summary.tsv                    - Per-SNP TF disruption summary
#   TF_summary_new.tsv                 - Per-TF affected SNP summary
#   pvalue_strong_hits.tsv             - Hits filtered by p-value strength
#
# Software:
#   R with motifbreakR, MotifDb, BSgenome, BiocParallel
# ==========================================================


# ----------------------------------------------------------
# Set working directory
# ----------------------------------------------------------

setwd("/home/srishti-sharma/Desktop/Srishti/motif/final_analysis/")


# ----------------------------------------------------------
# Load libraries
# ----------------------------------------------------------

library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(MotifDb)
library(BiocParallel)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(dplyr)


# ==========================================================
# Step 1: Load SNPs from BED file
# ==========================================================
# Reads variants from a BED file and maps them to dbSNP155
# rsIDs where possible using check.unnamed.for.rsid = TRUE.

snps <- snps.from.file(
  file    = "motif_2026/final_variants/input_variants.bed",
  format  = "bed",
  dbSNP   = SNPlocs.Hsapiens.dbSNP155.GRCh38,
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  check.unnamed.for.rsid = TRUE
)

# Save SNP details
snps_df <- as.data.frame(snps, row.names = NULL)

write.table(snps_df,
            file      = "motif_2026/final_variants/snps_output.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ==========================================================
# Step 2: Run motifbreakR (without p-values)
# ==========================================================
# Scores each SNP against TF binding motifs from HOCOMOCO
# (v11 core A/B/C) and JASPAR 2024.
# method = "ic" : information content scoring
# threshold = 1e-4 : minimum motif score threshold
# filterp = TRUE : pre-filter weak hits before scoring

results <- motifbreakR(
  snpList  = snps,
  filterp  = TRUE,
  pwmList  = subset(MotifDb,
                    dataSource %in% c("HOCOMOCOv11-core-A",
                                      "HOCOMOCOv11-core-B",
                                      "HOCOMOCOv11-core-C",
                                      "jaspar2024")),
  threshold = 1e-4,
  method    = "ic",
  bkg       = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  BPPARAM   = BiocParallel::SerialParam()
)

# Save results without p-values
results_df_r <- as.data.frame(results, row.names = NULL)
results_df_r$motifPos <- sapply(results_df_r$motifPos, paste, collapse = ",")

write.table(results_df_r,
            file      = "motif_2026/final_variants/motifbreakR_results_no_pval.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ==========================================================
# Step 3: Calculate p-values
# ==========================================================
# Computes empirical p-values for each motif disruption score.
# Uses multicore parallel processing for speed.

ncores <- parallel::detectCores() - 1

results_with_pval <- calculatePvalue(
  results,
  BPPARAM = MulticoreParam(workers = ncores)
)

# Save results with p-values
results_df <- as.data.frame(results_with_pval, row.names = NULL)
results_df$motifPos <- sapply(results_df$motifPos, paste, collapse = ",")

write.table(results_df,
            file      = "motif_2026/final_variants/motifbreakR_results_with_pval.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ==========================================================
# Step 4: Filter strong effect hits
# ==========================================================
# Retains only SNPs with BOTH strong allele effect size
# AND strong p-value effect — highest confidence hits.

final_hits <- results_with_pval[
  results_with_pval$effect == "strong" &
    results_with_pval$pvalueEffect == "strong"
]

# Flatten list columns for export
final_hits_df <- as.data.frame(final_hits, row.names = NULL)
final_hits_df[] <- lapply(final_hits_df, function(x) {
  if (is.list(x)) sapply(x, paste, collapse = ",") else x
})

write.table(final_hits_df,
            file      = "motif_2026/final_variants/motifbreakR_final_hits_strong.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ==========================================================
# Step 5: Per-SNP summary
# ==========================================================
# For each SNP, counts how many TFs are disrupted and
# reports mean/max allele effect size.

snp_summary <- final_hits_df %>%
  group_by(SNP_id) %>%
  summarise(
    n_TFs       = n(),
    TFs         = paste(unique(geneSymbol), collapse = ","),
    mean_effect = mean(alleleEffectSize, na.rm = TRUE),
    max_effect  = max(abs(alleleEffectSize), na.rm = TRUE)
  ) %>%
  arrange(desc(n_TFs))

write.table(snp_summary,
            file      = "motif_2026/final_variants/SNP_summary.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ==========================================================
# Step 6: Per-TF summary
# ==========================================================
# For each transcription factor, counts the number of
# distinct SNPs that disrupt its binding motif.

tf_summary <- final_hits_df %>%
  group_by(geneSymbol) %>%
  summarise(
    n_SNPs = n_distinct(SNP_id),
    SNPs   = paste(unique(SNP_id), collapse = ",")
  ) %>%
  arrange(desc(n_SNPs))

write.table(tf_summary,
            file      = "motif_2026/final_variants/TF_summary_new.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ==========================================================
# Step 7: Filter by p-value strength only
# ==========================================================
# Alternative filter — retains hits with strong p-value
# regardless of allele effect size. Useful for comparison
# against the combined strong-effect filter in Step 4.

pval_strong_hits <- results_with_pval[
  results_with_pval$pvalueEffect == "strong"
]

pval_strong_df <- as.data.frame(pval_strong_hits, row.names = NULL)
pval_strong_df[] <- lapply(pval_strong_df, function(x) {
  if (is.list(x)) sapply(x, paste, collapse = ",") else x
})

write.table(pval_strong_df,
            file      = "motif_2026/pvalue_strong_hits.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)


# ==========================================================
# Step 8: Optional — Examine / Visualize a Single Variant
# ==========================================================
# Subsets results for a specific SNP by rsID, recalculates
# its p-value individually, and plots the motif disruption.
# plotMB shows the reference vs. alternate allele motif logo
# and highlights the effect of the SNP on TF binding.
# Useful for validating or presenting a specific variant
# of interest from the results.

rsV <- results[names(results) %in% "rs1479475149"]
rsV <- calculatePvalue(rsV)

plotMB(results = results,
       rsid    = "rs1479475149",
       effect  = "strong")


# ==========================================================
# End of motifbreakR Pipeline
# ==========================================================
