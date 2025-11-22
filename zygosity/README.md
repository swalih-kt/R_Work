🧬 Figure 6 — Genotype Count Heatmaps (Zygosity per Variant / Gene)

This pipeline generates three high-resolution heatmaps representing genotype distribution of variants and genes based on heterozygous (het) and homozygous (hom) status.
The plots are combined into a single publication-quality figure.

📌 Overview of Panels
Panel	Description
A	Variant-level genotype counts across samples
B	Gene-level heterozygous/homozygous counts
C	Gene-level unique heterozygous/homozygous counts
📁 Required Input Files

Place all files in the same working directory before running the script.

Filename	Required Columns	Description
f6_var.txt	variant, zygosity, Individuals_count	Variant-level genotype counts
f6_gene_hethom.txt	gene, zygosity, Individuals_count	Het/Hom counts per gene
f6_gene_uhet_hom.txt	gene, zygosity, Individuals_count	Unique Het/Hom counts per gene

⚠ Column names must match exactly.

📂 Directory Structure
/home/treesa/R_files/input_files/cr_plot_md/updated_x/
  ├── f6_var.txt
  ├── f6_gene_hethom.txt
  ├── f6_gene_uhet_hom.txt
  └── (script will generate) Figure6_Heatmaps.png

▶️ R Script to Generate Figure 6
# Install if needed
# install.packages("ggpubr")

library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

setwd("/home/treesa/R_files/input_files/cr_plot_md/updated_x/")

# ---------- 1. Read input files ----------
a <- read.table("./f6_var.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- read.table("./f6_gene_hethom.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
c <- read.table("./f6_gene_uhet_hom.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---------- 2. Function to shorten long names ----------
shorten_name <- function(x, maxlen = 32) {
  sapply(x, function(name) {
    if (nchar(name) > maxlen) {
      trimmed <- substr(name, 1, maxlen)
      paste0(trimmed, "_", nchar(name) - maxlen, "nts")
    } else {
      name
    }
  })
}

a$variant <- shorten_name(a$variant)
b$gene    <- shorten_name(b$gene)
c$gene    <- shorten_name(c$gene)

# ---------- 3. Function to create a heatmap ----------
create_heatmap <- function(df, x_col, y_col, fill_col, title_text) {
  ggplot(df, aes_string(x = x_col, y = y_col, fill = fill_col)) +
    geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
    scale_fill_gradient(low = "#f3eded", high = "#4e9fe5", na.value = "#f3eded") +
    geom_text(data = dplyr::filter(df, !!as.name(fill_col) != 0),
              aes_string(label = fill_col), color = "black", size = 4) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 0.9, size = 14),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 12),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(0.5, "cm")
    ) +
    labs(title = title_text, x = x_col, y = y_col, fill = "Individuals count")
}

# ---------- 4. Create heatmap panels ----------
a1 <- create_heatmap(a, "zygosity", "variant", "Individuals_count",
                     "A) Variant-level Genotype Counts")
b1 <- create_heatmap(b, "zygosity", "gene", "Individuals_count",
                     "B) Gene-level Het/Hom Counts")
c1 <- create_heatmap(c, "zygosity", "gene", "Individuals_count",
                     "C) Gene-level Unique Counts")

# ---------- 5. Combine heatmaps ----------
combined <- ggarrange(
  a1,
  ggarrange(b1, c1, ncol = 2, nrow = 1, labels = c("B", "C"), font.label = list(size = 14)),
  ncol = 2, widths = c(2, 2), labels = c("A", ""), font.label = list(size = 14)
)

# ---------- 6. Save final figure ----------
png("./Figure6_Heatmaps.png", pointsize = 18, res = 300, width = 6000, height = 4000)
annotate_figure(
  combined,
  top = text_grob("Figure 6: Genotype Count Heatmaps", x = 0, hjust = 0, face = "bold", size = 18)
)
dev.off()

🎯 Output

The script generates:

Figure6_Heatmaps.png


A high-resolution (300 dpi) figure suitable for journal publication.
