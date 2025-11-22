install.packages("ggpubr")

library(ggpubr)
setwd("/home/treesa/R_files/input_files/cr_plot_md/updated_x/")


# ============================
# Figure 6: Genotype count heatmaps (zygosity per variant/gene)
# ============================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)  # for unit()

# ---- 1. Read input files ----
a <- read.table("./f6_var.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- read.table("./f6_gene_hethom.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
c <- read.table("./f6_gene_uhet_hom.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- 2. Function to shorten long names ----
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

# Apply trimming
a$variant <- shorten_name(a$variant)
b$gene    <- shorten_name(b$gene)
c$gene    <- shorten_name(c$gene)

# ---- 3. Function to create heatmap ----
create_heatmap <- function(df, x_col, y_col, fill_col, title_text) {
  ggplot(df, aes_string(x = x_col, y = y_col, fill = fill_col)) +
    geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
    scale_fill_gradient(low = "#f3eded", high = "red", na.value = "#f3eded") +
    geom_text(data = dplyr::filter(df, !!as.name(fill_col) != 0),
              aes_string(label = fill_col), color = "black", size = 3) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 0.9, size = 6),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(face = "bold", size = 10)) +
    labs(title = title_text, x = x_col, y = y_col)
}

# ---- 4. Generate heatmaps ----
a1 <- create_heatmap(a, "zygosity", "variant", "Individuals_count", "A) Variant-level Genotype Counts")
b1 <- create_heatmap(b, "zygosity", "gene", "Individuals_count", "B) Gene-level Het/Hom Counts")
c1 <- create_heatmap(c, "zygosity", "gene", "Individuals_count", "C) Gene-level Unique Counts")

# ---- 5. Combine panels ----
combined <- ggarrange(
  a1,
  ggarrange(b1, c1, ncol = 2, nrow = 1, labels = c("B", "C"), font.label = list(size = 12)),
  ncol = 2, widths = c(2, 2), labels = c("A", ""), font.label = list(size = 12)
)

# ---- 6. Save Figure 6 ----
png("./paper_fig.png", pointsize = 14, res = 300, width = 6000, height = 4000)
annotate_figure(combined,
                top = text_grob("Figure 6: Genotype Count Heatmaps", x = 0, hjust = 0, face = "bold", size = 16))
dev.off()

















# ============================
# Figure 6: Genotype count heatmaps (zygosity per variant/gene)
# ============================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)  # for unit()

# ---- 1. Read input files ----
a <- read.table("./f6_var.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- read.table("./f6_gene_hethom.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
c <- read.table("./f6_gene_uhet_hom.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---- 2. Function to shorten long names ----
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

# Apply trimming
a$variant <- shorten_name(a$variant)
b$gene    <- shorten_name(b$gene)
c$gene    <- shorten_name(c$gene)

# ---- 3. Function to create heatmap ----
create_heatmap <- function(df, x_col, y_col, fill_col, title_text) {
  ggplot(df, aes_string(x = x_col, y = y_col, fill = fill_col)) +
    geom_tile(color = "grey", size = 0.1, width = 0.98, height = 0.99) +
    scale_fill_gradient(low = "#f3eded", high = "#4e9fe5", na.value = "#f3eded") +
    geom_text(data = dplyr::filter(df, !!as.name(fill_col) != 0),
              aes_string(label = fill_col), color = "black", size = 4) +   # text size for counts
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 0.9, size = 14),  # increased x-axis label size
      axis.text.y = element_text(size = 12),  # increased gene/variant label size
      plot.title = element_text(face = "bold", size = 12),
      legend.title = element_text(size = 14, face = "bold"),   # increased legend title
      legend.text  = element_text(size = 12),                  # increased legend text
      legend.key.height = unit(1.2, "cm"),                     # larger legend color bar
      legend.key.width  = unit(0.5, "cm")
    ) +
    labs(title = title_text, x = x_col, y = y_col, fill = "Individuals count")
}

# ---- 4. Generate heatmaps ----
a1 <- create_heatmap(a, "zygosity", "variant", "Individuals_count", "A) Variant-level Genotype Counts")
b1 <- create_heatmap(b, "zygosity", "gene", "Individuals_count", "B) Gene-level Het/Hom Counts")
c1 <- create_heatmap(c, "zygosity", "gene", "Individuals_count", "C) Gene-level Unique Counts")

# ---- 5. Combine panels ----
combined <- ggarrange(
  a1,
  ggarrange(b1, c1, ncol = 2, nrow = 1, labels = c("B", "C"), font.label = list(size = 14)),
  ncol = 2, widths = c(2, 2), labels = c("A", ""), font.label = list(size = 14)
)

# ---- 6. Save Figure 6 ----
png("./papper.png", pointsize = 18, res = 300, width = 6000, height = 4000)
annotate_figure(
  combined,
  top = text_grob("Figure 6: Genotype Count Heatmaps", x = 0, hjust = 0, face = "bold", size = 18)
)
dev.off()


