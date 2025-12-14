# =============================================
# 📌 Stacked Horizontal Bar Plot for ACMG by Gene
# =============================================

setwd("/home/treesa/R_files/input_files/cr_plot_md/novel")

library(ggplot2)
library(dplyr)
library(readr)

# Read the TSV file
data <- read_tsv("Swalih_MD_plot_input - novel_bar.tsv")

# Rename columns for simplicity
colnames(data) <- c("Gene", "VCF", "ACMG", "Database_Status")

# Summarize counts per Gene × ACMG classification
summary_data <- data %>%
  group_by(Gene, ACMG) %>%
  summarise(count = n(), .groups = "drop")

# Save the figure as high-resolution PNG
png("./Figure3.png", pointsize = 18, res = 300, width = 2500, height = 3000)

# Stacked horizontal bar plot
ggplot(summary_data, aes(x = count, y = reorder(Gene, count), fill = ACMG)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  labs(
    title = "Variant Counts per Gene by ACMG Classification",
    x = "Count",
    y = "Gene",
    fill = "ACMG Classification"
  ) +
  theme_minimal(base_size = 13) +
  scale_x_continuous(breaks = seq(2, 28, by = 2), limits = c(0, 28)) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right"
  )

dev.off()
