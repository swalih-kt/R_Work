
## 🔷 Figure — HLCS Gene Lollipop Plot (Mutation Distribution Across Protein Domains)

This script generates a **publication-quality lollipop plot** showing mutation positions across the **HLCS gene**, annotated by mutation type and protein domains.

---

### 📌 Input File Format (TSV)

| Column Name | Description |
|-------------|-------------|
| `POS` | Amino acid / nucleotide position of variant |
| `AACHANGE` | Amino-acid change annotation |
| `frequency` | Frequency of mutation |
| `RANK` | Grouping index representing mutation class (e.g., Frameshift, Nonsynonymous, etc.) |

Example input file:  
`HLCS_lollipop_Final.tsv`

---

### 📊 Output

| Output | Description |
|--------|-------------|
| *(screen plot)* | HLCS mutation lollipop plot with mutation types and protein domains |
| *(optional PNG/PDF export can be added if required)* | — |

---

### 🧬 R Code

```r
# ========================================================
# 📌 HLCS Mutation Lollipop Plot — Protein Domain Annotation
# ========================================================

setwd("/home/treesa/R_files/input_files/HLCS")

library(trackViewer)
library(grid)  # for manual legend placement

# Load mutation data
data <- read.table("HLCS_lollipop_Final.tsv", header = TRUE, sep = "\t")

# Mutation positions
sele <- data$POS

# Assign label only to mutation with highest frequency
names <- rep("", length(sele))
names[which(data$frequency == max(data$frequency))] <- data$AACHANGE[which(data$frequency == max(data$frequency))]

# GRanges object for mutation events
y <- GRanges("chr1", IRanges(sele, width = 1))
y$score <- data$frequency
y$label.parameter.rot <- 45

# Mutation-type colors
group_colors <- c("#ff4000", "#990099", "#ff00d6", "#118ab2", "#66c2a5")
group <- data$RANK
y$color <- group_colors[group]
y$cex <- 0.6
y$border <- sample(c("gray80"))

# GRanges object for protein domains
feature <- GRanges(
  "chr1",
  IRanges(c(601, 816), width = c(198, 47))
)
feature$fill <- c("#66c2a5", "#fc8d62")

# Draw lollipop plot
lolliplot(
  y, feature,
  legend = NULL,
  ranges = GRanges("chr1", IRanges(1, 873)),
  yaxis = FALSE
)

# Axis positions
grid.text("1", x = 0.1, y = 0.073)
grid.text("|", x = 0.101, y = 0.105)
grid.text("873", x = 0.942, y = 0.073)
grid.text("|", x = 0.944, y = 0.105)

# Title & subtitles
grid.text("NM_001352514.2", x = .18, y = .85, gp = gpar(cex = 2))
grid.text("Mutation_Type", x = .198, y = .72, gp = gpar(cex = 2))
grid.text("HLCS", x = .5, y = .98, gp = gpar(cex = 1.5, fontface = "bold"))

# Legend — Mutation types
legend_labels <- c(
  "Frameshift_deletion", "Frameshift_insertion", "Frameshift_substitution",
  "Nonsynonymous", "Stopgain"
)
legend_colors <- c("#ff4000", "#990099", "#ff00d6", "#118ab2", "#66c2a5")

legend_x <- unit(0.3, "npc")
legend_y <- unit(0.62, "npc")
legend_grob <- legendGrob(labels = legend_labels, pch = 20,
                          gp = gpar(col = legend_colors, cex = 2),
                          ncol = 2)
grid.draw(editGrob(legend_grob, vp = viewport(x = legend_x, y = legend_y)))

# Legend — Protein domains
legend_labels <- c("BPL_LPL_Catalytic_Domain", "BPL_C_Terminal_Domain")
legend_colors <- c("#66c2a5", "#fc8d62")

grid.text("Domain", x = .72, y = .72, gp = gpar(cex = 2))

legend_x <- unit(0.8, "npc")
legend_y <- unit(0.65, "npc")
legend_grob <- legendGrob(labels = legend_labels, pch = 15,
                          gp = gpar(col = legend_colors, cex = 2),
                          ncol = 1)
grid.draw(editGrob(legend_grob, vp = viewport(x = legend_x, y = legend_y)))
