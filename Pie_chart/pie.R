install.packages("ggpubr")
install.packages("cowplot")


🚀 Workflow
🔹 Section 1 — Variant Classification Pie Chart
# =============================================
# 1️⃣ Publication Pie Chart — Variant Class Types
# =============================================

counts <- c(339, 70, 19, 17, 44)
labels <- c("Missense", "Frameshift", "Splice Acceptor", "Splice Donor", "Stop Gained")

# Convert counts → percentage labels
pct <- round(counts / sum(counts) * 100, 1)
labels <- paste0(pct, "%")

colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

pie(
  counts,
  labels = labels,
  radius = 1,
  main = paste0("Variant Classification (Total = ", sum(counts), ")"),
  col = colors,
  clockwise = TRUE,
  cex = 1.8,
  font = 2
)

legend(
  x = 1.2, y = 0.5,
  legend = paste0(c("Missense", "Frameshift", "Splice Acceptor", "Splice Donor", "Stop Gained"), " (", counts, ")"),
  fill = colors,
  cex = 1,
  bty = "n",
  y.intersp = 0.8,
  text.font = 2
)

🔹 Section 2 — REVEL Score Pie Chart
# =============================================
# 2️⃣ Publication Pie Chart — REVEL Category Distribution
# =============================================

library(ggpubr)
library(cowplot)

counts <- c(310, 13, 11, 5)
labels <- expression(
  bold("310 (91.4%)"),
  bold("13 (3.8%)"),
  bold("11 (3.2%)"),
  bold("5 (1.5%)")
)

colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

png("./Figure2.png", pointsize = 18, res = 300, width = 6000, height = 4000)
par(mar = c(4, 4, 4, 0))

pie(
  counts,
  labels = labels,
  col = colors,
  cex = 1.5,
  cex.main = 3,
  cex.lab = 2,
  border = "white"
)

legend(
  x = 0.9, y = 1,
  legend = c("VUS", "LP (Iv or V)", "P/LP (II or III)", "LB"),
  fill = colors,
  cex = 2,
  bty = "n"
)

dev.off()

🔹 Section 3 — LoF Variant Classification Pie Chart
# =============================================
# 3️⃣ Publication Pie Chart — Loss-of-Function Variant Classes
# =============================================

library(ggpubr)
library(cowplot)

counts <- c(112, 36, 2)
labels <- expression(
  bold("112 (74.7%)"),
  bold("36 (24%)"),
  bold("2 (1.3%)")
)

colors <- c("#1f77b4", "#2ca02c", "#d62728")

png("./Figure3_pie.png", pointsize = 18, res = 300, width = 6000, height = 4000)
par(mar = c(4, 4, 4, 0))

pie(
  counts,
  labels = labels,
  col = colors,
  cex = 1.5,
  cex.main = 3,
  cex.lab = 2,
  border = "white"
)

legend(
  x = 0.9, y = 1,
  legend = c("VUS", "P", "LB"),
  fill = colors,
  cex = 2,
  bty = "n"
)

dev.off()

📌 Output Files Generated
Figure File	Description
(screen display only)	Variant consequence distribution
Figure2.png	REVEL score category distribution
Figure3_pie.png	LoF variant class distribution
