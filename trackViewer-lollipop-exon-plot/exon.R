library(trackViewer)

# Read input file
data <- read.table("final_lollipop/PFIC_final/EXONE/Lollipop_PFIC_ABCB11_input_exon.tsv", 
                   header = TRUE, sep = "\t")



# ---- Continue as before ----
sel <- as.numeric(gsub("exon", "", data$Exon))
sel
data$Exon[which(data$freq == max(data$freq))]

names <- rep("", length(sel))
names[which(data$freq == max(data$freq))] <- data$Exon[which(data$freq == max(data$freq))]

y <- GRanges("chr1", IRanges(sel, width = 1))
y$score <- data$freq
y$label.parameter.rot <- 0.09

group_colors <- c("#ef476f","#8E44AD", "#808000","#FFB6C1", "#00BFFF", "#073b4c", "#1ABC9C","#33FF57","#3357FF","#FF8C00")
group <- data$RANK
y$color <- group_colors[group]
y$cex <- 0.4
y$border <- sample(c("gray80"))
feature.gr <- GRanges("chr1", IRanges(unique(floor(sel)), width = rep(1.9, length(unique(floor(sel))))))

feature.gr$fill <- c(
  "#FF0000", "#00FF00", "#0000FF", "#FFFF00", 
  "#FF00FF", "#00FFFF", "#FFB6C1", "#E6E6FA", "#E0FFFF",
  "#F5FFFA", "#F08080", "#FFFACD", "#808000", "#008080",
  "#800000", "#000080", "#A0522D", "#708090", "#800080",
  "#4B0082", "#DC143C", "#FF8C00", "#228B22", "#00BFFF", 
  "#808080", "#C0C0C0", "#36454F"
)

# ---- Custom x-axis ----
#xaxis <- c(1, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
xaxis <- 1:27
# Plot
lolliplot(y, feature.gr, 
          legend = NULL, 
          ranges = GRanges("chr1", IRanges(1, 27)), 
          xaxis =  xaxis, 
          yaxis = c(FALSE))



# Add custom labels
grid.text("NM_003742.4", 
          x = 0.1, y = 0.9, 
          just = "bottom",
          gp = gpar(fontsize = 16, fontface = "bold"))


# Define the legend data
legend_labels <- c("Frameshift_deletion", "Frameshift_insertion", "Nonframeshift_deletion",
                   "Nonsynonymous_SNV","splicing","startloss","stopgain","Synonymous_SNV")


legend_labels <- c("Frameshift_deletion", "Frameshift_insertion", "Frameshift_substitution","Nonframeshift_deletion","Nonframeshift_substitution",
                   "Nonsynonymous_SNV","splicing","startloss","stopgain")

legend_colors <- c("#ef476f","#8E44AD", "#808000","#FFB6C1", "#00BFFF", "#073b4c", "#1ABC9C","#33FF57","#3357FF","#FF8C00")
#legend_colors <- c("#ef476f", "#ffd166","#06d6a0","#118ab2","#808080","#FF5733", "#33FF57", "#3357FF")
grid.text("Mutation_type", x = .25, y = .82, just = "bottom", gp = gpar(cex = 1.7))

grid.text("ABCB4", x = .5, y = .98, just = "top", gp = gpar(cex = 1.5, fontface = "bold"))

# Create a custom two-column legend at the top-right using grid
legend_x <- unit(0.3, "npc")  # Position on the right side (near 85% of the plot width)
legend_y <- unit(0.68, "npc")  # Position at the top (near 95% of the plot height)

# Draw the legend on the top-right with a two-column layout
legend_grob <- legendGrob(labels = legend_labels, pch = 20, gp = gpar(col = legend_colors, cex = 1.7), ncol = 2)
grid.draw(editGrob(legend_grob, vp = viewport(x = legend_x, y = legend_y)))
