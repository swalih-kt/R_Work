
library(trackViewer)

# Read input file
data <- read.table("final_coverage/shimal_alg11/final_exon/ALG-PLOT - ALG12_Lollipop.tsv", 
                   header = TRUE, sep = "\t")

# ---- Add new variant between exon2 and exon3 ----
# Create a new row with the same columns as your data
new_variant <- data.frame(
  HGVS_VALIDATED   = "NM_024105.4:c.500A>G",          # random example
  ACMG             = "Likely pathogenic (II)",        # random ACMG
  freq             = 6,                               # your chosen frequency
  Exon             = "exon2.5",                       # between exon2 and exon3
  AA.change        = "Lys167Arg",                     # random amino acid change
  AAchange         = "K167R",                         # short form
  Mutation.Effect  = "Missense",                      # effect type
  RANK             = 2                                # your chosen rank
)

# Append it
data <- rbind(data, new_variant)


# ---- Continue as before ----
sel <- as.numeric(gsub("exon", "", gsub("exon2.5", "2.5", data$Exon)))  # allow fractional exon
data$Exon[which(data$freq == max(data$freq))]

names <- rep("", length(sel))
names[which(data$freq == max(data$freq))] <- data$Exon[which(data$freq == max(data$freq))]

Population <- GRanges("chr1", IRanges(sel, width = 1))
Population$score <- data$freq
Population$label.parameter.rot <- 0.09

group_colors <- c("#ef476f","#8E44AD")
group <- data$RANK
Population$color <- group_colors[group]
Population$cex <- 0.8

feature.gr <- GRanges("chr1", IRanges(unique(floor(sel)), width = rep(1.9, length(unique(floor(sel))))))

feature.gr$fill <- c("#FF5733", "#33FF57", "#3357FF",
                     "#F1C40F", "#8E44AD",
                     "#d3dc96","#F1B6DA",
                     "#895129","#073b4c", "#1ABC9C",
                     "#f4a261")

# ---- Custom x-axis ----
xaxis <- c(1, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

# Plot
lolliplot(Population, feature.gr, 
          legend = NULL, 
          ranges = GRanges("chr1", IRanges(1, 13)), 
          xaxis =  c(FALSE), 
          yaxis = c(FALSE))







# Add custom labels
grid.text("NM_024105.4", x = .26, y = .9, just = "bottom")


# Define the legend data
#legend_labels <- c( "Missense","Frameshift_deletion")
legend_labels <- c( "Frameshift_deletion","Missense")
                     #"Nonframeshift_insertion","Nonframeshift_substitution","Missense","Startloss","Stopgain")
legend_colors <-  c( "#ef476f","#8E44AD")
#"#3357FF","#8E44AD","#00DBFF","#f4a261")
grid.text("Type", x = .32, y = .8, just = "bottom", gp = gpar(cex = 1.2))

grid.text("ALG12", x = .6, y = .98, just = "top", gp = gpar(cex = 1.5, fontface = "bold"))

# Create a custom two-column legend at the top-right using grid
legend_x <- unit(0.33, "npc")  # Position on the right side (near 85% of the plot width)
legend_y <- unit(0.74, "npc")  # Position at the top (near 95% of the plot height)

# Draw the legend on the top-right with a two-column layout
legend_grob <- legendGrob(labels = legend_labels, pch = 20, gp = gpar(col = legend_colors, cex = 1.2), ncol = 1)
grid.draw(editGrob(legend_grob, vp = viewport(x = legend_x, y = legend_y)))




# Define the legend data
legend_labels <- c( "Exon3", "Exon4", "Exon5", "Exon6",  "Exon8", "Exon10", "Exon11")
#legend_labels <- c("Exon3", "Exon4")
                   #"Exon5", "Exon6", "Exon7", "Exon8", "Exon9", "Exon10", "Exon11"

legend_colors <- c( "#F1C40F","#FF5733","#8E44AD","#3357FF", "#33FF57","#d3dc96","#F1B6DA")
#legend_colors <- c(  "#f4a261", "#1ABC9C","#073b4c", "#895129","#FF5733", "#F1B6DA","#3357FF","#d3dc96"
                    # ,"#8E44AD","#33FF57", "#F1C40F")
#legend_colors <- c("#00DBFF" ,"#FF5733","#f4a261", "#1ABC9C","#ef476f","#073b4c", "#3357FF", "#895129", "#8E44AD"
                 #  ,"#33FF57","#F1C40F" )
#legend_colors <- c( "#1ABC9C","#073b4c","#8E44AD","#00DBFF","#f4a261","#FF00B1", "#3357FF","#ef476f","#A7ABDE", "#895129", "#F1C40F",
                   # "#FF5733","#33FF57" )


#legend_colors <- c("#33FF57", "#F1C40F","#073b4c","#FF5733","#f4a261", "#3357FF","#1ABC9C","#8E44AD","#895129" )



grid.text("Exon", x = .65, y = .8, just = "bottom", gp = gpar(cex = 1.2))


# Create a custom two-column legend at the top-right using grid
legend_x <- unit(0.63, "npc")  # Position on the right side (near 85% of the plot width)
legend_y <- unit(0.7, "npc")  # Position at the top (near 95% of the plot height)

# Draw the legend on the top-right with a two-column layout
legend_grob <- legendGrob(labels = legend_labels, pch = 15, gp = gpar(col = legend_colors, cex = 1.2), ncol = 2)
grid.draw(editGrob(legend_grob, vp = viewport(x = legend_x, y = legend_y)))






