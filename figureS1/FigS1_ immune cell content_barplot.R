
input <- "CIBERSORT.filter.txt"
group_file <- "risk.txt"
outpdf <- "barplot2.pdf"

# Read the data and group information
data <- read.table(input, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
data <- t(data)

group <- read.table(group_file, header=TRUE, sep="\t", row.names=1)
  group <- group[colnames(data), , drop=FALSE]  # Ensure consistency in the sequence.

  # Draw a bar chart using colors
  colors <- c(
    "#1b9e77",  # deep teal
    "#d95f02",  # warm orange
    "#7570b3",  # muted indigo
    "#e7298a",  # deep pink
    "#66a61e",  # olive green
    "#e6ab02",  # gold ochre
    "#a6761d",  # earth brown
    "#666666",  # medium gray
    "#66c2a5",  # light teal
    "#fc8d62",  # coral
    "#8da0cb",  # steel blue
    "#e78ac3",  # light pink
    "#a6d854",  # lime green
    "#ffd92f",  # sunflower
    "#e5c494",  # sand
    "#b3b3b3",  # light gray
    "#1f78b4",  # strong blue
    "#33a02c",  # vibrant green
    "#fb9a99",  # soft red
    "#cab2d6",  # lavender
    "#fdbf6f",  # peach
    "#b15928"   # brown ochre
)
  

pdf(outpdf, height=10, width=25)
par(las=1, mar=c(8,4,4,15))

a1 <- barplot(data, col=colors, yaxt="n", ylab="Relative Percent", xaxt="n",width = 0.0001,space = 0.2)
a2 <- axis(2, tick=FALSE, labels=FALSE)
axis(2, a2, paste0(a2*100, "%"))
axis(1, a1, labels=FALSE)

# Display the "high" and "low" labels (instead of the sample names)
group_labels <- as.character(group$risk)
group_colors <- ifelse(group_labels == "high", "#ff0000", "#00ffff")

bar_width <- mean(diff(a1))  # Estimate the column width
y_min <- -0.05  
y_max <- -0.01  

for (i in seq_along(a1)) {
  rect(xleft = a1[i] - bar_width / 2,
       xright = a1[i] + bar_width / 2,
       ybottom = y_min,
       ytop = y_max,
       col = group_colors[i],
       border = NA, xpd = TRUE)
}

par(srt=60, xpd=TRUE)
text(a1, -0.02, group_labels, adj=1, cex=0.6, col=group_colors)
par(srt=0)

# Add legend
legend(par('usr')[2]*0.98, par('usr')[4], legend=rownames(data), col=col, pch=15, bty="n", cex=1.3)

dev.off()
