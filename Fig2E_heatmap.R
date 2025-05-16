
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# Load the expression data and annotation data
DEG <- read.table("diffExp.txt",header=T, sep="\t", check.names=F,row.names = 1)

annotation_col <- read.table("annotation_col.txt",header=T, sep="\t", check.names=F,row.names = 1)

annotation_row <- read.table("annotation_row.txt",header=T, sep="\t", check.names=F,row.names = 1)
#Convert to factor type
annotation_col[] <- lapply(annotation_col, as.factor)

annotation_row[] <- lapply(annotation_row, as.factor)
#Set the name of the row for the comment
rownames(annotation_row) <- rownames(DEG)

rownames(annotation_col) <- colnames(DEG)

#Set annotation color
annotation_colors <- list(
  Regulation=c("Up"="#FC9F5B","Down"="#25998F"),
  Gender = c("Male"="#0AA0BF","Female"="#F36E98"),
  Fustat = c("Alive"="#78B177","Dead"="#F6114A"),
  Group = c("Metastases"="#CF0BF1","No_metastases"="#C3C4E9")
)

# Draw a heat map
p <- pheatmap(DEG,
              scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(50),
         show_rownames = TRUE,
         show_colnames =FALSE,
         treeheight_row = 0,
         treeheight_col = 0 )
p

dev.off()
