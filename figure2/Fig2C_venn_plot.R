
library(VennDiagram)
#input file
set1<- read.table("ARG.txt",header = T,row.names = 1,sep = "\t")
set2<- read.table("TARGET.txt",header = T,row.names = 1,sep = "\t")
set3 <- read.table("GSE21257.txt",header = T,row.names = 1,sep = "\t")

set1 <- row.names(set1)
set2 <- row.names(set2)
set3 <- row.names(set3)

# Ensure that the dataset is in the form of character vectors.
set1 <- as.character(set1)
set2 <- as.character(set2)
set3 <- as.character(set3)


plot <- venn.diagram(
  x = list(set1,set2,set3),
  category.names = c("ARG", "TARGET","GSE21257"),
  filename =NULL,  
  output=FALSE,
  # Circle Properties:
  col = "black", 
  lty = 1, 
  lwd = 1, 
  fill = c("#C1C1C1","#7AAACB","#E89874"),
  alpha = 0.60, 
  label.col = "black",
  cex = .5, 
  fontfamily = "serif",
  fontface = "bold",
  
  # Collection Name Attribute:
  cat.col = c("#C1C1C1","#7AAACB","#E89874"),
  cat.cex = .6,
  cat.fontfamily = "serif"
)

pdf("venn_diagram.pdf", width=6, height=5)
grid.draw(plot)
dev.off()


