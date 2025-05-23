  ####Immune cell infiltration heatmap####
library(pheatmap)
library(tidyverse)
results <- read.table("risk-violin.txt",header = T,sep = "\t",quote = "",row.names = 1)
#Start drawing after reading the grouped file
group_info<-read.table('riskinfo.txt',header=T,row.names = 1)
group_info<-group_info[order(match(group_info$sample,rownames(results))),]
results<- results[rownames(group_info), ]


annotation_row<-group_info[,"risk",drop=FALSE]

rownames(results)==rownames(annotation_row)

# Assign a color to each group
annotation_colors <- list(
  risk = c("high" = "#BB2222", "low" = "skyblue")
)

pdf("immune_cell_abundance_heatmap2.pdf",width=5,height=22)
pheatmap(results[,1:22],scale="row",
         clustering_distance_cols="euclidean",
         cluster_rows=F,
         cluster_cols=T,
         color=colorRampPalette(c("lightblue","white","#BB2222"))(100),
         show_rownames=F,
         show_colnames=TRUE,
         main="Immune Cell Abundance Heatmap",
         annotation_row=annotation_row,
         annotation_colors = annotation_colors,
         annotation_names_row=F,
         fontsize_col=10,
         cellwidth=10,
         cellheight=10
)
dev.off()

