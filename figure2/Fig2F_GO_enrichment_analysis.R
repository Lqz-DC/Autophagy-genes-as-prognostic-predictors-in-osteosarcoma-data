library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

rt=read.table("id.txt",sep="\t",header=T,check.names=F)           
rt=rt[is.na(rt[,"entrezID"])==F,]                                 
gene=rt$entrezID

#GO enrichment analysis
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)               


pdf(file="barplot.pdf",width = 10,height = 8)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
