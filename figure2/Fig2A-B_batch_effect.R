
library(sva)
library(limma)
library(FactoMineR)
library(factoextra)

rt=read.table("mergegroup.txt",sep="\t",header=T,check.names=F)
rt <- t(rt)
rt=as.matrix(rt)

exp <- rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

expr <- data
batch <- c(rep("GSE21257",53),rep("TARGET",86))
tissue <- c(rep("Non_metastasis",19),rep("Metastasis",34),rep("Non_metastasis",65),rep("Metastasis",21))
mode <- model.matrix(~as.factor(tissue))


####remove Batch Effect
limma_expr <- removeBatchEffect(expr,batch = batch,design = mode)

#PCA analysis without removing batch effects
pre.pca <- PCA(t(expr),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             title = "Before Batchremoving",
             col.ind = batch,
             addEllipses = TRUE,
             legend.title="Group"  )

#Batch Effect-Corrected PCA Analysis
combat.pca <- PCA(t(limma_expr),graph = FALSE)
fviz_pca_ind(combat.pca,
             geom= "point",
             title = "After Batchremoving",
             col.ind = batch,
             addEllipses = TRUE,
             legend.title="Group"  )
write.csv(limma_expr,"limma_expr.csv")
