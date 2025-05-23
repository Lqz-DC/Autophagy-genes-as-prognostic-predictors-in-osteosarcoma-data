rm(list=ls())
#加载R包
library(ggplot2)
library(tibble)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(patchwork)
library(clusterProfiler)
df = read.table( "input.txt",header = T,as.is = T)

df1<-df[order(df$foldChange, decreasing = T),] #Sort in descending order
rownames(df1) <- df1[,1]
df2 <- df1[,-1]
names(df2) <- rownames(df1)
head(df2)
####HALLMARK Gene Library####
library(msigdbr)
hallmark_gene_sets <- msigdbr(species ="Homo sapiens", category ="H")
geneset <- data.frame(term = hallmark_gene_sets$gs_name,
                      gene = hallmark_gene_sets$gene_symbol)
# Remove the "HALLMARK_" prefix to make the display more concise.
geneset$term <- gsub(pattern ="HALLMARK_","", geneset$term)
# Run GSEA analysis
egmt <- GSEA(df2,
             TERM2GENE = geneset,
             pvalueCutoff =0.05, 
             minGSSize =1,
             maxGSSize =500000)
head(df2)
head(geneset)

# 查看结果的前几行
head(egmt@result[,1:6])
# Extract the required data columns
data1 <- egmt@result[, c("ID","NES","setSize","pvalue")]
data <- egmt@result
# Sort by NES value
data <- data[order(data$NES, decreasing =TRUE), ]

write.table(data, "HALLMARK_GSEAresult.txt",sep = "\t",quote = F,col.names = T,row.names = F)

egmt <- read.table("HALLMARK_GSEAresult.txt",header = T,sep = "\t",quote = "",row.names =1 )
#Drawing plot
data<-egmt[,c("ID","NES","setSize","pvalue")]
data$setSize_1<-data$setSize/10
head(data)

#Set the name of the path on the Y-axis as "factor" and sort it.
data<-data[order(data$NES,decreasing=T),]
data$ID<-factor(data$ID,levels=data$ID)
data$xlab<-1:29
head(data)
summary(data$NES)
summary(data$setSize_1)

#The names of the pathways marked in the figure
label<-c("E2F_TARGETS","G2M_CHECKPOINT","MYC_TARGETS_V1","MYC_TARGETS_V2","CHOLESTEROL_HOMEOSTASIS",
         "KRAS_SIGNALING_DN","WNT_BETA_CATENIN_SIGNALING","GLYCOLYSIS","HYPOXIA",
         "INTERFERON_GAMMA_RESPONSE","INFLAMMATORY_RESPONSE","INTERFERON_ALPHA_RESPONSE",
         "IL6_JAK_STAT3_SIGNALING","EPITHELIAL_MESENCHYMAL_TRANSITION","KRAS_SIGNALING_UP",
         "TNFA_SIGNALING_VIA_NFKB","ANGIOGENESIS","P53_PATHWAY")

#Extract the corresponding pathway
data_label<-data[data$ID%in%label,]
data_label
#The color of the path in the picture
data_label$col<-c("#BA2121","#BA2121","#BA2121","#BA2121","#BA2121","#BA2121","#BA2121","#BA2121","#BA2121",
                  "#86CEEB","#86CEEB","#86CEEB","#86CEEB","#86CEEB","#86CEEB","#86CEEB","#86CEEB","#86CEEB")

p<-ggplot(data=data,aes(x=xlab,y=NES))+
  geom_point(aes(size=setSize_1,alpha=-log10(pvalue)),shape=21,stroke=0.7,fill="#0000ff",colour="black")+#stroke：Set the border width of the point.
  scale_size_continuous(range=c(0.2,6))+
  xlab(label="Hallmark gene sets")+
  ylab(label="Normalized enrichment score (NES)")+
  theme_classic(base_size=15)+
  scale_x_continuous(breaks=seq(0,50,by=10),labels=seq(0,50,by=10))+#Set the scale lines and labels for the x-axis
  scale_y_continuous(breaks=seq(-3,3,by=1),labels=seq(-3,3,by=1))+
  guides(size=guide_legend(title="setSize"),
         alpha=guide_legend(title="-log10(pvalue)"))+
  theme(
    axis.line=element_line(color="black",size=0.6),
    axis.text=element_text(face="bold"),
    axis.title=element_text(size=13)
  )

p


p3<-p+
  geom_text_repel(data=data_label,aes(x=xlab,y=NES,label=ID),size=3,color=data_label$col,
                  force=20,
                  point.padding=0.5,
                  min.segment.length=2,
                  hjust=1.2,
                  segment.color="grey20",
                  segment.size=0.3,
                  segment.alpha=1,
                  nudge_y=-0.2,
                  nudge_x = 0.1
  )
p3

#保存
ggsave(filename="HallmarkGSEA.pdf",plot=p3,width=9.2,height=4.6)

