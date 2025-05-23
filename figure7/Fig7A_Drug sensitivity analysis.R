
library(limma)
library(oncoPredict)
library(parallel)
library(limma)
library(ggplot2)
library(ggpubr)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)


#Random number seed
set.seed(999)

#Read the input file
data=read.table("exp.txt", header=T, sep="\t", check.names=F,row.names = 1)
#Convert to matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#Read the reference file
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

#Conduct a drug sensitivity analysis
calcPhenotype(trainingExprData = GDSC2_Expr,    
              trainingPtype = GDSC2_Res,        
              testExprData = data,              
              batchCorrect = 'eb',  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,    
              minNumSamples = 10,    
              printOutput = TRUE,    
              removeLowVaringGenesFrom = 'rawData')

#Read in the drug sensitivity file
senstivity=read.csv("calcPhenotype_Output/DrugPredictions.csv", header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))

#Read in the risk file
risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)

#Merge
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
senstivity[is.na(senstivity)] = 0
senstivity = log2(senstivity+1)
rt=cbind(risk, senstivity)

#Set up comparison groups
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#Extract the drugs with significant differences
sigGene=c()
for(i in colnames(rt)[2:(ncol(rt))]){
  if(sd(rt[,i])<0.05){next}
  wilcoxTest=wilcox.test(rt[,i] ~ rt[,"risk"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.001){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "risk")
rt=rt[,sigGene]

#Convert the data into a file format compatible with ggplot2
rt=melt(rt,id.vars=c("risk"))
colnames(rt)=c("risk","Gene","Expression")

#Set up comparison groups
group=levels(factor(rt$risk))
rt$risk=factor(rt$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

#Draw a box plot
boxplot=ggboxplot(rt, x="Gene", y="Expression", fill="risk",
                  xlab="",
                  ylab="Drug Senstivity",
                  legend.title="Risk",
                  width=0.8,
                  palette = c("DodgerBlue1","Firebrick2") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
  )

#Output image
pdf(file="drugSenstivity.pdf", width=20, height=8)
print(boxplot)
dev.off()
