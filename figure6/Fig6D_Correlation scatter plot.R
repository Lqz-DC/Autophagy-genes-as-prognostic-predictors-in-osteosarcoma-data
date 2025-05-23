
library(ggplot2)
library(ggpubr)
library(ggExtra)

inputFile="input2.txt"      
gene1="BNIP3"             #Gene name
gene2="B cells memory"   #Names of immune cells
  

#Read the input file and extract the gene expression levels
rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
x=as.numeric(rt[gene1,])
y=as.numeric(rt[gene2,])

#Correlation analysis
df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="pearson")
cor=corT$estimate
pValue=corT$p.value
p1=ggplot(df1, aes(x, y)) + 
			xlab(gene1)+ylab(gene2)+
			geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
			stat_cor(method = 'pearson', aes(x =x, y =y))

pdf(file="BNIP3-B cells memory.pdf",width=5,height=4.8)
print(p1)
dev.off()

#Other genes' correlations with cells are input according to the above content