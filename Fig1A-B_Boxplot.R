
library(ggpubr)           
inputFile="TARGET_clinical.txt"#input file
inputFile="GEO_clinical.txt"
outFile="TARGET-boxplot.pdf"      #output file
outFile="GEO-boxplot.pdf"


#reading a file
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
x=colnames(rt)[2]
y=colnames(rt)[3]
colnames(rt)=c("id","Group","Expression")

#Set up a comparison group
group=levels(factor(rt$Group))
rt$Group=factor(rt$Group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#Draw a boxplot
boxplot=ggboxplot(rt, x="Group", y="Expression", color="Group",
                  xlab=x,
                  ylab=y,
                  legend.title=x,
                  palette = c("#F1515E","#1DBDE6"),
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons,method =  "wilcox.test")+labs(x='',y='Futime(years)',title = "TARGET")
  #stat_compare_means(comparisons = my_comparisons,method =  "wilcox.test")+labs(x='',y='Futime(years)',title = "GES21257")
#输出图片
pdf(file=outFile,width=6,height=5)
print(boxplot)
dev.off()
