library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
#Prepare the data including gene expression and the IC values of the drugs you wish to study
dd1 <- read.table("input.txt",row.names = 1,header = T,sep = "\t")
gene= "MYC"#PEA15\SAR1A\MYC\BNIP3
####Vorinostat\Trametinib\Selumetinib\Ribociclib

##Grouping: The data were divided into the high-expression group (high) and the low-expression group (low) based on the median value of gene expression levels.
dd2<-dd1%>% mutate(Group=ifelse(dd1[,gene]>median(dd1[,gene]),'high','low'),.before=1)
table(dd2$Group)

####Draw a box plot####
pdf("Vorinostat-MYC.pdf",width = 5,height = 4)
ggplot(dd2,mapping=aes(x=Group,y=Vorinostat))+
  stat_boxplot(geom='errorbar',width=0.3,
               position=position_dodge(0.75)
  )+
  geom_boxplot(mapping=aes(fill=Group),position=position_dodge(0.75),size=0.8,outliers = F#outlier.shape = NA
  )+
  stat_compare_means(method='t.test',label ='p.format',
                     paired = F,
                     label.x.npc='centre',label.y.npc='top')+
  geom_jitter(mapping=aes(color=Group),width=0.1,size=2,height=0)+
  labs(y='Vorinostat sensitivity',x=gene)+
  guides(fill='none',color='none')+
  theme_bw()
dev.off()
