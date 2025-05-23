library(tidyverse)
library(dplyr)
results <- read.table("risk-violin.txt",header = T,sep = "\t",quote = "",row.names = 1)
risk_info <- read.table("risk.txt",header = T,sep = "\t",quote = "")
#Process data
data_long<-results%>%mutate(Sample=rownames(results))%>%gather(key="Cell_Type",value="Abundance",-Sample)#Convert to a long data frame
data_long<-merge(data_long,risk_info,by.x='Sample',by.y='id',all.x=T)

p_values<-data_long%>%group_by(Cell_Type)%>%summarise(p_value=wilcox.test(Abundance~risk)$p.value)
#Add the p-value to the data frame
data_long<-data_long%>%left_join(p_values,by="Cell_Type")%>%
                mutate(Significance=case_when(
                      p_value<0.001~"***",
                      p_value<0.01~"**",
                      p_value<0.05~"*",
                      TRUE~"ns"))
#Draw a box plot
plot<-ggplot(data_long,aes(x=Cell_Type,y=Abundance,fill=risk,color=risk))+
  geom_boxplot(fill=NA,outlier.shape=NA,size=1.2)+
  labs(y="",x="",title="Immune Cell Abundance",color='')+
  theme_bw()+
  scale_fill_manual(values=c("high" = "#BB2222", "low" = "skyblue"))+
  scale_color_manual(values=c("high" = "#BB2222", "low" = "skyblue"))+
  #Add significance labels for p-values
  geom_text(aes(label=Significance,y=max(Abundance)*1.1),size=5,color='black',show.legend=F)+
  theme(plot.title=element_text(hjust=0.5,size=30,face='bold'),
        legend.text=element_text(size=20),
        plot.margin=unit(c(3,3,3,3),'cm'),
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,size=20,color='black'),
        axis.text.y=element_text(size=20,color='black'))#Rotate the x-axis labels to prevent overlap
plot
ggsave(plot=plot,filename='diff.png',width=12,height=8,dpi=300)
ggsave(plot=plot,filename='risk——boxplot.pdf',width=12,height=8)
