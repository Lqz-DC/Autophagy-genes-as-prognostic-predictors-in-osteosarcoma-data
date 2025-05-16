
rm(list = ls())
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)

#input file
df <- read.csv("limmaTab.csv",header = T,row.names = 1)
head(df)

df$Group <- factor(ifelse(df$P.Value < 0.05 & abs(df$logFC) >= 0.2,
                          ifelse(df$logFC >= 0.2, 'Up','Down'),'Stable'))
df[1:10,1:7]

table(df$Group)
df$gene <- row.names(df)

p <- ggplot(df, aes(x = logFC, y = -log10(P.Value),,colour = Group))+
  geom_point( shape = 19, size=2.5,stroke = 0.5)+

  scale_color_manual(values=c( "#1874CD",'gray',"#CD2626"))+
  ylab('-log10 (Pvalue)')+
  xlab('log2 (Fold Change)')+
  labs(title = "No_metastases vs Metastases")+
  #The gene name of the added focus point
  geom_text_repel(
    data = df[df$P.Value < 0.05 & abs(df$logFC) > 0.2,],
    aes(label = gene),
    size = 3.5,
    segment.color = NA )+ 
  geom_vline(xintercept = c(-0.2,0.2),lty = 2, col = "black", lwd = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5)+
  theme_bw(
    base_line_size = 1  
  )+
  guides(fill = guide_legend(override.aes = list(size =3)))+
  theme_bw()+
  theme(
  axis.title.x = element_text(hjust = 0.5),
  legend.position = c(0.08, 0.86)
)


p + theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
  plot.title = element_text(
    hjust = 0.5, 
    size = 14,   
    face = "bold", 
    vjust = 1.5
  )
)



