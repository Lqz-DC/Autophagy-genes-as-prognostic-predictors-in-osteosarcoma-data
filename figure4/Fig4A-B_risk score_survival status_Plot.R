######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺2749657388@qq.com
######????΢??: 18520221056

#install.packages("pheatmap")

library(pheatmap)
setwd("C:\\Users\\ACER\\Desktop\\Autophagy\\13.risk")             #???ù???Ŀ¼
rt=read.table("riskTest.txt",sep="\t",header=T,row.names=1,check.names=F)       #??ȡ?????ļ?
rt=rt[order(rt$riskScore),]                                     #????riskScore????Ʒ????

#绘制风险曲线
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScore2.pdf",width = 10,height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Samples number",
     ylab="Risk score",
     col=c(rep("#4665D9",lowLength),
     rep("#F584F0",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "Low risk"),bty="n",pch=19,col=c("#F584F0","#4665D9"),cex=1.2)
dev.off()

#绘制生存状态图
color=as.vector(rt$fustat)
color[color==1]="#f8c828"
color[color==0]="#007e5d"
pdf(file="survStat.pdf",width = 10,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Samples number",
     ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#f8c828","#007e5d"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#绘制风险热图
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",width = 10,height = 4)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()

A <- rgb(245/255,132/255,240/255)
B <- rgb(70/255,101/255,217/255)
