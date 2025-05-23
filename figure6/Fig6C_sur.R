
picDir="picture"                                              
dir.create(picDir)

library(survival)
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F)
rt$futime=rt$futime/365                                        
outTab=data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  if (TRUE & FALSE %in% a){
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,3)

  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)

  pdf(file=paste(picDir,"\\",gene,".survival.pdf",sep=""),
       width = 6,            
       height =5,            
       bg="white")
  plot(fit, 
     lwd=2,
     col=c("#BA2121", "#86CEEB"),
     xlab="Time (year)",
     mark.time=T,
     ylab="Survival rate",
     main=paste(gene,"(p=", pValue ,")",sep="") )
  legend("topright", 
       c("High","Low"), 
       lwd=2, 
       col=c("#BA2121", "#86CEEB"))
  dev.off()
  }
}
write.table(outTab,file="survival.xls",sep="\t",row.names=F,quote=F)
