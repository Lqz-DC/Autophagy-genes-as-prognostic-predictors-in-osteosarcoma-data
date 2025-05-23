
library(survival)
library(forestplot)
                   
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")            
rt=read.table("input.txt",header=T,sep="\t",check.names=F,row.names=1)

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)


rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )         
pdf(file="forest-multi.pdf",onefile = FALSE,
       width = 6,            
       height = 4,            
       )
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
           )
dev.off()
