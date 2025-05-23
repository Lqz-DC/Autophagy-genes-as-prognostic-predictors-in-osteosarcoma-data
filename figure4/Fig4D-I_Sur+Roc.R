
library(survival)
library(survminer)
library(timeROC)

inputFile="riskTest.txt"       
survFile="test-survival.pdf"        
rocFile="test-ROC.pdf"             

inputFile="riskTrain.txt"       
survFile="train-survival.pdf"        
rocFile="train-ROC.pdf"  


rt=read.table(inputFile,header=T,sep="\t")

diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}

#Sur Plot
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=T,
		           pval=pValue,
		           pval.size=5,
		           risk.table=TRUE,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c('#EBBA37',"#8AD293"),
		           risk.table.height=.25,
                title = "Merged Test")
pdf(file=survFile,onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()

###Roc Plot
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file=rocFile,width=7,height=5)
plot(ROC_rt,time=1,col="#8AD293",title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='#EBBA37',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col="#4665D9",add=TRUE,title=FALSE,lwd=2)
title(main = "Merged Test", cex.main = 1.5)
legend('bottomright',
        c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],3)),
          paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],3)),
          paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],3))),
        col=c("#8AD293",'#EBBA37',"#4665D9"),lwd=2,bty = 'n')
dev.off()

