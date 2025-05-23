#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")

#???ð?
library(survival)
library(survminer)
library(timeROC)

inputFile="riskTest.txt"        #?????ļ?
survFile="test-survival.pdf"         #?????????ļ?
rocFile="test-ROC.pdf"               #ROC?????ļ?

inputFile="riskTrain.txt"        #?????ļ?
survFile="train-survival.pdf"         #?????????ļ?
rocFile="train-ROC.pdf"  

#??ȡ?????ļ?
rt=read.table(inputFile,header=T,sep="\t")

#?Ƚϸߵͷ????????????죬?õ???????pֵ
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
#????????????
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

a=rgb(138/255,210/255,147/255)
b=rgb(235/255,186/255,55/255)
c=rgb(70/255,101/255,217/255)
###ROC????
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


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ???䣺seqbio@foxmail.com
######?⿡??ʦ΢??: seqBio

"High"="#F1515E"
"Low"="#1DBDE6"


"Up"="#f0a202"
"Down"="#008751"

"Metastases"="#b597f6"
"No_metastases"="#96c6ea"


"Alive"="#b0db43"
"Dead"="#db504a"


"Male"="#f9c58d"
"Female"="#f492f0"






