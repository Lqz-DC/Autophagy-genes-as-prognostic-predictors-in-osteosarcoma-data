library(survival)
library(rms)
library(Hmisc)
library(rmda)
library(regplot)
library(ggplot2)
library(survcomp)

#Read the risk input file
rt=read.table("all.txt", header=T, sep="\t", check.names=F, row.names=1)

dd <- datadist(rt)
options(datadist = "dd")


#合并数据
colnames(rt)



# 模型构建

fit <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, data = rt)              # BM




# fit
pred <- predict(fit, type = "risk")
c2 <- concordance.index(x = pred2, surv.time = rt$futime, surv.event = rt$fustat, method = "noether")
print(c2$c.index)
print(c2$lower)
print(c2$upper)

pdf(file="BM_Nomogram.pdf", width=6, height=6)
nom=regplot(fit,
             clickable=F,
             title="",
             points=TRUE,
             droplines=T,
             observation=NULL,
             rank=NULL,
             failtime = c(1,3,5),
             showP = F,
             prfail = F) 
dev.off()
##### calibration curves####
####FMcalibration####
#1 year
pdf(file="FMcalibration.pdf", width=6, height=6)
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, x=T, y=T, surv=T, data=rt1, time.inc=1)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=1, m=28,B=1000)##  m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, x=T, y=T, surv=T, data=rt1, time.inc=3)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, x=T, y=T, surv=T, data=rt1, time.inc=5)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()



