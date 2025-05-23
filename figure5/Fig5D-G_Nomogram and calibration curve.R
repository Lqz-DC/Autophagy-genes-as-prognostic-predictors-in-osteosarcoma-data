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
rt1=rt[,-c(6,7,8,9,11)]#Metastasis,Age,Gender and risk Score
rt2=rt1[,-6]#Metastasis,Age and Gender
rt3 <- rt1[,-5]#Age,Gender and risk Score


# 模型构建
fit1 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group + riskScore, data = rt1)  # FM
fit2 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, data = rt2)              # BM
fit3 <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, data = rt3)          # PM




# fit1
pred1 <- predict(fit1, type = "risk")
c1 <- concordance.index(x = pred1, surv.time = rt1$futime, surv.event = rt1$fustat, method = "noether")
print(c1$c.index)
print(c1$lower)
print(c1$upper)

# fit2
pred2 <- predict(fit2, type = "risk")
c2 <- concordance.index(x = pred2, surv.time = rt2$futime, surv.event = rt2$fustat, method = "noether")
print(c2$c.index)
print(c2$lower)
print(c2$upper)
# fit3
pred3 <- predict(fit3, type = "risk")
c3 <- concordance.index(x = pred3, surv.time = rt3$futime, surv.event = rt3$fustat, method = "noether")
print(c3$c.index)
print(c3$lower)
print(c3$upper)
#Compare the merits of different models using the Likelihood Ratio Test
anova(fit2, fit1, test = "LRT")
anova(fit3, fit1, test = "LRT")
anova(fit2, fit3, test = "LRT")
AIC(fit1, fit2, fit3)


#Draw Nomograma Plot

nom1=regplot(fit1,
             clickable=F,
             title="",
             points=T,#Display the score line
             droplines=F,#Show the connection line
             observation=NULL,
             rank=NULL,
             failtime = c(1,3,5),
             showP = F,
             prfail = F) 

nom2=regplot(fit2,
             clickable=F,
             title="",
             points=TRUE,
             droplines=T,
             observation=NULL,
             rank=NULL,
             failtime = c(1,3,5),
             showP = F,
             prfail = F) 

nom3=regplot(fit3,
             clickable=F,
             title="",
             points=T,
             droplines=T,
             observation=NULL,
             rank=NULL,
             failtime = c(1,3,5),
             showP = F,
             prfail = F) 

#####Merge calibration curves####
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

##BMcalibration####
pdf(file="BMcalibration.pdf", width=6, height=6)
#1 year
fit2 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt2, time.inc=1)
cal <- calibrate(fit2, cmethod="KM", method="boot", u=1, m=28,B=1000)##   m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit2 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender + Group, x=T, y=T, surv=T, data=rt2, time.inc=3)
cal <- calibrate(fit2, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit2 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender + Group, x=T, y=T, surv=T, data=rt2, time.inc=5)
cal <- calibrate(fit2, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()

####PMcalibration####
pdf(file="PMcalibration.pdf", width=6, height=6)
#1 year
fit3 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, x=T, y=T, surv=T, data=rt3, time.inc=1)
cal <- calibrate(fit3, cmethod="KM", method="boot", u=1, m=28,B=1000)## m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit3 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, x=T, y=T, surv=T, data=rt3, time.inc=3)
cal <- calibrate(fit3, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit3 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+riskScore, x=T, y=T, surv=T, data=rt3, time.inc=5)
cal <- calibrate(fit3, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()

