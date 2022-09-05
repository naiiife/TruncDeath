## Data Analysis
library(tbd)
#setwd('D:/PAPER/SACE/codes')
dat = read.csv('leukemiaPKU.csv')
attach(dat)
Z <- as.numeric(dat$HLABHWD=='0')
RELAPSE <- as.numeric(dat$RELAPSE)
TRM <- as.numeric(dat$TRM)
N = length(Z)
V = dat$AGE
S = as.numeric(TRM==0|dat$TRMT>730)
Y = 1-(RELAPSE==0|dat$RELAPSET>730)
X1 = dat$MRD
X2 = dat$CR
X3 = dat$ALL
X = cbind(X1,X2,X3)
colnames(X) = c('MRD','CR','ALL')

## Estimation
tbd::sace(Z,S,Y,X,V)
sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
summary(lm(Y~Z+X+V,subset=(S==1)))

## Bootstrap
B = 1000
ate.aipw.b = rep(NA,B)
ate.reg.b = rep(NA,B)
for (b in 1:B){
  set.seed(b)
  X1 = dat$MRD
  X2 = dat$CR
  X3 = dat$ALL
  ss = sample(N,replace=TRUE)
  if (mean(X1[ss])==0 | mean(X1[ss])==1) X1 = NULL
  if (mean(X2[ss])==0 | mean(X2[ss])==1) X2 = NULL
  if (mean(X3[ss])==0 | mean(X3[ss])==1) X3 = NULL
  X0 = cbind(X1,X2,X3)
  if (is.null(X0)) next
  fit = sace_obs(Z,S,Y,X0,V,subset=ss,group='MNN',link='linear')
  ate.aipw.b[b] = fit$sace
  ate.reg.b[b] = fit$sacereg
  cat(b,'')
}
sd(ate.aipw.b)
sd(ate.reg.b)
boxplot(ate.aipw.b)

## Sensitivity analysis
V2 = as.numeric(V>27)
R = 39
ate.full = rep(0,R)
ate.dic = rep(0,R)
for (r in 1:R){
    rr = r/20-1
    fit.full = sace_obs(Z,S,Y,X,V,sensitivity=-log(1-rr)/log(2))
    fit.dic = sace_obs(Z,S,Y,X,V2,sensitivity=-log(1-rr)/log(2))
    ate.full[r] = fit.full$sace
    ate.dic[r] = fit.pos$sace
  cat(r,'')
}
plot(1:R/20-1,ate.full,type='l',lwd=2,col='black',ylim=c(-0.25,0.7),ylab='SACEC estimate',
     xlab=expression(eta),main='Sensitivity analysis')
points(1:R/20-1,ate.dic,type='l',lwd=2,lty=2,col='darkcyan')
abline(h=0,lty=2)
legend('topright',col=c('black','darkcyan'),lwd=c(2,2),lty=c(1,2),
       legend=c('Continuous V','Dichotomized V'),cex=0.8)
