library(tbd)
library(MASS)
expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))


g = c(-1,1,0,0)
a = c(0,-1,0,1,1)
b0 = c(-2,2,2,2,4)
b1 = c(2,-2,-2,-2,4)
c0 = c(-1,-1,1,-1,0)
c1 = c(-3,1,1,1,0)
d0 = c(0,1,2,2,-2,0,-2)
d1 = c(-1,1,2,2,-2,3,1)
RD = sum(d1-d0)

# Setting 1 General
cat('Setting 1 General \n')
N = 200
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 500
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 2000
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

# Setting 2 RCT
cat('Setting 2 RCT \n')
b1 = b0
c1 = c0
d0 = c(0,1,2,2,0,0,-2)
d1 = c(-1,1,2,2,0,3,1)
RD = sum(d1-d0)

N = 200
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 500
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 2000
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

# Setting 3 Principal ign
cat('Setting 3 Principal ign \n')
b0 = c(-2,2,2,2,4)
b1 = c(2,-2,-2,-2,4)
c0 = c(-1,-1,1,-1,0)
c1 = c(-3,1,1,1,0)
d0 = c(0,1,2,2,0,0,0)
d1 = c(-1,1,2,2,0,0,1)
RD = sum(d1-d0)

N = 200
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 500
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 2000
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

# Setting 4, General, Heterogeneous
cat('Setting 4 Heterogeneous \n')
generatedata <- function(size,g,a,b0,b1,c0,c1,d0,d1,seeds){
  set.seed(seeds)
  # Covariates simulated from a Multivariate Normal Distribution
  mu = c(1, 1)
  Sigma = matrix(c(1,0.5,0.5,1),2,2)
  X12 = mvrnorm(size, mu, Sigma)
  X3 = runif(size, min = 0, max = 2)
  X = cbind(rep(1, size), X12, X3)
  #X1 = rbinom(size, 1, 0.5)
  #X2 = rbinom(size, 1, 0.4)
  #X3 = rbinom(size, 1, 0.6)
  #X = cbind(rep(1, size), X1, X2, X3)
  colnames(X) = c('X0', 'X1', 'X2', 'X3')
  V = rnorm(size, X%*%g, 2)
  #V = rbinom(size, 1, expit(X%*%g/2))
  
  Z = rbinom(size, 1, expit(cbind(X,V)%*%a))
  #Z = as.numeric(cbind(X,V)%*%a+rnorm(N,0,1)>0)
  Nc = (Z==0)
  Nt = (Z==1)
  
  S0 = rep(0, size)
  S1 = rep(0, size)
  S0[Nc] = rbinom(sum(Nc), 1, expit(cbind(X,V)[Nc,]%*%b0))
  S1[Nc] = S0[Nc] + (1-S0[Nc])*rbinom(sum(Nc), 1, expit(cbind(X,V)[Nc,]%*%c0))
  S0[Nt] = rbinom(sum(Nt), 1, expit(cbind(X,V)[Nt,]%*%b1))
  S1[Nt] = S0[Nt] + (1-S0[Nt])*rbinom(sum(Nt), 1, expit(cbind(X,V)[Nt,]%*%c1))
  S = Z*S1 + (1-Z)*S0
  
  Y0 = rnorm(size, cbind(X,V,S0,S1)%*%d0, 1)
  Y1 = rnorm(size, cbind(X,V,S0,S1)%*%d1, 1) + 2*X[,4]
  Y = Z*Y1 + (1-Z)*Y0
  Y[S==0] = 0
  TrueRD = mean((Y1-Y0)[Z==0&S0==1])
  
  return(list(Z=Z,S=S,Y=Y,X=X[,-1],V=V,S1=S1,S0=S0,Y1=Y1,Y0=Y0,TrueRD=TrueRD))
}

g = c(-1,1,0,0)
a = c(0,-1,0,1,1)
b0 = c(-2,2,2,2,4)
b1 = c(2,-2,-2,-2,4)
c0 = c(-1,-1,1,-1,0)
c1 = c(-3,1,1,1,0)
d0 = c(0,1,2,2,0,0,-2)
d1 = c(-1,1,2,2,0,3,1)
RD = sum(d1-d0)

N = 200
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 500
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))

N = 2000
B = 1000
sace.dr = rep(NA,B)
sace.aipw = rep(NA,B)
sace.reg = rep(NA,B)
se.dr = rep(NA,B)
se.aipw = rep(NA,B)
se.reg = rep(NA,B)
RDb = rep(RD,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RDb[b] = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  sace.dr[b] = fit$sace
  sace.aipw[b] = fit$sace
  sace.reg[b] = fit$sacereg
  se.dr[b] = fit$sec
  se.aipw[b] = fit$se
  se.reg[b] = fit$sereg
  
}
RD = mean(RDb)
cv.dr = mean(abs(sace.dr-RDb)<1.96*se.dr)
cv.aipw = mean(abs(sace.aipw-RDb)<1.96*se.aipw)
cv.reg = mean(abs(sace.reg-RDb)<1.96*se.reg)
wd.dr = mean(2*1.96*se.dr)
wd.aipw = mean(2*1.96*se.aipw)
wd.reg = mean(2*1.96*se.reg)
cvbt.dr = mean(abs(sace.dr-RDb)<1.96*sd(sace.dr))
cvbt.aipw = mean(abs(sace.aipw-RDb)<1.96*sd(sace.aipw))
cvbt.reg = mean(abs(sace.reg-RDb)<1.96*sd(sace.reg))
wdbt.dr = mean(2*1.96*sd(sace.dr))
wdbt.aipw = mean(2*1.96*sd(sace.aipw))
wdbt.reg = mean(2*1.96*sd(sace.reg))
c(mean(cv.dr),mean(cv.aipw),mean(cv.reg))
c(mean(wd.dr),mean(wd.aipw),mean(wd.reg))
c(mean(cvbt.dr),mean(cvbt.aipw),mean(cvbt.reg))
c(mean(wdbt.dr),mean(wdbt.aipw),mean(wdbt.reg))
