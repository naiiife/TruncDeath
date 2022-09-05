library(tbd)
library(MASS)
expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# Simulation

generatedata <- function(size,g,a,b0,b1,c0,c1,d0,d1,seeds){
  set.seed(seeds)
  # Covariates simulated from a Multivariate Normal Distribution
  mu = c(1, 1)
  Sigma = matrix(c(1,0.5,0.5,1),2,2)
  X12 = mvrnorm(size, mu, Sigma)
  X3 = runif(size, min = 0, max = 2)
  X = cbind(rep(1, size), X12, X3)
  colnames(X) = c('X0', 'X1', 'X2', 'X3')
  V = rnorm(size, X%*%g, 2)
  
  Z = rbinom(size, 1, expit(cbind(X,V)%*%a))
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
  Y1 = rnorm(size, cbind(X,V,S0,S1)%*%d1, 1)
  Y = Z*Y1 + (1-Z)*Y0
  Y[S==0] = 0
  TrueRD = mean((Y1-Y0)[Z==0&S0==1])
  
  return(list(Z=Z,S=S,Y=Y,X=X[,-1],V=V,S1=S1,S0=S0,Y1=Y1,Y0=Y0,TrueRD=TrueRD))
}

# Estimation

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

N = 200
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 1, General, N=200\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 500
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 1, General, N=500\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 2000
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 1, General, N=2000\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

# Setting 2 RCT

b1 = b0
c1 = c0
d0 = c(0,1,2,2,0,0,-2)
d1 = c(-1,1,2,2,0,3,1)
RD = sum(d1-d0)

N = 200
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 2, RCT, N=200\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 500
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 2, RCT, N=500\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 2000
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 2, RCT, N=2000\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

# Setting 3 Principal ign

b0 = c(-2,2,2,2,4)
b1 = c(2,-2,-2,-2,4)
c0 = c(-1,-1,1,-1,0)
c1 = c(-3,1,1,1,0)
d0 = c(0,1,2,2,0,0,0)
d1 = c(-1,1,2,2,0,0,1)
RD = sum(d1-d0)

N = 200
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 3, Prinign, N=200\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 500
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 3, Prinign, N=500\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 2000
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 3, Prinign, N=2000\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

# Setting 4, General, Heterogeneous

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
b1 = b0
c0 = c(-1,-1,1,-1,0)
c1 = c0
d0 = c(0,1,2,2,-2,0,-2)
d1 = c(-1,1,2,2,-2,3,1)
RD = sum(d1-d0)

N = 200
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 4, hetero, N=200\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 500
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 4, hetero, N=500\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

N = 2000
B = 1000
bias.aipw = rep(0,B)
bias.reg = rep(0,B)
bias.dr = rep(0,B)
bias.sc = rep(0,B)
bias.wzr = rep(0,B)
for (b in 1:B){
  dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
  Z = dat$Z
  X = dat$X
  S = dat$S
  Y = dat$Y
  V = dat$V
  RD = dat$TrueRD
  fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
  bias.aipw[b] = fit$sace - RD
  bias.reg[b] = fit$sacereg - RD
  bias.dr[b] = fit$sacedr - RD
  bias.sc[b] = coef(lm(Y~Z+V+X,subset=(S==1)))['Z'] - RD
  bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  #cat(b,' ')
}
cat('Setting 4, hetero, N=2000\n')

c(mean(bias.aipw), median(bias.aipw), sqrt(mean(bias.aipw^2)))
c(mean(bias.reg), median(bias.reg), sqrt(mean(bias.reg^2)))
c(mean(na.omit(bias.wzr)), median(na.omit(bias.wzr)), sqrt(mean(na.omit(bias.wzr^2))))
c(mean(bias.sc), median(bias.sc), sqrt(mean(bias.sc^2)))
boxplot(na.omit(cbind(bias.aipw, bias.reg, bias.wzr, bias.sc)),ylim=c(-6,3))
abline(h=0, lty=2)

# Sensitivity analysis

g = c(-1,1,0,0)
a = c(0,-1,0,1,1)
b0 = c(-2,2,2,2,4)
b1 = c(2,-2,-2,-2,4)
c0 = c(-1,-1,1,-1,0)
c1 = c(-3,1,1,1,-4)
d0 = c(0,1,2,2,0,0,-2)
d1 = c(-1,1,2,2,0,3,1)
RD = sum(d1-d0)

N = 500
B = 1000
R = 39
bias.aipw.s = rep(0,R)
bias.reg.s = rep(0,R)
bias.sc.s = rep(0,R)
bias.wzr.s = rep(0,R)
for (r in 1:R){
  bias.aipw = rep(0,B)
  bias.reg = rep(0,B)
  bias.sc = rep(0,B)
  bias.wzr = rep(0,B)
  for (b in 1:B){
    dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
    Z = dat$Z
    X = dat$X
    S = dat$S
    Y = dat$Y
    V = dat$V
    RD = dat$TrueRD
    rr = r/20-1
    fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear',sensitivity=-log(1-rr)/log(2))
    bias.aipw[b] = fit$sace - RD
    bias.reg[b] = fit$sacereg - RD
    bias.sc[b] = mean(Y[S==1&Z==1]) - mean(Y[S==1&Z==0]) - RD
    bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  }
  bias.aipw.s[r] = mean(bias.aipw)
  bias.reg.s[r] = mean(bias.reg)
  bias.wzr.s[r] = mean(na.omit(bias.wzr))
  bias.sc.s[r] = mean(bias.sc)
  cat(r,'')
}
plot(1:R/20-1,bias.sc.s,type='l',lwd=2,ylim=c(-5,5),ylab='bias',
     xlab=expression(eta),main=expression(gamma==-4))
points(1:R/20-1,bias.wzr.s,type='l',lwd=2,col='red')
points(1:R/20-1,bias.aipw.s,type='l',lwd=2,col='blue')
points(1:R/20-1,bias.reg.s,type='l',lwd=2,lty=2,col='blue')
abline(h=0,lty=2)
legend(col=c('blue','blue','red','black'),lwd=c(2,2,2,2),lty=c(1,2,1,1),
       legend=c('AIPW','REG','WZR','SC'),'topright',cex=0.8)

c1 = c(-3,1,1,1,0)
bias.aipw.s = rep(0,R)
bias.reg.s = rep(0,R)
bias.sc.s = rep(0,R)
bias.wzr.s = rep(0,R)
for (r in 1:R){
  bias.aipw = rep(0,B)
  bias.reg = rep(0,B)
  bias.sc = rep(0,B)
  bias.wzr = rep(0,B)
  for (b in 1:B){
    dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
    Z = dat$Z
    X = dat$X
    S = dat$S
    Y = dat$Y
    V = dat$V
    RD = dat$TrueRD
    rr = r/20-1
    fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear',sensitivity=-log(1-rr)/log(2))
    bias.aipw[b] = fit$sace - RD
    bias.reg[b] = fit$sacereg - RD
    bias.sc[b] = mean(Y[S==1&Z==1]) - mean(Y[S==1&Z==0]) - RD
    bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  }
  bias.aipw.s[r] = mean(bias.aipw)
  bias.reg.s[r] = mean(bias.reg)
  bias.wzr.s[r] = mean(na.omit(bias.wzr))
  bias.sc.s[r] = mean(bias.sc)
  cat(r,'')
}
plot(1:R/20-1,bias.sc.s,type='l',lwd=2,ylim=c(-5,5),ylab='bias',
     xlab=expression(eta),main=expression(gamma==0))
points(1:R/20-1,bias.wzr.s,type='l',lwd=2,col='red')
points(1:R/20-1,bias.aipw.s,type='l',lwd=2,col='blue')
points(1:R/20-1,bias.reg.s,type='l',lwd=2,lty=2,col='blue')
abline(h=0,lty=2)
legend(col=c('blue','blue','red','black'),lwd=c(2,2,2,2),lty=c(1,2,1,1),
       legend=c('AIPW','REG','WZR','SC'),'topright',cex=0.8)

c1 = c(-3,1,1,1,4)
bias.aipw.s = rep(0,R)
bias.reg.s = rep(0,R)
bias.sc.s = rep(0,R)
bias.wzr.s = rep(0,R)
for (r in 1:R){
  bias.aipw = rep(0,B)
  bias.reg = rep(0,B)
  bias.sc = rep(0,B)
  bias.wzr = rep(0,B)
  for (b in 1:B){
    dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,seeds=b)
    Z = dat$Z
    X = dat$X
    S = dat$S
    Y = dat$Y
    V = dat$V
    RD = dat$TrueRD
    rr = r/20-1
    fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear',sensitivity=-log(1-rr)/log(2))
    bias.aipw[b] = fit$sace - RD
    bias.reg[b] = fit$sacereg - RD
    bias.sc[b] = mean(Y[S==1&Z==1]) - mean(Y[S==1&Z==0]) - RD
    bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
  }
  bias.aipw.s[r] = mean(bias.aipw)
  bias.reg.s[r] = mean(bias.reg)
  bias.wzr.s[r] = mean(na.omit(bias.wzr))
  bias.sc.s[r] = mean(bias.sc)
  cat(r,'')
}
plot(1:R/20-1,bias.sc.s,type='l',lwd=2,ylim=c(-5,5),ylab='bias',
     xlab=expression(eta),main=expression(gamma==4))
points(1:R/20-1,bias.wzr.s,type='l',lwd=2,col='red')
points(1:R/20-1,bias.aipw.s,type='l',lwd=2,col='blue')
points(1:R/20-1,bias.reg.s,type='l',lwd=2,lty=2,col='blue')
abline(h=0,lty=2)
legend(col=c('blue','blue','red','black'),lwd=c(2,2,2,2),lty=c(1,2,1,1),
       legend=c('AIPW','REG','WZR','SC'),'topright',cex=0.8)
