## Sensitivity analysis for monotonicity

generatedata <- function(size,g,a,b0,b1,c0,c1,d0,d1,varpi=1,seeds){
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
  S1[S0+S1==2] = rbinom(sum(S0+S1==2),1,varpi)
  S = Z*S1 + (1-Z)*S0
  
  Y0 = rnorm(size, cbind(X,V,S0,S1)%*%d0, 1)
  Y1 = rnorm(size, cbind(X,V,S0,S1)%*%d1, 1)
  Y = Z*Y1 + (1-Z)*Y0
  Y[S==0] = 0
  TrueRD = mean((Y1-Y0)[Z==0&S0==1&S1==1])
  
  return(list(Z=Z,S=S,Y=Y,X=X[,-1],V=V,S1=S1,S0=S0,Y1=Y1,Y0=Y0,TrueRD=TrueRD))
}

# Estimation

g = c(-1,1,0,0)
a = c(0,-1,0,1,1)
b0 = c(-2,2,2,2,4)
b1 = c(2,-2,-2,-2,4)
c0 = c(-1,-1,1,-1,0)
c1 = c(-3,1,1,1,0)
d0 = c(0,1,2,2,0,0,-2)
d1 = c(-1,1,2,2,0,3,1)
RD = sum(d1-d0)

N = 500
B = 100
meanbias.aipw = meanbias.reg = meanbias.sc = meanbias.wzr = NULL
mednbias.aipw = mednbias.reg = mednbias.sc = mednbias.wzr = NULL
for (vp in seq(1,0.5,length=11)){
  bias.aipw = rep(0,B)
  bias.reg = rep(0,B)
  bias.sc = rep(0,B)
  bias.wzr = rep(0,B)
  for (b in 1:B){
    dat = generatedata(N,g,a,b0,b1,c0,c1,d0,d1,varpi=vp,seeds=b)
    Z = dat$Z
    X = dat$X
    S = dat$S
    Y = dat$Y
    V = dat$V
    RD = dat$TrueRD
    fit = sace_obs(Z,S,Y,X,V,group='MNN',link='linear')
    bias.aipw[b] = fit$sace - RD
    bias.reg[b] = fit$sacereg - RD
    bias.sc[b] = mean(Y[S==1&Z==1]) - mean(Y[S==1&Z==0]) - RD
    bias.wzr[b] = tbd::sace(Z,S,Y,X,V,need.variance=FALSE)$sace - RD
    #cat(b,' ')
  }
  meanbias.aipw = append(meanbias.aipw, mean(bias.aipw))
  meanbias.reg = append(meanbias.reg, mean(bias.reg))
  meanbias.sc = append(meanbias.sc, mean(bias.sc))
  meanbias.wzr = append(meanbias.wzr, mean(na.omit(bias.wzr)))
  mednbias.aipw = append(mednbias.aipw, median(bias.aipw))
  mednbias.reg = append(mednbias.reg, median(bias.reg))
  mednbias.sc = append(mednbias.sc, median(bias.sc))
  mednbias.wzr = append(mednbias.wzr, median(na.omit(bias.wzr)))
  print(c(vp,mean(bias.aipw),median(bias.aipw)))
}
plot(seq(1,0.55,length=10), meanbias.aipw[1:10], type='l',col='blue', lty=1,ylab=c('Median bias'),
     xlab=expression(kappa),ylim=c(-4,4),lwd=1.5)
points(seq(1,0.55,length=10), meanbias.reg[1:10], type='l',col='blue', lty=2,lwd=1.5)
points(seq(1,0.55,length=10), meanbias.wzr[1:10], type='l',col='red', lty=1,lwd=1.5)
points(seq(1,0.55,length=10), meanbias.sc[1:10], type='l',col='black', lty=1,lwd=1.5)
abline(h=0,lty=2)
legend('topright',legend=c('AIPW','REG','WZR','SC'),col=c('blue','blue','red','black'),
       lty=c(1,2,1,1),cex=0.8)
