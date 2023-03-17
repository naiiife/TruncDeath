# Sensitivity analysis (nondifferential substitution)

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
    fit = sace_obs(Z,S,Y,X,V,link='linear',sensitivity=-log(1-rr)/log(2))
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
    fit = sace_obs(Z,S,Y,X,V,link='linear',sensitivity=-log(1-rr)/log(2))
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
    fit = sace_obs(Z,S,Y,X,V,link='linear',sensitivity=-log(1-rr)/log(2))
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
