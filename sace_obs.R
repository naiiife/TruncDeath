library(tbd)
library(MASS)
expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# Estimation Function

sace_obs <- function(Z,S,Y,X=NULL,V, subset=NULL, r0=NULL,
                     group='MNN', link='linear', sensitivity=0){
  ## robust=c('AIPW','regression','DR')
  ## link=c('iss','ss','oss','ross',linear')
  V = as.matrix(V)
  if (!is.null(subset)){
    Z = Z[subset]
    S = S[subset]
    Y = Y[subset]
    V = as.matrix(V[subset,])
    if (!is.null(X)) X = as.matrix(X)[subset,]
  }
  N = length(Z)
  Nt = which(Z==1)
  Nc = which(Z==0)
  N1 = sum(Z==1)
  N0 = sum(Z==0)
  Y[S==0] = 0
  if (is.null(r0)) {
    r0 = 1.5/log(N) #0.5/log(log(N)) #2.5/N^(1/3)
  }
  q = dim(V)[2]
  p = dim(X)[2]
  R = rep(1,N)
  if (!is.null(X)){
    for (j in 1:q){
      vj = V[,j]
      cv = lm(vj~.,data=data.frame(vj,X))
      R = R*(vj-predict(cv,newdata=data.frame(vj,X)))
    }
  } else {
    for (j in 1:q){
      vj = V[,j]
      R = R*(vj-mean(vj))
    }
  }
  V = as.numeric(R)
  s1 = glm(S~., family='binomial', data=data.frame(cbind(S,V,X))[Z==1,])
  betav = coef(s1)['V']
  rho = exp(betav*sensitivity*V)
  fs1 = as.numeric(predict(s1, newdata=data.frame(cbind(S,V,X)), type='response'))
  s0 = glm(S~., family='binomial', data=data.frame(cbind(S,V,X))[Z==0,])
  fs0 = as.numeric(predict(s0, newdata=data.frame(cbind(S,V,X)), type='response'))
  e = glm(Z~., family='binomial', data=data.frame(cbind(Z,V,X)))
  fe = as.numeric(predict(e, newdata=data.frame(cbind(Z,V,X)), type='response'))
  
  #fs1X = fs1*X
  m1 = glm(Y~., data=data.frame(cbind(Y,V,X))[Z==1&S==1,])
  fm1 = as.numeric(predict(m1,newdata=data.frame(cbind(Y,V,X)),type='response'))
  m0 = glm(Y~., data=data.frame(cbind(Y,V,X))[Z==0&S==1,])
  fm0 = as.numeric(predict(m0,newdata=data.frame(cbind(Y,V,X)),type='response'))
  
  Nt0 = 1:N
  neighbor = 1:N
  Ew = fm1 - fm0
  fy1 = Z*S*(Y-fm1)/fe/fs1
  fy0 = (1-Z)*S*(Y-fm0)/(1-fe)/fs0
  fw = Ew + fy1
  fwdr = Ew + fy1 - fy0
  w = v = wp = rep(1,N)
  bw0 = (1-fe)*fs0
  bw = bw0
  
  X = as.matrix(X)
  if (!is.null(X)) CX = solve(as.matrix(var(X)))
  for (i in 1:N){
    if (!is.null(X)){
    r = r0
    Xi = X - X + rep(1,N)%*%t(X[i,])
    M = diag((X-Xi)%*%CX%*%t(X-Xi))
    M = M/mean(M)
    for (j in 1:100){
      neighbor = which(M<r)
      if (length(neighbor)>10 & sd(V[neighbor])>0) break
      r = r*1.2
    }
    } else {
      neighbor = 1:N
    }
    fs1n = fs1[neighbor]
    rhon = rho[neighbor]
    Vn = V[neighbor]
    if (link=='iss') wup = wu = exp(betav*Vn)*rhon*(1/fs1n-mean(na.omit(1/fs1n)))
    if (link=='ss') wup = wu = exp(betav*Vn)*rhon*(fs1n-mean(na.omit(fs1n)))
    if (link=='oss') wup = wu = exp(betav*Vn)*rhon*(fs1n/(1-fs1n)-mean(na.omit(fs1n/(1-fs1n))))
    if (link=='ross') wup = wu = exp(betav*Vn)*rhon*(sqrt(fs1n/(1-fs1n))-mean(na.omit(sqrt(fs1n/(1-fs1n)))))
    if (link=='linear') {
      wu = exp(betav*Vn)*rhon*(Vn-mean(Vn))
      wup = exp(betav*Vn)*rhon*Vn*(Vn-mean(Vn))
    }
    w[i] = wu[neighbor==i]
    wp[i] = mean(wup)
    v[i] = mean(wu)
    bw[i] = mean(bw0[neighbor])
  }
  Gam1 = mean(w*(V*v-wp)/v^2*bw/mean((1-Z)*S)*(fm1-fm0))^2*summary(s1)$coefficients['V',2]^2
  w = w/v
  Rm = w*bw/mean((1-Z)*S)-(1-fe)*fs0/mean((1-Z)*S)
  G = t(Rm)%*%cbind(1,V,X)/N
  Gam0 = as.numeric(G%*%(summary(m0)$cov.unscaled*mean(residuals(m0)^2))%*%t(G))
  
  fY = w*fw*bw - (1-Z)*S*(Y-fm0)
  fYreg = w*Ew*bw
  n = length(fY)
  RD = mean(fY)/mean(S*(1-Z))
  vr = var(fY)/mean(S*(1-Z))^2/N
  vrc = vr + Gam1 + Gam0
  RDreg = mean(fYreg)/mean(S*(1-Z))
  vrreg = var(fYreg)/mean(S*(1-Z))^2/N
  
  return(list(sace=RD, se=sqrt(vr), sec=sqrt(vrc),
              sacereg=RDreg, sereg=sqrt(vrreg)))
}
