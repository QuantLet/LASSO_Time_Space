[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **LTS_lambda_dep** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

﻿
Name of Quantlet: LTS_lambda_dep

Published in: LASSO_Time_Space

Description: 'Compares the prediction norm with the penalty level selected equation by equation or jointly over equations in a regression system with dependent data.'

Keywords : 
- simulation
- lasso shrinkage
- L1-norm penalty
- regression
- error

Author:          Chen Huang

Submitted:        Tuesday, June 19 2018 by Chen Huang


```

### R Code
```r

libraries = c("mvtnorm", "hdm", "Matrix", "sandwich", "matrixStats", "doSNOW")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

c1 = makeCluster(20) # core numbers
registerDoSNOW(c1)

nboot = 5000
rep = 1000
T     = 100 #sample size
P     = c(50,100,150)  #number of equations
KK     = c(50,100,150) #number of covariates 
b = 10 #nonzero entries
d = 5
Bn = c(4,6,8,10)  #length of block
Rho = c(0.1,1)
df = 8 #4

sim_scene = function(nboot,T,p,K,rho,Bn,df,d,b){ 
  beta.matrix = matrix(0, p, K)
  lambda.ga = 2*1.1*qnorm(1-0.1/(2*K*p))*sqrt(T)
  Smatrix.boot  = array(0,dim=c(p,K,nboot))
  load = matrix(0,p,K)
  c = rep(0, nboot)
  for (i in 1:p){
    indx = ceiling(i/d)
    beta.matrix[i, (d*(indx-1)+1):(indx*d)] = b
  }
  Y.data = matrix(0, T, p)
  
  innov = numeric(0)
  for(i in 1:K){
    innov = cbind(innov, rt(T+1000+1, df)/sqrt(df/(df-2)))
  }
  eta = matrix(0,T+1000,K)
  for(t in (1+1):(T+1+1000)){
    eta[t-1,] = innov[t,]*sqrt(0.8*innov[t-1,]^2+0.2)
  }
  M = A = list()
  for(j in 1:(1000+1)){
    M[[j]] = matrix(rnorm(K*K), nrow=K, ncol=K)
    A[[j]] = (j-1+1)^(-rho-1)*M[[j]]
  }
  X = matrix(0,nrow=T,ncol=K)
  for(t in 1:(T)){
    for(tt in 1:(1000+1)){
      X[t,] = X[t,] + A[[tt]]%*%eta[t-tt+1000+1,]
    }
  }
  
  innov2 = numeric(0)
  for(i in 1:p){
    innov2 = cbind(innov2, rt(T+1+1000, df)/sqrt(df/(df-2)))
  }
  eta2 = matrix(0,T+1000,p)
  for(t in (1+1):(T+1+1000)){
    eta2[t-1,] = innov2[t,]*sqrt(0.8*innov2[t-1,]^2+0.2)
  }
  M2 = A2 = list()
  for(j in 1:(1000+1)){
    M2[[j]] = matrix(rnorm(p*p), nrow=p, ncol=p)
    A2[[j]] = (j-1+1)^(-rho-1)*M2[[j]]
  }
  eps = matrix(0,nrow=T,ncol=p)
  for(t in 1:(T)){
    for(tt in 1:(1000+1)){
      eps[t,] = eps[t,] + A2[[tt]]%*%eta2[t-tt+1000+1,]
    }
  }
  
  rmse =  rmse.joint = l2 = l2.joint = array(0,dim=c(p,length(Bn),2))
  ratio.median = ratio.median2 = matrix(0,nrow=length(Bn),ncol=2)
  for(bb in 1:length(Bn)){
    bn = Bn[bb]
    ln = floor(T/bn) #number of blocks
    for (i in 1:p){
      Y = X%*%beta.matrix[i,] + eps[,i] #rnorm(T, 0, 1) #sd=1 or 0.5
      Y.data[,i] = Y
      fit1 =  rlasso(X, Y, penalty = list(homoscedastic = "none", lambda.start = 2*0.5*sqrt(T)*qnorm(1-0.1/(2*K))), post=FALSE)
      beta = fit1$beta
      intercept = fit1$intercept
      res = Y - X %*% beta - intercept * rep(1, length(Y))
      for (j in 1:nboot){
        res.boot = rnorm(ln) 
        for (k in 1:K){
          Smatrix.boot[i,k,j] = sum(c(rep(res.boot, each=bn),rep(rnorm(1),T-ln*bn))*X[,k]*res)/sqrt(T)
        }
      }
      #### compute the penalty loadings
      for(k in 1:K){
        load[i,k] = sqrt(lrvar(X[,k]*res)*T)
      }
    }
    
    for(j in 1:nboot){
      c[j]=max(abs(Smatrix.boot[,,j]/load))
    }
    lambda.boot = 2*quantile(c, 0.9)*sqrt(T)*1.1
    
    # compare the in-sample fitting performance
    for (i in 1:p){
      Y = Y.data[,i]
      fit_hdm2 =  rlasso(X, Y, penalty = list(X.dependent.lambda = FALSE))
      rmse[i,bb,] = sqrt(mean((X%*%beta.matrix[i,]-X%*%fit_hdm2$beta)^2))
      l2[i,bb,] = sqrt(sum((fit_hdm2$beta-beta.matrix[i,])^2))
      fit_hdm.joint =  rlasso(X,Y,penalty = list(homoscedastic = "none", lambda.start = lambda.boot))
      fit_hdm.joint2 =  rlasso(X,Y,penalty = list(homoscedastic = "none", lambda.start = lambda.ga))
      rmse.joint[i,bb,1] = sqrt(mean((X%*%beta.matrix[i,]-X%*%fit_hdm.joint$beta)^2))
      l2.joint[i,bb,1] = sqrt(sum((fit_hdm.joint$beta-beta.matrix[i,])^2))
      rmse.joint[i,bb,2] = sqrt(mean((X%*%beta.matrix[i,]-X%*%fit_hdm.joint2$beta)^2))
      l2.joint[i,bb,2] = sqrt(sum((fit_hdm.joint2$beta-beta.matrix[i,])^2))
    }
    
    ratio.median[bb,] = colMedians(rmse.joint[,bb,])/colMedians(rmse[,bb,])
    ratio.median2[bb,] = colMedians(l2.joint[,bb,])/colMedians(l2[,bb,])
  }

  list(ratio.median = ratio.median, ratio.median2 = ratio.median2)#}
}

for(rr in 1:length(Rho)){
  rho = Rho[rr]
  for(kk in 1:length(P)){
    p = P[kk]
    K = KK[kk]
    results = foreach(l=1:rep, .packages=c("mvtnorm", "Matrix", "hdm", "sandwich", "matrixStats"), .inorder=FALSE) %dopar%{ 
      sim_scene(nboot=nboot,T=T,p=p,K=K,rho=rho,Bn=Bn,df=df,d=d,b=b)
    }
    save(results, file = paste("ratios_rho",rho,"_K",K,".dat", sep = ""))
  }
}

ratios.mean = ratios.median = ratios.mean2 = ratios.median2 = ratios.sd = ratios.sd2 = list()
for(rr in 1:length(Rho)){
  ratios.mean[[rr]] = ratios.median[[rr]] = ratios.mean2[[rr]] =  ratios.median2[[rr]] = ratios.sd[[rr]] = ratios.sd2[[rr]] = array(0,dim=c(length(Bn),length(P),2))
}
for(rr in 1:length(Rho)){
  rho = Rho[rr]
  for(kk in 1:length(P)){
    K = KK[kk]
    load(file = paste("ratios_rho",rho,"_K",K,".dat", sep = ""))
    ratio.median = ratio.median2 = array(0,dim=c(rep,length(Bn),2))
    for(r in 1:rep){
      ratio.median[r,,] = results[[r]]$ratio.median
      ratio.median2[r,,] = results[[r]]$ratio.median2
    }
    for(bb in 1:length(Bn)){
      ratios.mean[[rr]][bb,kk,] = colMeans(ratio.median[,bb,])
      ratios.median[[rr]][bb,kk,] = colMedians(ratio.median[,bb,])
      ratios.mean2[[rr]][bb,kk,] = colMeans(ratio.median2[,bb,])
      ratios.median2[[rr]][bb,kk,] = colMedians(ratio.median2[,bb,])
      ratios.sd[[rr]][bb,kk,] = colSds(ratio.median[,bb,])
      ratios.sd2[[rr]][bb,kk,] = colSds(ratio.median2[,bb,])
    }
  }   
}
      
round(ratios.mean[[1]][,,1], digits = 4)
round(ratios.median[[1]][,,1], digits = 4)
round(ratios.sd[[1]][,,1], digits = 4)
round(ratios.mean2[[1]][,,1], digits = 4)
round(ratios.median2[[1]][,,1], digits = 4)
round(ratios.sd2[[1]][,,1], digits = 4)

round(ratios.mean[[2]][,,1], digits = 4)
round(ratios.median[[2]][,,1], digits = 4)
round(ratios.sd[[2]][,,1], digits = 4)
round(ratios.mean2[[2]][,,1], digits = 4)
round(ratios.median2[[2]][,,1], digits = 4)
round(ratios.sd2[[2]][,,1], digits = 4)

```

automatically created on 2019-04-25