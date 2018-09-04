[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **LTS_inference** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

﻿
Name of Quantlet: LTS_inference

Published in: LASSO_Time_Space

Description: 'Compares the average rejection rate of the null hypotheses for the individual and simultaneous inference (confidence intervals constructed by asymptotics or bootstrap).'

Keywords : 
- simulation
- lasso shrinkage
- bootstrap
- regression
- hypothesis-testing

Author:          Chen Huang

Submitted:        Tuesday, June 19 2018 by Chen Huang

```

### R Code
```r

libraries = c("mvtnorm", "hdm", "pracma", "np", "sandwich", "flare", "doSNOW", "AER", "evd")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# single target variable in each equation
# med.iv = function(alpha, y, x, z){
#   u = y - alpha*x
#   phi = rep(0,length(u))
#   for(t in 1:length(u)){
#     if(u[t]<=0){phi[t]=-0.5}else{phi[t]=0.5}
#   }
#   4*mean(phi*z)^2/mean(z^2)
# }

rep = 100
nboot = 100
T = 100
K = 100
M = 50
bn = 25 #length of block
rho = 0.5
rho.var = 0.1
lag = 1
dx = 5
dd = 5
cd = 0.25
cy = 0.5
alpha0 = seq(0,1.5,length.out = 10)

c1 = makeCluster(6) # core numbers: cohen 20 cores
registerDoSNOW(c1)

#results = list()
sim_scene = function(nboot,T,K,M,rho,rho.var,lag,dx,dd,cd,cy,bn,alpha0){ 
  #for(r in 1:rep){
  lambda.ga = 2*1.1*qnorm(1-0.1/(2*K*M))
  dn = sqrt(2*log(K*M)) - (log(log(K*M))+log(4*pi))/(2*sqrt(2*log(K*M)))
  #alpha = c(rep(1,M/2),rep(0,M/2))
  alpha = rep(alpha0,M)
  ln = T/bn #number of blocks
  beta.matrix = matrix(0, M, K)
  for (m in 1:M){
    indx = ceiling(m/dd)
    beta.matrix[m, (dd*(indx-1)+1):(indx*dd)] = cy/c(1:dd)
  }
  theta.matrix = matrix(0, M, K)
  for (m in 1:M){
    indx = ceiling(m/dd)
    theta.matrix[m, (dd*(indx-1)+1):(indx*dd)] = cd/c(1:dd)
  }
  Phi = bdiag(matrix(rho.var,dx,dx))
  for(pp in 1:(K/dx-1)){
    Phi = bdiag(Phi,matrix(rho.var,dx,dx))
  }
  Phi = as.matrix(Phi)
  Y = matrix(0,T,M)
  d = matrix(0,T,M)
  v = matrix(0,T,M)
  psi = matrix(0,T,M)
  load = matrix(0,M,K)
  load1 = matrix(0,M,K+1)
  beta.matrix.hat = matrix(0,M,K)
  theta.matrix.hat = matrix(0,M,K)
  intercept.hat = rep(0,M)
  e = matrix(0, nrow = T, ncol = M)
  c = rep(0, nboot)
  Smatrix.boot  = array(0,dim=c(M,K,nboot))
  Smatrix.boot1  = array(0,dim=c(M,K+1,nboot))
  omega = rep(0, M)
  phi = rep(0, M)
  #upper_0.05 = rep(0, M)
  #lower_0.05 = rep(0, M)
  #upper.boot_0.05 = rep(0, M)
  #lower.boot_0.05 = rep(0, M)
  stat.boot = matrix(0, nboot, M)
  
  alpha.hat = rep(0,M)
  alpha0.hat = rep(0,M)
  sigma = rep(0,M)
  
  #length_0.05 = rep(0,M)
  #length.boot_0.05 = rep(0,M)
  
  reject_0.01 = rep(0,M)
  reject_0.05 = rep(0,M)
  reject_0.1 = rep(0,M)
  reject.simul_0.01 = 0
  reject.simul_0.05 = 0
  reject.simul_0.1 = 0
  reject.boot_0.01 = rep(0,M)
  reject.boot_0.05 = rep(0,M)
  reject.boot_0.1 = rep(0,M)
  #reject.boot.simul_0.05_power = 0
  #reject.boot.simul_0.05_size = 0
  reject.boot.simul_0.01 = 0
  reject.boot.simul_0.05 = 0
  reject.boot.simul_0.1 = 0
  
  ## errors~AR(1)
  innov1  = rmvnorm(T, rep(0,M), diag(M)) 
  innov2  = rmvnorm(T, rep(0,M), diag(M))
  for (t in 2:T) {
    innov1[t,] = rho * innov1[t-1,] + rnorm(1)
    innov2[t,] = rho * innov2[t-1,] + rnorm(1)
  }
  ## X~VAR(lag)
  X = matrix(0,nrow=T+lag,ncol=K)
  innov = rmvnorm((T+lag), rep(0,K), diag(K))
  X[1,] = innov[1,]
  for(t in 2:(T+lag)){
    X[t,] = Phi %*% X[t-lag,] + innov[t,]
  }
  X = X[-1,]
  
  ## joint penalty level (for step 1)
  for (m in 1:M){
    d[,m] = X %*% theta.matrix[m,] + innov1[,m]
    #d[,m] = scale(d[,m])
    Y[,m] = alpha[m]*d[,m] + X %*% beta.matrix[m,] + innov2[,m]
    #Y[,m] = scale(Y[,m])
    #X = scale(X)
    X.all = cbind(d[,m], X)
    fit0 = rlasso(X.all, Y[,m], penalty = list(homoscedastic = "none", lambda.start = 2*0.5*sqrt(T)*qnorm(1-0.1/(2*(K+1)))), post=FALSE)
    #load1[m,] = fit0$loadings
    e[,m] = Y[,m] - predict(fit0)
    for (j in 1:nboot){
      e.boot = rnorm(ln)
      for (k in 1:(K+1)){
        sum = 0
        for (l in 1:ln){
          sum = sum + sum(X.all[((l-1)*bn+1):(l*bn),k]*e[((l-1)*bn+1):(l*bn),m])*e.boot[l]
        }
        Smatrix.boot1[m,k,j] = sum/sqrt(T)
      }
    }
    #### compute the penalty loadings
    for(k in 1:(K+1)){
      load1[m,k] = sqrt(lrvar(X.all[,k]*e[,m])*T)
    }
    
  }
  for(j in 1:nboot){
    c[j]=max(abs(Smatrix.boot1[,,j]/load1))
  }
  lambda.boot1 = 2*quantile(c, 0.9)*sqrt(T)*1.1
  
  ## joint penalty level (for step 2)
  for (m in 1:M){
    fit0 =  rlasso(X, d[,m], penalty = list(homoscedastic = "none", lambda.start = 2*0.5*sqrt(T)*qnorm(1-0.1/(2*K))), post=FALSE)
    #load[m,] = fit0$loadings
    e[,m] = d[,m] - predict(fit0)
    for (j in 1:nboot){
      e.boot = rnorm(ln)
      for (k in 1:K){
        sum = 0
        for (l in 1:ln){
          sum = sum + sum(X[((l-1)*bn+1):(l*bn),k]*e[((l-1)*bn+1):(l*bn),m])*e.boot[l]
        }
        Smatrix.boot[m,k,j] = sum/sqrt(T)
      }
    }
    #### compute the penalty loadings
    for(k in 1:K){
      load[m,k] = sqrt(lrvar(X[,k]*e[,m])*T)
    }
    
  }
  for(j in 1:nboot){
    c[j]=max(abs(Smatrix.boot[,,j]/load))
  }
  lambda.boot = 2*quantile(c, 0.9)*sqrt(T)*1.1
  
  #reject_0.05 = rep(0,M)
  for (m in 1:M){
    # d[,m] = X %*% theta.matrix[m,] + innov1[,m]
    #d[,m] = scale(d[,m])
    # Y[,m] = alpha[m]*d[,m] + X %*% beta.matrix[m,] + innov2[,m]
    #Y[,m] = scale(Y[,m])
    #X = scale(X)
    
    X.all = cbind(d[,m], X)
    
    ## step 1
    fit1 =  rlasso(X.all, Y[,m], penalty = list(homoscedastic = "none", lambda.start = lambda.boot1))
    #fit1 =  rlasso(X.all, Y[,m])
    #fit1 =  rlasso(X.all, Y[,m], penalty = list(homoscedastic = "none", lambda.start = lambda.ga1))
    
    beta.matrix.hat[m,] = fit1$beta[-1]
    intercept.hat[m] = fit1$intercept
    
    ## step 2
    fit2 =  rlasso(X, d[,m], penalty = list(homoscedastic = "none", lambda.start = lambda.boot))
    #fit2 =  rlasso(X, d[,m])
    #fit2 =  rlasso(X, d[,m], penalty = list(homoscedastic = "none", lambda.start = lambda.ga))
    theta.matrix.hat[m,] = fit2$beta
    v[,m] = d[,m] - predict(fit2)
    #v[,m] = scale(v[,m])
    
    y = Y[,m] - X %*% beta.matrix.hat[m,] - rep(intercept.hat[m], dim(Y)[1])
    
    ### LS-IV
    fit.iv = summary(ivreg(y ~ 0 + d[,m] | v[,m]))
    alpha.hat[m] = fit.iv$coefficients[,"Estimate"]
    psi[,m] = v[,m]*(y - alpha.hat[m]*d[,m])
    phi[m] = -mean(v[,m]^2)
    ### LS-IV
    
    ### LAD-IV
    # alpha0.hat[m] = fit1$beta[1]
    # b = sqrt(mean(d[,m]^2))*log(T)
    # alpha.hat[m] = optimize(med.iv, c(alpha0.hat[m]-10/b, alpha0.hat[m]+10/b), y = y, x = d[,m], z = v[,m])$minimum
    # #alpha.hat[m] = gridSearch(med.iv, levels = list(seq(alpha0.hat[m]-10/b, alpha0.hat[m]+10/b, length.out = 10000)), y = y, x = d[,m], z = v[,m], printDetail = FALSE)$minlevels
    # 
    # for(t in 1:T){
    #   if((y[t] - alpha.hat[m]*d[t,m])<=0){psi[t,m]=-0.5*v[t,m]}else{psi[t,m]=0.5*v[t,m]}
    # }

    # error = y - alpha.hat[m]*d[,m]
    # ind    = which(density(error)$x == min(abs(density(error)$x)))
    # # if (length(ind) == 0) {
    # #   f = density(error)$y[which(density(error)$x == -min(abs(density(error)$x)))]
    # # } else {
    # #   f = density(error)$y[which(density(error)$x == min(abs(density(error)$x)))]
    # # }
    # f = npudens(tdat = error, edat = 0)$dens
    # phi[m] = -mean(v[,m]^2*f)
    ### LAD-IV
    
    #omega[m] = var(psi[,m])
    omega[m] = lrvar(psi[,m])*T

    sigma[m] = sqrt(omega[m]/phi[m]^2)
    
    #### asymptotic individual inference
    #lower_0.05[m] = alpha.hat[m]-sigma[m]*qnorm(0.975)/sqrt(T)
    #upper_0.05[m] = alpha.hat[m]+sigma[m]*qnorm(0.975)/sqrt(T)
    if(abs(alpha.hat[m]*sqrt(T)/sigma[m])>qnorm(0.995)){reject_0.01[m] = 1} 
    if(abs(alpha.hat[m]*sqrt(T)/sigma[m])>qnorm(0.975)){reject_0.05[m] = 1} 
    if(abs(alpha.hat[m]*sqrt(T)/sigma[m])>qnorm(0.95)){reject_0.1[m] = 1} 
    #### asymptotic individual inference
  }
  #length_0.05 = upper_0.05 - lower_0.05
  
  #### bootstrap individual inference  
  for(m in 1:M){
    inf = -psi[,m]/(phi[m]*sigma[m])
    for(j in 1:nboot){
      e.boot = rnorm(ln)
      sum = 0
      for (l in 1:ln){
        sum = sum + sum(inf[((l-1)*bn+1):(l*bn)])*e.boot[l]
      }
      stat.boot[j,m] = sum/sqrt(T)
    }
    
    
    #upper.boot_0.05[m] = alpha.hat[m] - quantile(stat.boot[,m],0.025)*sigma[m]/sqrt(T)
    #lower.boot_0.05[m] = alpha.hat[m] - quantile(stat.boot[,m],0.975)*sigma[m]/sqrt(T)
    if((sqrt(T)*alpha.hat[m]/sigma[m])<=quantile(stat.boot[,m],0.005) || (sqrt(T)*alpha.hat[m]/sigma[m])>=quantile(stat.boot[,m],0.995)){reject.boot_0.01[m] = 1}
    if((sqrt(T)*alpha.hat[m]/sigma[m])<=quantile(stat.boot[,m],0.025) || (sqrt(T)*alpha.hat[m]/sigma[m])>=quantile(stat.boot[,m],0.975)){reject.boot_0.05[m] = 1}
    if((sqrt(T)*alpha.hat[m]/sigma[m])<=quantile(stat.boot[,m],0.05) || (sqrt(T)*alpha.hat[m]/sigma[m])>=quantile(stat.boot[,m],0.95)){reject.boot_0.1[m] = 1}
    
  }
  #length.boot_0.05 = upper.boot_0.05 - lower.boot_0.05
  
  #### bootstrap individual inference
  
  #### bootstrap simultaneous inference
  stat.max.boot = 0
  for(j in 1:nboot){
    stat.max.boot = c(stat.max.boot, max(abs(stat.boot[j,1:(M)])))
  }
  stat.max.boot = stat.max.boot[-1]
  for(m in 1:(M)){
    if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.max.boot,0.99) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.99)){reject.boot.simul_0.01 = 1}
    if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.max.boot,0.95) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.95)){reject.boot.simul_0.05 = 1}
    if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.max.boot,0.90) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.90)){reject.boot.simul_0.1 = 1}
  }
  #### bootstrap simultaneous inference 
  
  #### asymptotic simultaneous inference
  for(m in 1:(M)){
    #if(sqrt(T)*alpha.hat[m]/sigma[m]<=qgumbel(0.005)/(2*log(K*M))+dn || sqrt(T)*alpha.hat[m]/sigma[m]>=qgumbel(0.995)/(2*log(K*M))+dn){reject.simul_0.01 = 1}
    #if(sqrt(T)*alpha.hat[m]/sigma[m]<=qgumbel(0.025)/(2*log(K*M))+dn || sqrt(T)*alpha.hat[m]/sigma[m]>=qgumbel(0.975)/(2*log(K*M))+dn){reject.simul_0.05 = 1}
    #if(sqrt(T)*alpha.hat[m]/sigma[m]<=qgumbel(0.05)/(2*log(K*M))+dn || sqrt(T)*alpha.hat[m]/sigma[m]>=qgumbel(0.95)/(2*log(K*M))+dn){reject.simul_0.1 = 1}
    if(sqrt(T)*alpha.hat[m]/sigma[m]<=qgumbel(0.005) || sqrt(T)*alpha.hat[m]/sigma[m]>=qgumbel(0.995)){reject.simul_0.01 = 1}
    if(sqrt(T)*alpha.hat[m]/sigma[m]<=qgumbel(0.025) || sqrt(T)*alpha.hat[m]/sigma[m]>=qgumbel(0.975)){reject.simul_0.05 = 1}
    if(sqrt(T)*alpha.hat[m]/sigma[m]<=qgumbel(0.05) || sqrt(T)*alpha.hat[m]/sigma[m]>=qgumbel(0.95)){reject.simul_0.1 = 1}
  }
  #### asymptotic simultaneous inference
  
  # #### bootstrap simultaneous inference (power)
  # stat.max.boot = 0
  # for(j in 1:nboot){
  #   stat.max.boot = c(stat.max.boot, max(abs(stat.boot[j,1:(M/2)])))
  # }
  # stat.max.boot = stat.max.boot[-1]
  # for(m in 1:(M/2)){
  #   if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.max.boot,0.95) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.95)){reject.boot.simul_0.05_power = 1}
  # }
  # #### bootstrap simultaneous inference (power)
  
  # #### bootstrap simultaneous inference (size)
  # stat.max.boot = 0
  # for(j in 1:nboot){
  #   stat.max.boot = c(stat.max.boot, max(abs(stat.boot[j,(M/2+1):M])))
  # }
  # stat.max.boot = stat.max.boot[-1]
  # for(m in (M/2+1):M){
  #   if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.max.boot,0.95) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.95)){reject.boot.simul_0.05_size = 1}
  # }
  # #### bootstrap simultaneous inference (size)
  
  
  list(#length_0.05 = length_0.05, length.boot_0.05 = length.boot_0.05,
    reject_0.05 = reject_0.05, reject.boot_0.05 = reject.boot_0.05,
    reject.simul_0.05 = reject.simul_0.05, reject.boot.simul_0.05 = reject.boot.simul_0.05,
    reject_0.01 = reject_0.01, reject.boot_0.01 = reject.boot_0.01,
    reject.simul_0.01 = reject.simul_0.01, reject.boot.simul_0.01 = reject.boot.simul_0.01,
    reject_0.1 = reject_0.1, reject.boot_0.1 = reject.boot_0.1,
    reject.simul_0.1 = reject.simul_0.1, reject.boot.simul_0.1 = reject.boot.simul_0.1)
  #reject.boot.simul_0.05_power = reject.boot.simul_0.05_power, reject.boot.simul_0.05_size = reject.boot.simul_0.05_size)
}

results = foreach(l=1:rep, .packages=c("mvtnorm", "hdm", "pracma", "np", "sandwich", "flare", "AER", "evd"), .inorder=FALSE) %dopar%{ 
  sim_scene(nboot=nboot,T=T,K=K,M=M,rho=rho,rho.var=rho.var,lag=lag,dx=dx,dd=dd,cd=cd,cy=cy,bn=bn,alpha0=alpha0)
}

save(results, file = "para.results_K100_M50_bn25.LS_a1.dat")
load(file = "para.results_K100_M50_bn25.LS_a1.dat")

for(a in 2:length(alpha0)){
  results = foreach(l=1:rep, .packages=c("mvtnorm", "hdm", "pracma", "np", "sandwich", "flare", "AER", "evd"), .inorder=FALSE) %dopar%{ 
    sim_scene(nboot=nboot,T=T,K=K,M=M,rho=rho,rho.var=rho.var,lag=lag,dx=dx,dd=dd,cd=cd,cy=cy,bn=bn,alpha0=alpha0[a])
  }
  save(results, file = paste("para.results_K100_M50_bn25.LS_a",a,".dat",sep=""))
}

par(mfrow = c(2, 2), mgp = c(1.5, 0.5, 0), oma = c(1, 1, 1, 1), mar = c(2.5, 1.8, 1, 1), xpd = NA, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

rep = 100
M = 200
alpha0 = seq(0,1.5,length.out = 10)

rej.asy = matrix(0,10,3)
rej.asy.simul = matrix(0,10,3)
rej.boot = matrix(0,10,3)
rej.boot.simul = matrix(0,10,3)

for(a in 1:10){
  load(file = paste("para.results_K200_M200_bn25.LS_a",a,".dat",sep=""))
  
  reject_0.01 = matrix(0,rep,M)
  reject_0.05 = matrix(0,rep,M)
  reject_0.1 = matrix(0,rep,M)
  reject.boot_0.01 = matrix(0,rep,M)
  reject.boot_0.05 = matrix(0,rep,M)
  reject.boot_0.1 = matrix(0,rep,M)
  
  for(r in 1:rep){
    reject_0.01[r,] = results[[r]]$reject_0.01
    reject_0.05[r,] = results[[r]]$reject_0.05
    reject_0.1[r,] = results[[r]]$reject_0.1
    reject.boot_0.01[r,] = results[[r]]$reject.boot_0.01
    reject.boot_0.05[r,] = results[[r]]$reject.boot_0.05
    reject.boot_0.1[r,] = results[[r]]$reject.boot_0.1
  }
  
  rej.asy[a,1] = mean(colMeans(reject_0.01)[1:(M)])
  rej.asy[a,2] = mean(colMeans(reject_0.05)[1:(M)]) 
  rej.asy[a,3] = mean(colMeans(reject_0.1)[1:(M)]) 
  
  rej.boot[a,1] = mean(colMeans(reject.boot_0.01)[1:(M)]) 
  rej.boot[a,2] = mean(colMeans(reject.boot_0.05)[1:(M)]) 
  rej.boot[a,3] = mean(colMeans(reject.boot_0.1)[1:(M)]) 
  
  reject.boot.simul_0.01 = rep(0,rep)
  for(r in 1:rep){
    reject.boot.simul_0.01[r] = results[[r]]$reject.boot.simul_0.01
  }
  rej.boot.simul[a,1] = mean(reject.boot.simul_0.01)  
  
  reject.boot.simul_0.05 = rep(0,rep)
  for(r in 1:rep){
    reject.boot.simul_0.05[r] = results[[r]]$reject.boot.simul_0.05
  }
  rej.boot.simul[a,2] = mean(reject.boot.simul_0.05) 
  
  reject.boot.simul_0.1 = rep(0,rep)
  for(r in 1:rep){
    reject.boot.simul_0.1[r] = results[[r]]$reject.boot.simul_0.1
  }
  rej.boot.simul[a,3] = mean(reject.boot.simul_0.1)
  
  reject.simul_0.01 = rep(0,rep)
  for(r in 1:rep){
    reject.simul_0.01[r] = results[[r]]$reject.simul_0.01
  }
  rej.asy.simul[a,1] = mean(reject.simul_0.01)
  
  reject.simul_0.05 = rep(0,rep)
  for(r in 1:rep){
    reject.simul_0.05[r] = results[[r]]$reject.simul_0.05
  }
  rej.asy.simul[a,2] = mean(reject.simul_0.05)
  
  reject.simul_0.1 = rep(0,rep)
  for(r in 1:rep){
    reject.simul_0.1[r] = results[[r]]$reject.simul_0.1
  }
  rej.asy.simul[a,3] = mean(reject.simul_0.1)
}

plot(alpha0, rej.asy[,2], xlab = expression(alpha^0), ylab = "Rejection rate", 
     type = "l", ylim = c(0, 1), lwd = 3)
lines(alpha0, rej.boot[, 2], lwd = 3, lty = 2)
lines(alpha0, rej.asy.simul[, 2], lwd = 3, lty = 3)
#lines(alpha0, rej.boot.simul[, 2], lwd = 3, lty = 4)
#lines(alpha0, rep(0.05,10), lwd = 3, lty = 5)
title(main = "K=200, J=200", line = 0.2, lwd = 1)
#legend(x=0.9, y=0.4, legend=c("Ind. Asym.","Ind. Boot","Simult. Asym.","Simult. Boot"), lty=c(1,2,3,4), border = "transparent", cex = 1, box.col = "transparent", bg="transparent")
legend(x=0.9, y=0.4, legend=c("Ind. Asym.","Ind. Boot","Simult. Boot"), lty=c(1,2,3), border = "transparent", cex = 1, box.col = "transparent", bg="transparent")


```

automatically created on 2018-09-04