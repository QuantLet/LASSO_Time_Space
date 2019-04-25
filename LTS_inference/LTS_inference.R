libraries = c("mvtnorm", "hdm", "pracma", "np", "sandwich", "flare", "doSNOW", "AER", "evd")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

med.iv = function(alpha, y, x, z){
  u = y - alpha*x
  phi = rep(0,length(u))
  for(t in 1:length(u)){
    if(u[t]<=0){phi[t]=-0.5}else{phi[t]=0.5}
  }
  4*mean(phi*z)^2/mean(z^2)
}

rep = 1000
nboot = 5000
T = 100
KK = c(50,100,150)
MM = c(50,100,150)
Bn = c(4,6,8,10)  #length of block
Rho = c(0.1,1)
dx = 5
dd = 5
df = 8 #4
cd = 0.25
cy = 5
alpha0 = c(0, 2.5, 5)

c1 = makeCluster(20) # core numbers
registerDoSNOW(c1)

sim_scene = function(nboot,T,K,M,rho,df,dx,dd,cd,cy,Bn,a0){ 
  alpha = runif(M,0,a0)  

  #lambda.ga = 2*1.1*qnorm(1-0.1/(2*K*M))*sqrt(T)
  beta.matrix = matrix(0, M, K)
  for (m in 1:M){
    indx = ceiling(m/dd)
    beta.matrix[m, (dd*(indx-1)+1):(indx*dd)] = runif(dd,0,cy)#cy#/c(1:dd)
  }
  theta.matrix = matrix(0, M, K)
  for (m in 1:M){
    indx = ceiling(m/dd)
    theta.matrix[m, (dd*(indx-1)+1):(indx*dd)] = runif(dd,0,cd)#cd/c(1:dd)
  }
  
  Y = matrix(0,T,M)
  d = matrix(0,T,M)
  v = matrix(0,T,M)
  
  load = matrix(0,M,K)
  load1 = matrix(0,M,K+1)
  e = matrix(0, nrow = T, ncol = M)
  c = rep(0, nboot)
  Smatrix.boot  = array(0,dim=c(M,K,nboot))
  Smatrix.boot1  = array(0,dim=c(M,K+1,nboot))
  
  alpha.hat = rep(0,M)
  sigma = rep(0,M)
  phi = rep(0, M)
  psi = matrix(0,T,M)
  stat = rep(0,M)
  pv = rep(0,M)
  stat.boot = matrix(0, nboot, M)
  
  p = M
  innov = numeric(0)
  for(i in 1:K){
    innov = cbind(innov, rt(T+1+1000, df)/sqrt(df/(df-2)))
  }
  eta = matrix(0,T+1000,K)
  for(t in (1+1):(T+1+1000)){
    eta[t-1,] = innov[t,]*sqrt(0.8*innov[t-1,]^2+0.2)
  }
  B = A = list()
  for(j in 1:(1000+1)){
    B[[j]] = matrix(rnorm(K*K), nrow=K, ncol=K)
    A[[j]] = (j-1+1)^(-rho-1)*B[[j]]
  }
  X = matrix(0,nrow=T,ncol=K)
  for(t in 1:T){
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
  B2 = A2 = list()
  for(j in 1:(1000+1)){
    B2[[j]] = matrix(rnorm(p*p), nrow=p, ncol=p)
    A2[[j]] = (j-1+1)^(-rho-1)*B2[[j]]
  }
  eps1 = matrix(0,nrow=T,ncol=p)
  for(t in 1:T){
    for(tt in 1:(1000+1)){
      eps1[t,] = eps1[t,] + A2[[tt]]%*%eta2[t-tt+1000+1,]
    }
  }
  
  innov3 = numeric(0)
  for(i in 1:p){
    innov3 = cbind(innov3, rt(T+1+1000, df)/sqrt(df/(df-2)))
  }
  eta3 = matrix(0,T+1000,p)
  for(t in (1+1):(T+1+1000)){
    eta3[t-1,] = innov3[t,]*sqrt(0.8*innov3[t-1,]^2+0.2)
  }
  B3 = A3 = list()
  for(j in 1:(1000+1)){
    B3[[j]] = matrix(rnorm(p*p), nrow=p, ncol=p)
    A3[[j]] = (j-1+1)^(-rho-1)*B3[[j]]
  }
  eps2 = matrix(0,nrow=T,ncol=p)
  for(t in 1:T){
    for(tt in 1:(1000+1)){
      eps2[t,] = eps2[t,] + A3[[tt]]%*%eta3[t-tt+1000+1,]
    }
  }
  
  reject = reject.boot = reject.boot.multi = array(0,dim=c(length(Bn),3,M))
  reject.boot.simul = matrix(0,nrow=length(Bn),ncol=3)
  
  for(bb in 1:length(Bn)){
    bn = Bn[bb]
    ln = floor(T/bn) #number of blocks
    ## joint penalty level (for step 1)
    for (m in 1:M){
      d[,m] = X %*% theta.matrix[m,] + eps1[,m]
      Y[,m] = alpha[m]*d[,m] + X %*% beta.matrix[m,] + eps2[,m]
      X.all = cbind(d[,m], X)
      fit0 = rlasso(X.all, Y[,m], penalty = list(homoscedastic = "none", lambda.start = 2*0.5*sqrt(T)*qnorm(1-0.1/(2*(K+1)))), post=FALSE)
      e[,m] = Y[,m] - predict(fit0)
      for (j in 1:nboot){
        e.boot = rnorm(ln)
        for (k in 1:(K+1)){
          Smatrix.boot1[m,k,j] = sum(c(rep(e.boot, each=bn),rep(rnorm(1),T-ln*bn))*X.all[,k]*e[,m])/sqrt(T)
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
      e[,m] = d[,m] - predict(fit0)
      for (j in 1:nboot){
        e.boot = rnorm(ln)
        for (k in 1:K){
          Smatrix.boot[m,k,j] = sum(c(rep(e.boot, each=bn),rep(rnorm(1),T-ln*bn))*X[,k]*e[,m])/sqrt(T)
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
    
    reject_0.01 = rep(0,M)
    reject_0.05 = rep(0,M)
    reject_0.1 = rep(0,M)
    reject.boot_0.01 = rep(0,M)
    reject.boot_0.05 = rep(0,M)
    reject.boot_0.1 = rep(0,M)
    reject.boot.simul_0.01 = 0
    reject.boot.simul_0.05 = 0
    reject.boot.simul_0.1 = 0
    reject.boot.multi_0.01 = rep(0,M)
    reject.boot.multi_0.05 = rep(0,M)
    reject.boot.multi_0.1 = rep(0,M)
    
    for (m in 1:M){
      fit1 =  rlasso(cbind(d[,m], X), Y[,m], penalty = list(homoscedastic = "none", lambda.start = lambda.boot1))
      fit2 =  rlasso(X, d[,m], penalty = list(homoscedastic = "none", lambda.start = lambda.boot))
      v[,m] = d[,m] - predict(fit2)
      y = Y[,m] - X %*% fit1$beta[-1] - rep(fit1$intercept, dim(Y)[1])
      b = sqrt(mean(d[,m]^2))*log(T)
      alpha.hat[m] = optimize(med.iv, c(fit1$beta[1]-10/b, fit1$beta[1]+10/b), y = y, x = d[,m], z = v[,m])$minimum
      for(t in 1:T){
        if((y[t] - alpha.hat[m]*d[t,m])<=0){psi[t,m]=-0.5*v[t,m]}else{psi[t,m]=0.5*v[t,m]}
      }
      error = y - alpha.hat[m]*d[,m]
      f = npudens(tdat = error, edat = 0)$dens
      phi[m] = -mean(v[,m]^2*f)
      sigma[m] = sqrt(lrvar(psi[,m])*T/phi[m]^2)
      
      #### asymptotic individual inference
      if(abs(alpha.hat[m]*sqrt(T)/sigma[m])>qnorm(0.995)){reject_0.01[m] = 1}
      if(abs(alpha.hat[m]*sqrt(T)/sigma[m])>qnorm(0.975)){reject_0.05[m] = 1}
      if(abs(alpha.hat[m]*sqrt(T)/sigma[m])>qnorm(0.95)){reject_0.1[m] = 1}
      #### asymptotic individual inference
    }
    
    #### bootstrap individual inference  
    for(m in 1:M){
      inf = -psi[,m]/(phi[m]*sigma[m])
      for(j in 1:nboot){
        e.boot = rnorm(ln)
        stat.boot[j,m] = sum(inf*c(rep(e.boot, each=bn),rep(rnorm(1),T-ln*bn)))/sqrt(T)
      }
      
      if((sqrt(T)*alpha.hat[m]/sigma[m])<=-quantile(abs(stat.boot[,m]),0.99) || (sqrt(T)*alpha.hat[m]/sigma[m])>=quantile(abs(stat.boot[,m]),0.99)){reject.boot_0.01[m] = 1}
      if((sqrt(T)*alpha.hat[m]/sigma[m])<=-quantile(abs(stat.boot[,m]),0.95) || (sqrt(T)*alpha.hat[m]/sigma[m])>=quantile(abs(stat.boot[,m]),0.95)){reject.boot_0.05[m] = 1}
      if((sqrt(T)*alpha.hat[m]/sigma[m])<=-quantile(abs(stat.boot[,m]),0.90) || (sqrt(T)*alpha.hat[m]/sigma[m])>=quantile(abs(stat.boot[,m]),0.90)){reject.boot_0.1[m] = 1}
    }
    #### bootstrap individual inference
    
    #### bootstrap simultaneous inference
    #stat.max.boot = numeric(0)
    #stat.min.boot = numeric(0)
    stat.abs.boot = numeric(0)
    for(j in 1:nboot){
      #stat.max.boot = c(stat.max.boot, max(stat.boot[j,1:(M)]))
      #stat.min.boot = c(stat.min.boot, min(stat.boot[j,1:(M)]))
      stat.abs.boot = c(stat.abs.boot, max(abs(stat.boot[j,1:(M)])))
    }
    for(m in 1:(M)){
      if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.abs.boot,0.99) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.abs.boot,0.99)){
        reject.boot.simul_0.01 = 1
        reject.boot.multi_0.01[m] = 1
      }
      if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.abs.boot,0.95) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.abs.boot,0.95)){
        reject.boot.simul_0.05 = 1
        reject.boot.multi_0.05[m] = 1
      }
      if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.abs.boot,0.90) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.abs.boot,0.90)){
        reject.boot.simul_0.1 = 1
        reject.boot.multi_0.1[m] = 1
      }
      
      #if(sqrt(T)*alpha.hat[m]/sigma[m]<=quantile(stat.min.boot,0.005) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.995)){reject.boot.simul_0.01 = 1}
      #if(sqrt(T)*alpha.hat[m]/sigma[m]<=quantile(stat.min.boot,0.025) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.975)){reject.boot.simul_0.05 = 1}
      #if(sqrt(T)*alpha.hat[m]/sigma[m]<=quantile(stat.min.boot,0.05) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.max.boot,0.95)){reject.boot.simul_0.1 = 1}
    }
    if(length(which(reject.boot.multi_0.01!=0))>0 && length(which(reject.boot.multi_0.01!=0))<M){
      reject.boot.multi_0.01_temp = reject.boot.multi_0.01
      l=2
      while(l<=M){
        set = which(reject.boot.multi_0.01_temp==0)
        stat.abs.boot = numeric(0)
        for(j in 1:nboot){
          stat.abs.boot = c(stat.abs.boot, max(abs(stat.boot[j,set])))
        }
        for(m in set){
          if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.abs.boot,0.99) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.abs.boot,0.99)){reject.boot.multi_0.01_temp[m] = 1}
        }
        set.new = which(reject.boot.multi_0.01_temp==0)
        if(length(set.new)==0 || length(set)==length(set.new)){
          reject.boot.multi_0.01 = reject.boot.multi_0.01_temp
          l = M + 1
        }else{
          l = l + 1
        }
      }
    }
    if(length(which(reject.boot.multi_0.05!=0))>0 && length(which(reject.boot.multi_0.05!=0))<M){
      reject.boot.multi_0.05_temp = reject.boot.multi_0.05
      l=2
      while(l<=M){
        set = which(reject.boot.multi_0.05_temp==0)
        stat.abs.boot = numeric(0)
        for(j in 1:nboot){
          stat.abs.boot = c(stat.abs.boot, max(abs(stat.boot[j,set])))
        }
        for(m in set){
          if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.abs.boot,0.95) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.abs.boot,0.95)){reject.boot.multi_0.05_temp[m] = 1}
        }
        set.new = which(reject.boot.multi_0.05_temp==0)
        if(length(set.new)==0 || length(set)==length(set.new)){
          reject.boot.multi_0.05 = reject.boot.multi_0.05_temp
          l = M + 1
        }else{
          l = l + 1
        }
      }
    }
    if(length(which(reject.boot.multi_0.1!=0))>0 && length(which(reject.boot.multi_0.1!=0))<M){
      reject.boot.multi_0.1_temp = reject.boot.multi_0.1
      l=2
      while(l<=M){
        set = which(reject.boot.multi_0.1_temp==0)
        stat.abs.boot = numeric(0)
        for(j in 1:nboot){
          stat.abs.boot = c(stat.abs.boot, max(abs(stat.boot[j,set])))
        }
        for(m in set){
          if(sqrt(T)*alpha.hat[m]/sigma[m]<=-quantile(stat.abs.boot,0.90) || sqrt(T)*alpha.hat[m]/sigma[m]>=quantile(stat.abs.boot,0.90)){reject.boot.multi_0.1_temp[m] = 1}
        }
        set.new = which(reject.boot.multi_0.1_temp==0)
        if(length(set.new)==0 || length(set)==length(set.new)){
          reject.boot.multi_0.1 = reject.boot.multi_0.1_temp
          l = M + 1
        }else{
          l = l + 1
        }
      }
    }
    #### bootstrap simultaneous inference 
    reject[bb,1,] = reject_0.01
    reject[bb,2,] = reject_0.05
    reject[bb,3,] = reject_0.1
    reject.boot[bb,1,] = reject.boot_0.01
    reject.boot[bb,2,] = reject.boot_0.05
    reject.boot[bb,3,] = reject.boot_0.1
    reject.boot.multi[bb,1,] = reject.boot.multi_0.01
    reject.boot.multi[bb,2,] = reject.boot.multi_0.05
    reject.boot.multi[bb,3,] = reject.boot.multi_0.1
    reject.boot.simul[bb,1] = reject.boot.simul_0.01
    reject.boot.simul[bb,2] = reject.boot.simul_0.05
    reject.boot.simul[bb,3] = reject.boot.simul_0.1
  }
  
  list(reject = reject, reject.boot = reject.boot, reject.boot.simul = reject.boot.simul, reject.boot.multi = reject.boot.multi)
}

for(rr in 1:length(Rho)){
  rho = Rho[rr]
  for(kk in 1:length(KK)){
    K = KK[kk]
    M = MM[kk]
    for(a in 1:length(alpha0)){
      a0 = alpha0[a]
      results = foreach(l=1:rep, .packages=c("mvtnorm", "hdm", "pracma", "np", "sandwich", "flare", "AER", "evd"), .inorder=FALSE) %dopar%{ 
        sim_scene(nboot=nboot,T=T,K=K,M=M,rho=rho,df=df,dx=dx,dd=dd,cd=cd,cy=cy,Bn=Bn,a0=a0)
      }
      save(results, file = paste("test_rho",rho,"_K",K,"_a0",a0,".dat", sep = ""))
    }
  }
}

rej.asy = rej.boot = rej.boot.simul = rej.boot.multi = list()
for(rr in 1:length(Rho)){
  rej.asy[[rr]] = rej.boot[[rr]] = rej.boot.simul[[rr]] = rej.boot.multi[[rr]] = array(0,dim=c(length(Bn),length(KK),3))
}

a0 = alpha0[a]
for(rr in 1:length(Rho)){
  rho = Rho[rr]
  for(kk in 1:length(KK)){
    K = KK[kk]
    load(file = paste("test_rho",rho,"_K",K,"_a0",a0,".dat", sep = ""))
    reject = reject.boot = reject.boot.multi = array(0,dim=c(rep,length(Bn),3,K))
    reject.boot.simul = array(0,dim=c(rep,length(Bn),3))
    for(r in 1:rep){
      reject[r,,,] = results[[r]]$reject
      reject.boot[r,,,] = results[[r]]$reject.boot
      reject.boot.simul[r,,] = results[[r]]$reject.boot.simul
      reject.boot.multi[r,,,] = results[[r]]$reject.boot.multi
    }
    for(bb in 1:length(Bn)){
      rej.asy[[rr]][bb,kk,] = rowMeans(colMeans(reject[,bb,,]))
      rej.boot[[rr]][bb,kk,] = rowMeans(colMeans(reject.boot[,bb,,]))
      rej.boot.simul[[rr]][bb,kk,] = mean(colMeans(reject.boot.simul[,bb,]))
      rej.boot.multi[[rr]][bb,kk,] = rowMeans(colMeans(reject.boot.multi[,bb,,]))
    }
  }
}

round(rej.asy[[1]][,,2], digits = 4)
round(rej.boot[[1]][,,2], digits = 4)
round(rej.boot.simul[[1]][,,2], digits = 4)
round(rej.boot.multi[[1]][,,2], digits = 4)

round(rej.asy[[2]][,,2], digits = 4)
round(rej.boot[[2]][,,2], digits = 4)
round(rej.boot.simul[[2]][,,2], digits = 4)
round(rej.boot.multi[[2]][,,2], digits = 4)
