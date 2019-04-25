libraries = c("mvtnorm", "hdm", "sandwich", "matrixStats","doSNOW")
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

sim_scene = function(nboot,T,p,K,d,b){ 
  cov = matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      cov[i,j] = 0.5^(abs(i-j))
    }
  }
  beta.matrix = matrix(0, p, K)
  lambda.ga = 2*1.1*qnorm(1-0.1/(2*K*p))*sqrt(T)
  for (i in 1:p){
    indx = ceiling(i/d)
    beta.matrix[i, (d*(indx-1)+1):(indx*d)] = b
  }
  Y.data = matrix(0, T, p)
  load = matrix(0,p,K)
  Smatrix.boot = array(0,dim=c(p,K,nboot))
  c = rep(0, nboot)
  
  X = rmvnorm(T, rep(0,K), cov)
  for (i in 1:p){
    Y = X%*%beta.matrix[i,] + rnorm(T, 0, 1) #sd=1 or 0.5
    Y.data[,i] = Y
    fit1 =  rlasso(X, Y, penalty = list(homoscedastic = "none", lambda.start = 2*0.5*sqrt(T)*qnorm(1-0.1/(2*K))), post=FALSE)
    beta = fit1$beta
    intercept = fit1$intercept
    res = Y - X %*% beta - intercept * rep(1, length(Y))
    for (j in 1:nboot){
      res.boot = rnorm(T) 
      for (k in 1:K){
        Smatrix.boot[i,k,j] = sum(res.boot*X[,k]*res)/sqrt(T)
      }
    }
    #### compute the penalty loadings
    for(k in 1:K){
      load[i,k] = sd(X[,k]*res)
    }
  }
  
  for(j in 1:nboot){
    c[j]=max(abs(Smatrix.boot[,,j]/load))
  }
  lambda.boot = 2*quantile(c, 0.9)*sqrt(T)*1.1
  
  rmse = rmse.joint = l2 = l2.joint = matrix(0,p,2)
  # compare the in-sample fitting performance
  for (i in 1:p){
    Y = Y.data[,i]
    fit_hdm =  rlasso(X, Y, penalty = list(X.dependent.lambda = TRUE, numSim = nboot))
    fit_hdm2 =  rlasso(X, Y, penalty = list(X.dependent.lambda = FALSE))
    rmse[i,1] = sqrt(mean((X%*%beta.matrix[i,]-X%*%fit_hdm$beta)^2))
    l2[i,1] = sqrt(sum((fit_hdm$beta-beta.matrix[i,])^2))
    rmse[i,2] = sqrt(mean((X%*%beta.matrix[i,]-X%*%fit_hdm2$beta)^2))
    l2[i,2] = sqrt(sum((fit_hdm2$beta-beta.matrix[i,])^2))
    fit_hdm.joint =  rlasso(X,Y,penalty = list(homoscedastic = "none", lambda.start = lambda.boot))
    fit_hdm.joint2 =  rlasso(X,Y,penalty = list(homoscedastic = "none", lambda.start = lambda.ga))
    rmse.joint[i,1] = sqrt(mean((X%*%beta.matrix[i,]-X%*%fit_hdm.joint$beta)^2))
    l2.joint[i,1] = sqrt(sum((fit_hdm.joint$beta-beta.matrix[i,])^2))
    rmse.joint[i,2] = sqrt(mean((X%*%beta.matrix[i,]-X%*%fit_hdm.joint2$beta)^2))
    l2.joint[i,2] = sqrt(sum((fit_hdm.joint2$beta-beta.matrix[i,])^2))
  }
  
  list(ratio.median = colMedians(rmse.joint)/colMedians(rmse), ratio.median2 = colMedians(l2.joint)/colMedians(l2))
}

for(pp in 1:length(P)){
  p = P[pp]
  K = KK[pp]
  results = foreach(l=1:rep, .packages=c("mvtnorm", "Matrix", "hdm", "sandwich", "matrixStats"), .inorder=FALSE) %dopar%{ 
    sim_scene(nboot=nboot,T=T,p=p,K=K,d=d,b=b)
  }
  save(results, file = paste("ratios_K",K,".dat", sep = ""))
}

ratios.mean = ratios.median = ratios.mean2 = ratios.median2 = ratios.sd = ratios.sd2 = matrix(0, length(P),2)
for(pp in 1:length(P)){
  K = KK[pp]
  load(file = paste("ratios_K",K,".dat", sep = ""))
  ratio.median = ratio.median2 = numeric(0)
  for(r in 1:rep){
    ratio.median = rbind(ratio.median,results[[r]]$ratio.median)
    ratio.median2 = rbind(ratio.median2,results[[r]]$ratio.median2)
  }
  ratios.mean[pp,] = colMeans(ratio.median)
  ratios.median[pp,] = colMedians(ratio.median)
  ratios.mean2[pp,] = colMeans(ratio.median2)
  ratios.median2[pp,] = colMedians(ratio.median2)
  ratios.sd[pp,] = colSds(ratio.median)
  ratios.sd2[pp,] = colSds(ratio.median2)
}

round(ratios.mean[,1], digits = 4)
round(ratios.median[,1], digits = 4)
round(ratios.sd[,1], digits = 4)
round(ratios.mean2[,1], digits = 4)
round(ratios.median2[,1], digits = 4)
round(ratios.sd2[,1], digits = 4)
