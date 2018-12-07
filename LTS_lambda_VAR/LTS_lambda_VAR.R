libraries = c("mvtnorm", "hdm", "Matrix", "sandwich", "matrixStats", "doSNOW")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
    install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

#c1 = makeCluster(4)  # core numbers
#registerDoSNOW(c1)
set.seed(2017)

# data generating (VAR(1))
T = 100  #sample size
P = c(50, 100, 150)  #number of variables
Bn = c(4, 10, 20, 25)  #length of block
lag = 1
nboot = 500
c = rep(0, nboot)
rep = 100
Rho = c(0.05, 0.15)
d = 5
ratios.mean = list()
ratios.median = list()
ratios.mean2 = list()
ratios.median2 = list()
for (rr in 1:length(Rho)) {
    ratios.mean[[rr]] = array(0, dim = c(length(Bn), length(P), 2))
    ratios.median[[rr]] = array(0, dim = c(length(Bn), length(P), 2))
    ratios.mean2[[rr]] = array(0, dim = c(length(Bn), length(P), 2))
    ratios.median2[[rr]] = array(0, dim = c(length(Bn), length(P), 2))
}

for (rr in 1:length(Rho)) {
    rho = Rho[rr]
    for (bb in 1:length(Bn)) {
        bn = Bn[bb]
        ln = T/bn  #number of blocks
        for (kk in 1:length(P)) {
            p = P[kk]
            K = p
            lambda.ga = 2 * 1.1 * qnorm(1 - 0.1/(2 * K * p)) * sqrt(T)
            Phi = bdiag(matrix(rho, d, d))
            for (pp in 1:(K/d - 1)) {
                Phi = bdiag(Phi, matrix(rho, d, d))
            }
            Phi = as.matrix(Phi)
            
            Smatrix.boot = array(0, dim = c(p, K, nboot))
            load = matrix(0, p, K)
            Sample = matrix(0, nrow = (T + lag), ncol = p)
            ratio.mean = 0
            ratio.median = 0
            ratio.mean2 = 0
            ratio.median2 = 0
            rmse = matrix(0, p, 2)
            rmse.joint = matrix(0, p, 2)
            l2 = matrix(0, p, 2)
            l2.joint = matrix(0, p, 2)
            for (r in 1:rep) {
                # foreach(r=1:rep,.packages=c('mvtnorm', 'Matrix', 'hdm', 'sandwich',
                # 'matrixStats')) %dopar% {
                innov = rmvnorm((T + lag), rep(0, p), diag(p))
                Sample[1, ] = innov[1, ]
                for (t in 2:(T + lag)) {
                  Sample[t, ] = Phi %*% Sample[t - lag, ] + innov[t, ]
                }
                
                # heteroscadastic case
                for (i in 1:p) {
                  X = Sample[1:T, ]
                  Y = Sample[(1 + lag):(T + lag), i]
                  fit1 = rlasso(X, Y, penalty = list(homoscedastic = "none", lambda.start = 2 * 
                    0.5 * sqrt(T) * qnorm(1 - 0.1/(2 * K))), post = FALSE)
                  beta = fit1$beta
                  res = Y - X %*% beta - fit1$intercept
                  for (j in 1:nboot) {
                    res.boot = rnorm(ln)
                    for (k in 1:K) {
                      sum = 0
                      for (l in 1:ln) {
                        sum = sum + sum(X[((l - 1) * bn + 1):(l * bn), k] * res[((l - 
                          1) * bn + 1):(l * bn)]) * res.boot[l]
                      }
                      Smatrix.boot[i, k, j] = sum/sqrt(T)
                    }
                  }
                  #### compute the penalty loadings
                  for (k in 1:K) {
                    load[i, k] = sqrt(lrvar(X[, k] * res) * T)
                  }
                }
                
                for (j in 1:nboot) {
                  c[j] = max(abs(Smatrix.boot[, , j]/load))
                }
                lambda.boot = 2 * quantile(c, 0.9) * sqrt(T) * 1.1
                
                # compare the in-sample fitting performance
                for (i in 1:p) {
                  X = Sample[1:T, ]
                  Y = Sample[(1 + lag):(T + lag), i]
                  fit_hdm = rlasso(X, Y, penalty = list(X.dependent.lambda = TRUE, 
                    numSim = nboot))
                  fit_hdm2 = rlasso(X, Y, penalty = list(X.dependent.lambda = FALSE))
                  rmse[i, 1] = sqrt(mean((X %*% Phi[i, ] - predict(fit_hdm))^2))
                  l2[i, 1] = sqrt(sum((fit_hdm$beta - Phi[i, ])^2))
                  rmse[i, 2] = sqrt(mean((X %*% Phi[i, ] - predict(fit_hdm2))^2))
                  l2[i, 2] = sqrt(sum((fit_hdm2$beta - Phi[i, ])^2))
                  fit_hdm.joint = rlasso(X, Y, penalty = list(homoscedastic = "none", 
                    lambda.start = lambda.boot))
                  fit_hdm.joint2 = rlasso(X, Y, penalty = list(homoscedastic = "none", 
                    lambda.start = lambda.ga))
                  rmse.joint[i, 1] = sqrt(mean((X %*% Phi[i, ] - predict(fit_hdm.joint))^2))
                  l2.joint[i, 1] = sqrt(sum((fit_hdm.joint$beta - Phi[i, ])^2))
                  rmse.joint[i, 2] = sqrt(mean((X %*% Phi[i, ] - predict(fit_hdm.joint2))^2))
                  l2.joint[i, 2] = sqrt(sum((fit_hdm.joint2$beta - Phi[i, ])^2))
                }
                ratio.mean = rbind(ratio.mean, colMeans(rmse.joint)/colMeans(rmse))
                ratio.median = rbind(ratio.median, colMedians(rmse.joint)/colMedians(rmse))
                ratio.mean2 = rbind(ratio.mean2, colMeans(l2.joint)/colMeans(l2))
                ratio.median2 = rbind(ratio.median2, colMedians(l2.joint)/colMedians(l2))
            }
            ratios.mean[[rr]][bb, kk, ] = colMeans(ratio.mean[-1, ])
            ratios.median[[rr]][bb, kk, ] = colMeans(ratio.median[-1, ])
            ratios.mean2[[rr]][bb, kk, ] = colMeans(ratio.mean2[-1, ])
            ratios.median2[[rr]][bb, kk, ] = colMeans(ratio.median2[-1, ])
        }
    }
}

round(ratios.mean[[1]], digits = 4)
round(ratios.median[[1]], digits = 4)
round(ratios.mean2[[1]], digits = 4)
round(ratios.median2[[1]], digits = 4)

round(ratios.mean[[2]], digits = 4)
round(ratios.median[[2]], digits = 4)
round(ratios.mean2[[2]], digits = 4)
round(ratios.median2[[2]], digits = 4)

results = list(ratios.mean = ratios.mean, ratios.median = ratios.median, ratios.mean2 = ratios.mean2, 
    ratios.median2 = ratios.median2)
save(results, "IS-VAR.dat")
