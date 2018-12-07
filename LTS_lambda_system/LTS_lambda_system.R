libraries = c("mvtnorm", "hdm", "sandwich", "matrixStats", "doSNOW")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
    install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

#c1 = makeCluster(4)  # core numbers
#registerDoSNOW(c1)
set.seed(2017)

# data generating (systems)
T = 100  #sample size
P = c(50, 100, 150)  #number of equations
KK = c(50, 100, 150)  #number of covariates 
b = 10  #nonzero entries
d = 5
nboot = 500
c = rep(0, nboot)
rep = 100
ratios.mean = matrix(0, length(P), 2)
ratios.median = matrix(0, length(P), 2)
ratios.mean2 = matrix(0, length(P), 2)
ratios.median2 = matrix(0, length(P), 2)

for (pp in 1:length(P)) {
    p = P[pp]
    K = KK[pp]
    cov = matrix(0, K, K)
    for (i in 1:K) {
        for (j in 1:K) {
            cov[i, j] = 0.5^(abs(i - j))
        }
    }
    beta.matrix = matrix(0, p, K)
    lambda.ga = 2 * 1.1 * qnorm(1 - 0.1/(2 * K * p)) * sqrt(T)
    for (i in 1:p) {
        indx = ceiling(i/d)
        beta.matrix[i, (d * (indx - 1) + 1):(indx * d)] = b
    }
    Y.data = matrix(0, T, p)
    load = matrix(0, p, K)
    Smatrix.boot = array(0, dim = c(p, K, nboot))
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
        X = rmvnorm(T, rep(0, K), cov)
        # X = matrix(rnorm(T*K), ncol=K)
        
        # heteroscadastic case
        for (i in 1:p) {
            Y = X %*% beta.matrix[i, ] + rnorm(T, 0, 1)  #sd=1 or 0.5
            Y.data[, i] = Y
            fit1 = rlasso(X, Y, penalty = list(homoscedastic = "none", lambda.start = 2 * 
                0.5 * sqrt(T) * qnorm(1 - 0.1/(2 * K))), post = FALSE)
            beta = fit1$beta
            intercept = fit1$intercept
            res = Y - X %*% beta - intercept * rep(1, length(Y))
            for (j in 1:nboot) {
                res.boot = rnorm(T)
                for (k in 1:K) {
                  Smatrix.boot[i, k, j] = sum(res.boot * X[, k] * res)/sqrt(T)
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
            Y = Y.data[, i]
            fit_hdm = rlasso(X, Y, penalty = list(X.dependent.lambda = TRUE, numSim = nboot))
            fit_hdm2 = rlasso(X, Y, penalty = list(X.dependent.lambda = FALSE))
            rmse[i, 1] = sqrt(mean((X %*% beta.matrix[i, ] - predict(fit_hdm))^2))
            l2[i, 1] = sqrt(sum((fit_hdm$beta - beta.matrix[i, ])^2))
            rmse[i, 2] = sqrt(mean((X %*% beta.matrix[i, ] - predict(fit_hdm2))^2))
            l2[i, 2] = sqrt(sum((fit_hdm2$beta - beta.matrix[i, ])^2))
            fit_hdm.joint = rlasso(X, Y, penalty = list(homoscedastic = "none", lambda.start = lambda.boot))
            fit_hdm.joint2 = rlasso(X, Y, penalty = list(homoscedastic = "none", 
                lambda.start = lambda.ga))
            rmse.joint[i, 1] = sqrt(mean((X %*% beta.matrix[i, ] - predict(fit_hdm.joint))^2))
            l2.joint[i, 1] = sqrt(sum((fit_hdm.joint$beta - beta.matrix[i, ])^2))
            rmse.joint[i, 2] = sqrt(mean((X %*% beta.matrix[i, ] - predict(fit_hdm.joint2))^2))
            l2.joint[i, 2] = sqrt(sum((fit_hdm.joint2$beta - beta.matrix[i, ])^2))
        }
        ratio.mean = rbind(ratio.mean, colMeans(rmse.joint)/colMeans(rmse))
        ratio.median = rbind(ratio.median, colMedians(rmse.joint)/colMedians(rmse))
        ratio.mean2 = rbind(ratio.mean2, colMeans(l2.joint)/colMeans(l2))
        ratio.median2 = rbind(ratio.median2, colMedians(l2.joint)/colMedians(l2))
    }
    ratios.mean[pp, ] = colMeans(ratio.mean[-1, ])
    ratios.median[pp, ] = colMeans(ratio.median[-1, ])
    ratios.mean2[pp, ] = colMeans(ratio.mean2[-1, ])
    ratios.median2[pp, ] = colMeans(ratio.median2[-1, ])
}
ratios.mean
ratios.median
ratios.mean2
ratios.median2

results = list(ratios.mean = ratios.mean, ratios.median = ratios.median, ratios.mean2 = ratios.mean2, 
    ratios.median2 = ratios.median2)
save(results, "IS-system.dat")
