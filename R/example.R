#--------------------------------------
## example simulation run 
##  data-generating functions in sim_fns.R
#--------------------------------------
#----------------------------------------------
## setting parameters
n <- 500
p <- 30
m <- 4
rho <- 0.3
k <- 10
n.ptb <- 20
beta.m <- matrix(c(rep(c(1, 0), c(20, 10)), 
                   rep(c(0.5, 0), c(16, 14)), 
                   rep(c(1, 0), c(12, 18)),
                   rep(c(0.5, 0), c(8, 22))), ncol=4)
#----------------------------------------------

#----------------------------------------------
## generating data
x <- generate.x(n = n, p.low = 0.15, p.high = 0.15, p.x = p)
y.mat <- generate.y(x = x, beta.m = beta.m, rho = rho, n = n, m = m, k = k)
y <- lapply(data.frame(y.mat), function(x) x)
#----------------------------------------------

#----------------------------------------------
## estimation
psd <- fit.ini(x = x, y = y.mat, xname = 'snp', p = p)
betatilde <- psd$betahat
inv.info <- psd$inv.info
fit <- with(psd, fit.adaptive.hlasso(x.t = pseudoX.t, y = pseudoY, beta.ini = betahat, m = m, n = n, p = p, info.half = pseudoX, BIC = FALSE))
betahat <- fit$hlasso.beta
#----------------------------------------------

#-----------------------------------------------
## perturb
beta.star <- perturb.beta(num.ptb = n.ptb, y = y.mat, x = x, psd = psd)
#-----------------------------------------------

#-----------------------------------------------
## re-fit
betahat.star <- matrix(-9, n.ptb, m*p)
for (kk in 1:n.ptb) {
  print(kk)
  psX <- with(psd, pseudoX)
  psY.new <- psX %*% beta.star[kk,]
  psX.tnew <- psX %*% diag(abs(beta.star[kk,]))
  fit.new <- fit.adaptive.hlasso(x.t = psX.tnew, y = psY.new, 
                                 info.half = psX,
                                 beta.ini = beta.star[kk,], 
                                 m = m, n = n, p = p)
  betahat.star[kk,] <- fit.new$hlasso.beta
}
#-----------------------------------------------

#---------------------------------------------------
## test
rejected.tests <- smrt.stepdown(beta = betahat,
                                betastar = betahat.star,
                                inv.info = inv.info,
                                z.th = 0.95,
                                p = p,
                                m = m)
