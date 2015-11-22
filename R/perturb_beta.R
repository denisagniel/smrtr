## perturb beta
#'
#' Get perturbed version of the nonparametric maximum likelihood estimator
#'
#' @param num.ptb Number of perturbations to perform
#' @param x n x p original design matrix
#' @param y List of length m containing outcomes
#' @param psd Pseudo-data generated from x and y based on marginal models
#' @param influence.fn Function to obtain the influence function of beta based on psd

#' @details
#'
#' Perturbed estimates are computed based on a beta distribution chosen to have 
#' first three moments equal to 1 (and then mean-shifted).
#'
#'@return n x num.ptb matrix with perturbed estimates
#'
#'@export
perturb.beta <- function(num.ptb, x, y, psd, influence.fn = po.influence) {
  n <- nrow(x)
  m <- length(y)
  c <- 4; a <- (c-2)/c; b <- (c-1)*a
  V <- matrix(rbeta(n*num.ptb, shape1=a, shape2=b)*c, nrow=num.ptb, ncol=n) - 1
  U <- influence.fn(psd, x, y)
  betahat <- psd$betahat
  beta.star <- t(betahat + t(V %*% U))
  return(beta.star)
}

po.influence <- function(psd, x, y, A.c = NULL) {
  k <- with(psd, sapply(alpha.inds, length))
  k.b0 <- with(psd, sapply(b0.inds, length))
  ps <- with(psd, sapply(beta.inds, length))
  m <- length(k)
  U <- NULL
  for (mm in 1:m) {
    this.k <- k[mm]
    p <- ps[mm]
    this.k.b0 <- k.b0[mm] + length(A.c[[mm]])
    p.sub <- p
    this.b0.ind <- NULL
    if (this.k.b0 > 0) {
      this.b0.ind <- with(psd, b0.inds[[mm]]) - this.k
    }
    this.b0.ind <- c(this.b0.ind, A.c[[mm]])
    this.beta.ind <- with(psd, beta.inds[[mm]]) - this.k
    this.beta.ind <- this.beta.ind[!this.beta.ind %in% this.b0.ind]

    U.a <- matrix(rep(0, (this.k + this.k.b0)*n),
                  ncol = this.k + this.k.b0)
    U.b <- matrix(rep(0, n*p.sub), ncol=p.sub)
    y.column <- y[,mm]
    beta.c <- with(psd, c(coefs[[mm]][b0.inds[[mm]]],
                          betahat[(mm-1)*p + 1:p]))
    alpha.c <- with(psd, alphahat[[mm]])
    ordered.y <- sort(unique(y.column))
    prob.gteq <- t(apply(x, 1, function(x) return(c(1, expit(alpha.c[1:this.k] + x %*% beta.c), 0))))
    for (j in 2:(this.k + 1)) {
      U.a[,j-1] <- (1 + exp(alpha.c[j-1] + x %*% beta.c))^-2 *
        exp(alpha.c[j-1] + x %*% beta.c) *
        ((y.column == ordered.y[j]) * (prob.gteq[,j] - prob.gteq[,j + 1])^-1 -
           (y.column == ordered.y[(j - 1)]) * (prob.gteq[,j - 1] - prob.gteq[,j])^-1)
    }

    if (this.k.b0 > 0) {
      for (ind in 1:this.k.b0) {
        for (ii in 1:n) {
          U.a[ii,(this.k + ind)] <- sum(
            (y.column[ii] == ordered.y) *
              (prob.gteq[ii,1:(this.k+1)] - prob.gteq[ii,2:(this.k+2)])^-1 *
              (c(0,
                 x[ii,this.b0.ind[ind]] *
                   exp(alpha.c[1:this.k] + x[ii,] %*% beta.c)/
                   (1 + exp(alpha.c[1:this.k] + x[ii,] %*% beta.c))^2) -
                 c(
                   x[ii,this.b0.ind[ind]] *
                     exp(alpha.c[1:this.k] + x[ii,] %*% beta.c)/
                     (1 + exp(alpha.c[1:this.k] + x[ii,] %*% beta.c))^2, 0)))
        }
      }
    }
    for (pp in 1:p.sub) {
      for (ii in 1:n) {
        U.b[ii,pp] <- sum((y.column[ii] == ordered.y) *
                            (prob.gteq[ii,1:(this.k+1)] -
                               prob.gteq[ii,2:(this.k+2)])^-1 *
                            (c(0,
                               x[ii,this.beta.ind[pp]] *
                                 exp(alpha.c[1:this.k] + x[ii,] %*% beta.c)/
                                 (1 + exp(alpha.c[1:this.k] + x[ii,] %*% beta.c))^2) -
                               c(x[ii,this.beta.ind[pp]] *
                                   exp(alpha.c[1:this.k] + x[ii,] %*% beta.c)/
                                   (1 + exp(alpha.c[1:this.k] + x[ii,] %*% beta.c))^2, 0)))
      }
    }
    B <- with(psd, inv.info[[mm]])
    A.inds <- with(psd, beta.inds[[mm]][!beta.inds[[mm]] %in% (A.c[[mm]] + this.k)])
    Ac.inds <- (1:(p + this.k))[-A.inds]
    B21 <- B[A.inds, Ac.inds]
    B22 <- B[A.inds, A.inds]
    if (ncol(U.a) != ncol(B21) | ncol(U.b) != nrow(B22)) browser()
    U <- cbind(U, U.a %*% t(B21) + U.b %*% B22)
  }
  U
}

expit <- function(x) exp(x)/(1 + exp(x))
