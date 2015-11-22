# sim fns
#' Generate x values in \{0,1,2\}.
#' 
#' @param n Sample size.
#' @param p.low Lower bound for minor allele frequency.
#' @param p.high Upper bound for minor allele frequency.
#' @param p.x Number of predictors.
#' @return x matrix
#'
#'@export

generate.x <- function(n, p.low, p.high, p.x) {
  x.mat <- c()
  for (i in 1:p.x) {
    if (p.low == p.high) {
      p <- p.low
    } else p <- runif(1)*(p.high-p.low) + p.low
    xx <- rep(1, n)
    xx <- sample(2:0, prob=c(p^2, 2*p*(1-p), (1-p)^2),
                 size = n, replace=TRUE)
    x.mat <- cbind(x.mat, xx)
  }
  colnames(x.mat) <- paste('snp', 1:p.x, sep='.')
  return(x.mat)
}

rank.unique.vec <- function(v, k, ties = FALSE, ...) {
  if (!ties) {
    n <- length(v)
    inds <- rep(floor(n/k), k)
    inds.2 <- c(rep(1, n %%k), rep(0, k - n %% k))
    r <- rep(1:k, inds + inds.2)
    return(r[order(order(v))])
  } else {
    q.v <- quantile(v, (1:(k-1))/k, ...)
    return(sapply(v, function(x) sum(x > q.v))+ 1)
  }
}

#' Generate x values in {0,1,2}.
#' 
#' @param x Matrix of predictor values.
#' @param beta.m Effect size matrix.
#' @param rho Correlation between outcomes.
#' @param n Sample size.
#' @param m Number of outcomes.
#' @param k Number of levels to coarsen y into.
#' @return y matrix
#'
#'@export

generate.y <- function(x, beta.m, rho, n, m, k = n) {
  p <- dim(beta.m)[1]
  z <- matrix(rnorm(n*m), n, m)*sqrt((1-rho)) + sqrt(rho)*(rnorm(n))
  u <- apply(z, 2, pnorm)
  e <- log(u) - log(1-u)
  y <- exp(x %*% beta.m + e)
  yy <- apply(y, 2, rank.unique.vec, k = k)
}
