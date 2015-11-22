#  fit_adaptive_hlasso
#'
#' Fit adaptive hierarchical lasso
#'
#'
#' @param x.t Design matrix or pseudo-design matrix where each column has been
#'  multiplied by abs(beta) corresponding to it
#' @param y Outcome or pseudo-outcome matrix
#' @param m Number of outcomes
#' @param n Number of observations
#' @param p Number of predictors
#' @param info.half Half Fisher information matrix for beta
#' @param BIC Logical, if FALSE, modified BIC criterion is used
#' @param max.iter Maximum number of iterations in iterative procedure
#' @param print Logical indicating whether to print intermediate updates
#'
#' @details
#'
#' Fits an adaptive version of the hierarchical lasso
#' suitable for multivariate mixed outcomes. Variable selection is performed
#' to both remove uninformative predictors and remove unimportant predictor-outcome
#' relationships. Tuning parameter selection is done via a BIC-style criterion.
#'
#'@return A list with the following elements:\itemize{
#'   \item \code{lasso.beta}: the non-hierarchical adaptive lasso estimate
#'   \item \code{hlasso.beta}: the adaptive hierarchical lasso estimate
#'   \item \code{num.iter}: the number of iterations used
#'   \item \code{d}: the estimate of d
#'   \item \code{alpha}: the estimate of alpha
#' }
#'
#'@seealso \code{\link[lars]{lars}}
#'
#'
#'@import lars
#'
#'@export

fit.adaptive.hlasso <- function(x.t, y, beta.ini, m, n, p, info.half, BIC = FALSE, max.iter = 50, print = FALSE) {
  #--------------------------------
  ## fit adaptive hierarchical lasso
  ##
  ## x.t : design matrix or pseudo-design matrix where each column has been
  ##  multiplied by abs(beta) corresponding to it
  ## y : outcome or pseudo-outcome matrix
  ## beta.ini : initial vector of estimates coming from NPMLE
  ## m : number of outcomes
  ## n : number of observations
  ## p : number of predictors
  ## info.half : half Fisher information matrix for beta
  ## BIC : logical, if FALSE, modified BIC criterion is used
  ## max.iter : maximum number of iterations in iterative procedure
  ## print : logical indicating whether to print intermediate updates
  #------------------------------------
  d <- rep(1, p)
  alpha <- get.alpha(x = x.t, y = y, d = d, n = n, beta.ini = beta.ini, info.half = info.half, BIC = BIC)
  lasso.beta <- alpha*abs(beta.ini)
  d <- get.d(x = x.t, y = y, alpha = alpha, m = m, p.x = p)

  phi.hat <- alpha*d
  phi.new <- rep(0, m*p)

  iter <- 0
  while(1/p/m*sum(abs(phi.hat - phi.new)) > 1/n & !all(d == 0) & iter < max.iter) {
    phi.hat <- phi.new
    phi.new <- 0
    alpha <- get.alpha(x = x.t, y = y, d = d, n = n, beta.ini = beta.ini, info.half = info.half, BIC = BIC)
    d <- get.d(x = x.t, y = y, alpha = alpha, m = m, p.x = p)
    phi.new <- alpha*d
    if (print == TRUE) print(sum(abs(phi.hat - phi.new)))
    iter <- iter + 1
  }
  return(list(lasso.beta = lasso.beta, hlasso.beta = d*alpha*abs(beta.ini),
              num.iter = iter, d = d, alpha = alpha))
}

get.alpha <- function(x, y, d, n, beta.ini, info.half, BIC = FALSE) {

  #########################################################
  ## x : mp x mp columns
  ## y : mp x 1 vector
  ## d : p x 1 vector
  ## n : scalar, sample size
  ## beta.ini : mp x 1 vector of initial coefficients
  ## info.half : half Fisher information matrix for beta
  ## BIC : logical, if FALSE, modified BIC criterion is used
  ##########################################################
  # browser()
  m <- length(y)/length(d)
  if (!((m %% 1) == 0)) stop("Dimensions of y and d do not match.")
  d.alpha <- diag(rep(d, m))       ## d matrix to find alpha
  x.alpha <- x %*% d.alpha         ## x * d
  # browser()
  fit.alpha <- lars::lars(x = x.alpha, y = y, type = 'lasso', intercept = FALSE, normalize = FALSE)

  grid.alpha <- seq(from=1, to=max(fit.alpha$df), length=500)
  df.alpha <- fit.alpha$df[floor(grid.alpha)]
  coef.alpha <- lars::predict.lars(fit.alpha, x.alpha, s=grid.alpha, type="coefficients", mode="step")
  neg2loglik.alpha <- apply(coef.alpha$coefficients, 1, function(alpha) {
    t(beta.ini - d*alpha*abs(beta.ini)) %*% info.half %*% info.half %*%
      (beta.ini - d*alpha*abs(beta.ini))
  })
  bic.alpha <- neg2loglik.alpha + ifelse(BIC, log(n), n^0.1) * df.alpha
  t.alpha <- grid.alpha[bic.alpha == min(bic.alpha)]

  return(alpha = coef(fit.alpha, s = t.alpha, mode='step'))
}

get.d <- function(x, y, alpha, m, p.x) {
  if (any(alpha != 0)) {
    alpha.d <- diag(alpha[1:p.x])
    for (i in 2:m) {
      alpha.d <- rbind(alpha.d, diag(alpha[((i-1)*p.x + 1):(i*p.x)]))
    }
    x.d <- x %*% alpha.d
    fit.d <- lars(x = x.d, y = y, type = 'lasso', intercept = FALSE, normalize = FALSE)

    return(coef(fit.d, s=1, mode='lambda'))
  }
  else return(rep(0, p.SNP))
}

