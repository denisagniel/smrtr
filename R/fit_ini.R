#  fit_ini
#'
#' Produce pseudo x and y based on marginal proportion odds model fits and quadratic
##   expansion of the (profile) likelihood
#'
#'
#' @param x Named matrix of predictors
#' @param y n x m matrix of outcomes
#' @param xname String to identify predictors in x
#' @param yname String to identify outcomes in y
#' @param p Number of predictors (doesn't need to be same as number of columns of x)
#'
#'
#'
#' #'@return A list with the following elements:\itemize{
#'   \item \code{pseudoY}: pseudo-outcome for quadratic expansion of the likelihood
#'   \item \code{betahat}: coefficients corresponding to the NPMLE of finite-dimensional parameters
#'   \item \code{alphahat}: coefficients corresponding to the NPMLE of infinite dimensional parameters
#'   \item \code{pseudoX}: pseudo-predictors for quadratic expansion
#'   \item \code{pseudoX.t}: adaptive version of pseudo-predictors
#'   \item \code{inv.info}: inverse Fisher information from  marginal models
#'   \item \code{inv.info.beta}: just the part of the inverse information corresponding to beta
#'   \item \code{beta.inds}: indices in the coefficient vectors corresponding to beta
#'   \item \code{alpha.inds}: indices in the coefficient vectors corresponding to the infinite-dimensional parameter
#'   \item \code{b0.inds}: indices in the coefficient vectors corresponding to nuisance covariates
#'   \item \code{coefs}: full coefficient vector
#' }
#'
#'
#'
#'@importFrom rms lrm
#'
#'@export
#'
fit.ini <- function(x, y, xname, yname = 'y', p) {
  pseudoY = pseudoX = pseudoX.t = betahat = alphahat = inv.info = coefs =
    beta.inds = alpha.inds = b0.inds = inv.info.beta = NULL
  m <- ncol(y)
  for (mm in 1:m) {
    this.y <- y[,mm]
    fit <- lrm(as.formula('this.y ~ x'), data = data.frame(this.y, x))
    coef.ini <- fit$coef
    beta.ind <- grep(xname, names(coef.ini))
    if ('Intercept' %in% names(coef.ini)) {
      alpha.ind <- grep('Intercept', names(coef.ini))
    } else {
      alpha.ind <- grep(yname, names(coef.ini))
    }
    b0.ind <- (1:length(coef.ini))[-c(beta.ind, alpha.ind)]
    tmpalpha = coef.ini[alpha.ind]
    beta0 <- coef.ini[b0.ind]
    tmpbeta = coef.ini[beta.ind]
    tmpinfo = fit$info.matrix
    tmpind <- c(alpha.ind, b0.ind)
    info = tmpinfo[-tmpind,-tmpind] - tmpinfo[-tmpind,tmpind]%*%
      solve(tmpinfo[tmpind,tmpind])%*%tmpinfo[tmpind,-tmpind]
    junk = svd(info)
    info.half = junk$u%*%diag(sqrt(junk$d))%*%t(junk$u)
    tmpY = info.half%*%tmpbeta
    tmpX = info.half
    tmpX.t = tmpX %*% diag(abs(tmpbeta))
    pseudoY   = rbind(pseudoY, tmpY)
    betahat  = c(betahat, tmpbeta)
    alphahat[[mm]] = tmpalpha
    pseudoX   = cbind(pseudoX,  kronecker(diag(m),tmpX)[  ,1:p + (mm-1)*p])
    pseudoX.t = cbind(pseudoX.t,kronecker(diag(m),tmpX.t)[,1:p + (mm-1)*p])
    inv.info[[mm]] = solve(tmpinfo)
    inv.info.beta[[mm]] <- solve(info)
    beta.inds[[mm]] <- beta.ind
    alpha.inds[[mm]] <- alpha.ind
    b0.inds[[mm]] <- b0.ind
    coefs[[mm]] <- coef.ini
  }
  list(pseudoY = pseudoY, betahat = betahat, alphahat = alphahat,
       pseudoX = pseudoX, pseudoX.t = pseudoX.t, inv.info = inv.info,
       inv.info.beta = inv.info.beta, beta.inds = beta.inds,
       alpha.inds = alpha.inds, b0.inds = b0.inds, coefs = coefs)
}
