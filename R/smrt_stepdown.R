# smrt.stepdown
#'
#' Perform stepdown testing on parameter estimates for each predictor
#' 
#' @param beta Vector of estimates
#' @param betastar Matrix of perturbed estimates
#' @param inv.info Full Fisher information matrix
#' @param inv.info.beta Sub-matrix of Fisher information corresponding to \code{beta}.
#' @param z.th Threshold for rejection, either scalar or p-length vector.
#' @param get.u Function to obtain the influence function of beta based on \code{psd}.

#' @details
#'
#' Estimates should be arranged such that all estimates for first outcome
#' come first, then all estimates for second outcomes, and so on. Only one of
#' inv.info and inv.info.beta need to be supplied. 
#'
#'@return Indices of beta corresponding to rejected tests.
#'
#'@export


smrt.stepdown <- function(beta, betastar, inv.info = NULL, inv.info.beta = NULL, z.th, p, m) {
  if (is.null(inv.info) & is.null(inv.info.beta)) {
    stop('please provide inverse information.')
  } else if (!is.null(inv.info.beta)) {
    inv.info <- inv.info.beta
  } else {
    inv.info <- lapply(inv.info, function(x) {
      p.prime <- dim(x)[1]
      x[(p.prime - p + 1):p.prime, (p.prime - p + 1):p.prime]
    })
  }
  
  if (length(z.th) == 1) {
    z.th <- rep(z.th, p)
  } else if (length(z.th) != p) {
    stop('number of thresholds specified must be equal to p')
  }
  rej.set <- c()
  sigma <- sqrt(unlist(lapply(inv.info, diag)))
  z <- abs(beta/sigma)
  z.star <- abs(t((t(betastar) - beta)/sigma))
  
  for (i in 1:p) {
    this.snp <- matrix(z, ncol = m)[i,]
    this.snpstar <- t(apply(z.star, 1, function(x) matrix(x, ncol = m)[i,]))
    this.snp.rej <- stepdown(this.snp, this.snpstar, z.th[i])
    if (!is.null(this.snp.rej)) {
      rej.set <- c(rej.set, get.x.null(i, p, m)[this.snp.rej])
    }
  }
  rej.set
}

#' Basic stepdown function
#' 
#' @param z Index vector of estimates
#' @param zstar Matrix of perturbed estimates
#' @param thr Threshold for rejection
#' @return Indices corresponding to rejected tests.
#'
#'@export

stepdown <- function(z, zstar, thr) {
  rej.set <- c()
  names(z) <- 1:length(z)
  for (i in 1:length(z)) {
    max.ind <- which.max(z)
    if (length(z) > 1) {
      max.z <- max(z)
      max.zstar <- apply(zstar, 1, max)
    } else {
      max.z <- z
      max.zstar <- zstar
    }
    if (max.z <= quantile(max.zstar, thr, na.rm = TRUE)) {
      return(rej.set)
    } else {
      rej.set <- c(rej.set, as.numeric(names(z)[z == max.z]))
      if (length(z) > 1) {
        z <- z[-max.ind]
        zstar <- zstar[,-max.ind]
      }
    }
  }
  rej.set
}

get.x.null <- function(inds, p, m) {
  if (all(inds == 0)) {
    NULL
  } else if (any(inds == 0)) {
    stop('if 0 is an index, it must be the only index.')
  } else {  
    sort(c(sapply(inds, function(x) x + 0:(m-1)*p)))
  }
}

