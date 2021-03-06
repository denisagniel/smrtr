% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/sim_fns.R
\name{generate.y}
\alias{generate.y}
\title{Generate x values in {0,1,2}.}
\usage{
generate.y(x, beta.m, rho, n, m, k = n)
}
\arguments{
\item{x}{Matrix of predictor values.}

\item{beta.m}{Effect size matrix.}

\item{rho}{Correlation between outcomes.}

\item{n}{Sample size.}

\item{m}{Number of outcomes.}

\item{k}{Number of levels to coarsen y into.}
}
\value{
y matrix
}
\description{
Generate x values in {0,1,2}.
}

