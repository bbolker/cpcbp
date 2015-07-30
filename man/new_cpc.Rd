\name{new_cpc}
\alias{new_cpc}
\title{FG Algorithm}
\usage{
  new_cpc(matList, n = rep(1, length(matList)),
    verbose = FALSE, eps.F = 1e-06, eps.G = 1e-06,
    Bstart = diag(p))
}
\arguments{
  \item{matList}{list of covariance matrices}

  \item{n}{numeric vector of sample sizes}

  \item{eps.F}{accuracy constant for determining when to
  stop the F-level}

  \item{eps.G}{accuracy constant for determining when to
  stop the G-level}

  \item{Bstart}{initial B approximation (by default
  diag(p))}
}
\description{
  Applies the FG Algorithm (Flury and Gautschi, 1986) to a
  set of matrices to be simultaeneously orthogonally
  diagonalized by a common matrix B.  The matrix B contains
  the common principal components (CPC).
}
\examples{
source("test_examples.R")
new_cpc(test_iris)
}
\references{
  Flury B.N.; Gautschi, W.  An algorithm for simultaneous
  orthogonal transformation of several positive definite
  symmetric matrices to nearly diagonal form.  SIAM J. Sci.
  Stat. Comput.  1986, 7, 169-84
}

