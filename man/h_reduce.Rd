\name{h_reduce}
\alias{h_reduce}
\title{Reducing Parameter Space}
\usage{
  h_reduce(matList, n, q, f, alpha = 0.05)
}
\arguments{
  \item{matList}{list of covariance matrices}

  \item{n}{numeric vector of sample sizes}

  \item{q}{number of eigenvectors to be kept}

  \item{f}{a fraction value between 0 and 1 which
  determines the relative contribution of an eigenvalue to
  CPC}

  \item{alpha}{confidence level}
}
\description{
  Tests whether to reject the NH that the last p-q
  eigenvalues are very small and approximately zero
  variance, based on approximation 4.3 (Flury, 1986).
}
\references{
  Flury B.N.  Asymptotic theory for common principal
  component analysis.  Annals Statist.  1986.  14, 418-430.
}

