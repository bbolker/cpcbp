\name{h_sphere}
\alias{h_sphere}
\title{Hypothesis of Sphericity}
\usage{
  h_sphere(matList, n, q)
}
\arguments{
  \item{matList}{list of covariance matrices}

  \item{n}{numeric vector of sample sizes}

  \item{q}{number of eigenvectors to be kept}
}
\description{
  Tests whether to reject null hypothesis that the last p-q
  eigenvalues are equal, based on approximation 4.9
  (Flury,1986). The value p is the dimension of the
  matrices and the value q is the number of eigenvectors to
  be kept.
}
\examples{
source("test_examples.R")
h_sphere(test_turtle,n_turtle,q=1)
}
\references{
  Flury B.N.  Asymptotic theory for common principal
  component analysis.  Annals Statist.  1986.  14, 418-430.
}

