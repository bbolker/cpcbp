\name{h_one}
\alias{h_one}
\title{Comparing first CPC to one hypothetical eigenvector}
\usage{
  h_one(matList, n, v)
}
\arguments{
  \item{matList}{list of matrices}

  \item{n}{numeric vector of sample sizes}

  \item{v}{eigenvector to be compared to}
}
\description{
  Tests whether to reject the null hypothesis that the CPC
  with the largest eigenvalue is similar to a hypothetical
  eigenvector, based on equation 3.4 (Flury, 1986).  If
  comparing 2 or more eigenvectors, see the function
  h_many.
}
\examples{
source("test_examples.R")
v=rep((1/sqrt(3)),3)
h_one(test_turtle,n_turtle,v)
}
\references{
  Flury B.N.  Asymptotic theory for common principal
  component analysis.  Annals Statist.  1986.  14, 418-430.
}

