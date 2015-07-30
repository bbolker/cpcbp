\name{h_many}
\alias{h_many}
\title{Comparing CPC to set of eigenvectors}
\usage{
  h_many(matList, n, betaList)
}
\arguments{
  \item{matList}{list of matrices}

  \item{n}{numeric vector of sample sizes}

  \item{betaList}{list of eigenvectors (organized from
  largest to smallest eigenvalue) to be compared to7}
}
\description{
  Tests whether to reject the null hypothesis that a set of
  q hypothetical eigenvectors are equal to the first q
  likelihood eigenvectors in CPC, based on equation 3.3
  (Flury, 1986). The value q must be a value of 2 or
  greater.  For q=1, see the function h_one.
}
\references{
  Flury B.N.  Asymptotic theory for common principal
  component analysis.  Annals Statist.  1986.  14, 418-430.
}

