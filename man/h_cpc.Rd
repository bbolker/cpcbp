\name{h_cpc}
\alias{h_cpc}
\title{Test whether all CPC are identical}
\usage{
  h_cpc(matList, n, B = new_cpc(matList, n))
}
\arguments{
  \item{matList}{list of matrices}

  \item{n}{numeric vector of sample sizes}

  \item{B}{CPC matrix}
}
\description{
  Tests whether the matrices share CPC based on equation
  2.18 (Flury, 1984). The null hypothesis is that all
  groups have the same PC.
}
\examples{
source("test_examples.R")
h_cpc(test_turtle,n_turtle)
}
\references{
  Flury, B.N.  Common principal components in K groups.
  J.Amer.Statist.Assoc.  1984, 79, 892-898.
}

