\name{pooled.cpc}
\alias{pooled.cpc}
\title{Compute CPC by mean-correcting each group}
\description{
Subtracts the mean from each group in the data set to shift
it to the origin, lumps the data, and computes principal components:
\emph{if} PCs for all groups are the same, this gives an estimate
of CPC
}
\usage{
pooled.cpc(x, f, use="complete.obs")
}
\arguments{
  \item{x}{a numeric matrix (or data frame with all numeric values,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a factor describing the group structure of the data}
  \item{use}{method for missing observations when computing covariances
    (see \code{\link{cov}} for details)}
}
\value{
  A matrix of the common principal components (eigenvector of the
  variance-covariance matrix) of the pooled data
}
\author{Ben Bolker}
\seealso{\code{\link{phillips.cpc}}}
\examples{
  X = simdata(npts=1000)
  pooled.cpc(X)
}
\keyword{array}
