\name{calc.cpcerr}
\alias{calc.cpcerr}
\title{Calculate errors of CPC eigenvectors}
\description{
  Use Flury's formulas to calculate variance-covariance matrix of first
  CPC eigenvector components (or variances of all CPC eigenvector components)
}
\usage{
calc.cpcerr(x, f, cpcmat = NULL,  calc.cov = TRUE, use="complete.obs")
}
\arguments{
  \item{x}{a numeric matrix (or data frame) with all numeric values,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a factor describing the group structure of the data}
  \item{cpcmat}{CPC matrix}
  \item{calc.cov}{(logical) calculate variance-covariance matrix of
    elements of first eigenvector?}
  \item{use}{method for missing observations when computing covariances
    (see \code{\link{cov}} for details)}

}
\value{
  Either (if \code{calc.cov==TRUE}) the variance-covariance matrix of the
  elements of the first common principal component eigenvector, or (if
  \code{calc.cov==FALSE}) the variances of the elements of all of the
  CPC values
}
\author{Ben Bolker}
\examples{
  X=simdata(npts=200)
  calc.cpcerr(X)
  calc.cpcerr(X,calc.cov=FALSE)
}
\keyword{misc}
