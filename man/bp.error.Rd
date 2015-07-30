\name{bp.error}
\alias{bp.error}
\title{calculate back-projection errors}
\description{
  Calculates the size of back-projection errors, based on errors of
  common principal components
}
\usage{
bp.error(x,f,cpcmat, m, eigvar,
use="complete.obs",debug=FALSE,
center=TRUE)
}

\arguments{
  \item{x}{a numeric matrix (or data frame with all numeric values,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a factor describing the group structure of the data}
  \item{cpcmat}{common principal components matrix (if already
    computed)}
  \item{m}{matrix of mean values of variables by group (if already computed)}
  \item{eigvar}{variances of eigenvalues (if already computed)}
  \item{use}{method for missing observations when computing covariances
    (see \code{\link{cov}} for details)}
  \item{debug}{(logical) print debugging information?}
  \item{center}{center group means on grand mean?}
}
\details{
  Uses CPC, variances of CPC (eigenvector) components, and means of
  groups to calculate total back-projection errors in variables.
  If \code{cpcmat} is not provided, calls \code{cpcvecfun} to
  compute CPC; if \code{eigvar} is not provided, calls
  \code{calc.cpcerr} to calculate variances on eigenvalues.
}
\value{
  If \code{new=TRUE}, a matrix of back-projection errors by
  variable and group; otherwise, 
  a vector of total back-projection errors in each variable.

}
\references{Flury}
\author{Ben Bolker}
\examples{
   set.seed(1001)
   X <- simdata(npts=100)
   bp.error(X)
}
\keyword{misc}
