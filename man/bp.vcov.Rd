\name{bp.vcov}
\alias{bp.vcov}
\title{Variance-covariance matrices for back-projection coefficients}
\description{
  Calculates variance-covariance matrices for back-projection
  coefficients for each variable for a CPC/Flury back-projection analysis
}
\usage{
bp.vcov(x, f, cpcmat, eigvar)
}
\arguments{
  \item{x}{a numeric matrix (or data frame) with all numeric values,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a factor describing the group structure of the data}
  \item{cpcmat}{CPC matrix (calculated if not specified)}
  \item{eigvar}{variances of eigenvector elements (calculated if not
    specified)}
}
\value{
  An n by n by n array (where n is the number of variables/columns of
  x), where [i,j,k] gives the jk element of the covariance matrix for
  the ith variable
}
\author{Ben Bolker}
\examples{
bp.vcov(simdata())
}
\keyword{misc}
