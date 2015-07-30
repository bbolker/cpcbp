\name{covmat}
\alias{covmat}
\title{construct variance-covariance matrix}
\description{
  Generates a variance-covariance matrix from variances and correlations
}
\usage{
covmat(vars, cors, tol=1e-7)
}
\arguments{
  \item{vars}{a numeric vector of variances}
  \item{cors}{a numeric vector of correlations: replicated to the appropriate
    length if necessary}
  \item{tol}{numeric tolerance for checking positive definiteness}
}
\details{
  Creates a variance-covariance matrix with the specified variances and
  correlations.
}
\value{
  A square numeric matrix.
}
\author{Ben Bolker}
\note{Warns, but does not stop, if the matrix is non-positive-definite.}
\examples{
  covmat(c(1,2,3),c(0.4,-0.2,-0.1))
}
\keyword{array}
