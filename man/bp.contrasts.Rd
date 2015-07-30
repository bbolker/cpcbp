\name{bp.contrasts}
\alias{bp.contrasts}
\title{Estimate back-projected contrasts}
\description{
  Estimates contrasts with back-projection error incorporated;
  this can include estimates of the mean values of each group
}
\usage{
bp.contrasts(formula)
}
\arguments{
  \item{formula}{a formula}
}
\value{
  A table of estimated values, standard errors, t values, and p values
  for each variable, with the standard errors etc. adjusted to allow for
  back-projection error
}
\note{
  Use \code{x~f-1} (i.e., fit a model without an intercept) to get
  the values of each group rather than differences among them.
  It is the USER'S RESPONSIBILITY to check that the groups in the
  data actually share a common principal component!}
\references{Flury}
\author{Ben Bolker}
\examples{
X = simdata(offset=100)
## commented out because hanging in checks
\dontrun{
with(X,bp.contrasts(data~f))
with(X,bp.contrasts(data~f-1))
}
}
\keyword{misc}
