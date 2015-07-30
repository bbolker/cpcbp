\name{sim.theor}
\alias{sim.theor}
\title{Generate theoretical values for back-projection exercises}
\description{
  Returns a list of the theoretical group means, common
  variance-covariance matrix, and common eigenvectors for a specified
  set of back-projection simulation parameters
}
\usage{
sim.theor(vars = c(10, 10, 10), cors = 0.8, m1a = rep(0, length(vars)), offset = 1, offset2 = 0)
}
\arguments{
  \item{vars}{within-group variances}
  \item{cors}{within-group correlations}
  \item{m1a}{means of first group}
  \item{offset}{offset(s) in the first PC direction}
  \item{offset2}{offset(s) in the second PC direction}
}
\value{
  \item{mean}{the mean values of each group}
  \item{varcov}{the common variance-covariance matrix}
  \item{eigs}{the common eigenvector matrix}
}
\seealso{\code{\link{simdata}}}
\author{Ben Bolker}
\examples{
   sim.theor(offset=10,offset2=3)
}
\keyword{misc}
