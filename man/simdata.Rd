\name{simdata}
\alias{simdata}
\title{Simulate data for back-projection exercises}
\description{
  Creates identically shaped but offset multivariate normal groups
}
\usage{
simdata(vars = c(10, 10, 10), cors = 0.8, npts = 200, seed, m1a = rep(0, length(vars)), offset = 1, offset2 = 0)
}
\arguments{
  \item{vars}{within-group variances}
  \item{cors}{within-group correlations}
  \item{npts}{number of points in each group (will be repeated as necessary)}
  \item{seed}{random-number seed}
  \item{m1a}{means of first group}
  \item{offset}{offset(s) in the first PC direction}
  \item{offset2}{offset(s) in the second PC direction}
}
\value{
  a data frame containing the "data" for each group, followed
  by a factor specifying the group
}
\author{Ben Bolker}
\examples{
   X = simdata(offset=10,offset2=3)
   plot_multigrp(X)
}
\keyword{misc}
