\name{plot_multigrp}
\alias{plot_multigrp}
\title{plot grouped data}
\description{
  produces pair plots of multivariate, grouped data
}
\usage{
plot_multigrp(x, f, vars = 1:ncol(x),
cols,pchs,
    eqsc=FALSE,
    xlim=NULL, ylim=NULL, ...)
}
\arguments{
  \item{x}{a numeric matrix (or data frame with all numeric values,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a factor describing the group structure of the data}
  \item{vars}{variables to plot}
  \item{cols}{colors to represent groups}
  \item{pchs}{plotting characters to represent groups}
  \item{eqsc}{attempt equal x and y scales?}
  \item{xlim}{x limits}
  \item{ylim}{y limits}
  \item{\dots}{further arguments to \code{plot} and \code{points}}
}
\value{
  Produces a plot on the current device.
}
\author{Ben Bolker}
\examples{
   X = simdata(npts=20,offset=20,offset2=10)
   plot_multigrp(X, pch=16)
## this can also be done with lattice ...
   library(lattice)
   splom(~data,groups=f,data=X)
}
\keyword{hplot}
