\name{nullsim}
\alias{nullsim}
\alias{nullFstats}
\title{Simulate null hypothesis (size differences only)}
\description{
  Analyze multivariate data and simulate (a) null set(s) of
  data assuming that the different groups have a common
  first principal component and that their means differ only
  along the CPC1 axis; extract F statistics
}
\usage{
nullsim(x, f, P1, plot.it = FALSE)
nullFstats(x,f,n=1,bpanova=FALSE)
}
\arguments{
  \item{x}{a numeric data matrix}
  \item{f}{grouping variable (factor)}
  \item{P1}{the result of running \code{\link{phillips.cpc}}
    on the data; if missing, will be recomputed}
  \item{plot.it}{produce a descriptive plot? (STUB)}
  \item{n}{number of replicates}
  \item{bpanova}{(logical) use \code{\link{bp.anova}} to compute F
    statistics (alternatively use \code{\link{anova}} on back-projected
    data, without BP error correction)}
}
\value{
  \code{nullsim}: a numeric matrix the same shape as \code{x}.
  \code{nullFstats}: a matrix of replicate F statistics for
  the null case.
  
}
\author{Ben Bolker}
\examples{
  s1 = simdata(offset=10,offset2=10)
  plot_multigrp(s1)
  s1n = nullsim(s1$data,s1$f)
  plot_multigrp(s1n,s1$f)
}
\keyword{misc}
