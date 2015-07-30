\name{bp.anova}
\alias{bp.anova}
\title{Analysis of variance incorporating back-projection error}
\description{
  Does an one-way analysis of variance of multiple variables,
incorporating a   term for back-projection variance in the error term
}
\usage{
bp.anova(x, f, verbose = FALSE, use="complete.obs",debug=FALSE,
center=TRUE, nsim=10000)
}
\arguments{
  \item{x}{a numeric matrix (or data frame with all numeric values,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a factor describing the group structure of the data}
  \item{verbose}{(logical) return more information about ANOVA?}
  \item{use}{method for missing observations when computing covariances
    (see \code{\link{cov}} for details)}
  \item{debug}{(logical) print debugging information?}
  \item{center}{center groups before back-projecting?}
  \item{nsim}{number of simulations for $p$-value calculation}
}
\details{
  Back-projects the data; runs ANOVA on each variable; then
  recomputes the $p$-value by simulating the distribution.
}
\value{
  If \code{verbose=FALSE}, a list of \code{anova} objects (one for each
  variable) representing a one-way ANOVA for each variable.  Otherwise,
  a list:
  \item{alist}{ANOVA list}
  \item{bp}{back-projected data}
  \item{cpc.evecs}{CPC eigenvectors}
}
\note{It is the USER'S RESPONSIBILITY to check that the groups in the
  data actually share a common principal component!}
\author{Ben Bolker}
\examples{
   X = simdata(npts=200)
   plot_multigrp(X)
   a1 <- bp.anova(X)
   bp.pvals(a1)
   Y = simdata(npts=200,offset2=2)
   plot_multigrp(Y)
   a2 <- bp.anova(Y,verbose=TRUE,nsim=20000)
   bp.pvals(a2$alist)
   plot_multigrp(a2$bp,Y$f)
}
\keyword{misc}
