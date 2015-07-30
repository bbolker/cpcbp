\name{bp.fstats}
\alias{bp.fstats}
\alias{bp.pvals}
\title{Extract statistics from ANOVA lists}
\description{
  Extracts F statistics, or p values, from lists of
  univariate ANOVAs done on size-corrected data
}
\usage{
bp.fstats(obj)
bp.pvals(obj)
}
\arguments{
  \item{obj}{a list of ANOVA fits (such as produced by \code{\link{bp.anova}})}
}
\value{
  A vector of F statistics or p values, one for each variable
}
\author{Ben Bolker}
\examples{
s1 = simdata()
b1 = bp.anova(s1$data,s1$f)
bp.fstats(b1)
bp.pvals(b1)
}
\keyword{misc}
