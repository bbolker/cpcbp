\name{bp.means}
\alias{bp.means}
\title{Estimate back-projected means and standard deviations}
\description{
  Estimates the group means (or differences between groups),
  raw standard deviations, and
  back-projection error-corrected standard deviations for a set of data
}
\usage{
bp.means(x, f, center = FALSE)
}
\arguments{
  \item{x}{a numeric data matrix or data frame, or
    (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a grouping factor the same length as x}
  \item{center}{center data set at zero?}
}
\value{
  \item{mean or meandiff}{group means or differences between group means
    for each variable}
  \item{sd.raw}{standard deviations of the means or of the differences
    between the means, without correction for back-projection error}
  \item{sd.corr}{standard deviations corrected for back-projection error}
}
\note{It is the USER'S RESPONSIBILITY to check that the groups in the
  data actually share a common principal component!}
\references{Flury}
\author{Ben Bolker}
\examples{
x = simdata()
bp.means(x)
}
\keyword{misc}
