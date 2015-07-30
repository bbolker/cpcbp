\name{coverfun}
\alias{coverfun}
\title{Calculate error coverage}
\description{
  Given computed p-values, compute the coverage of a statistical test
}
\usage{
coverfun(p, alpha = 0.05, one.tailed = FALSE, na.rm = TRUE)
}
\arguments{
  \item{p}{A vector of p-values}
  \item{alpha}{Type-I error level}
  \item{one.tailed}{(logical) one-tailed test?}
  \item{na.rm}{Remove NAs from p-value vector?}
}
\value{
  The coverage (fraction of p-values between 0 and alpha for 1-tailed,
  or alpha/2 and 1-alpha/2 for 2-tailed)
}
\author{Ben Bolker}
\examples{
x = matrix(rnorm(1000),ncol=10)
p = pnorm(x)
apply(p,2,coverfun)
}
\keyword{misc}

