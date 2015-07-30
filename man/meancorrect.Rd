\name{meancorrect}
\alias{meancorrect}
\title{mean-correct a data matrix}
\description{
   Given a numeric matrix, subtracts
the grand means to set the grand mean for each variable to zero
}
\usage{
meancorrect(x)
}
\arguments{
  \item{x}{a numeric matrix or all-numeric data frame}
}
\value{
  A matrix with the grand means of the columns removed
}
\author{Ben Bolker}
\examples{
  X = matrix(1:30,ncol=3)
  meancorrect(X)
}
\keyword{array}
