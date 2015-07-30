\name{makeSymmetric}
\alias{makeSymmetric}
\title{make a symmetrix matrix}
\description{
  takes a matrix and copies its lower or upper triangle
  to generate a symmetric matrix
}
\usage{
makeSymmetric(m, upper = TRUE)
}
\arguments{
  \item{m}{a square matrix}
  \item{upper}{(logical) copy upper triangle to lower triangle?
  (Alternatively, copy upper triangle to lower triangle.)}
}
\value{
  a matrix with the same dimensions as \code{m}, but made symmetric
}
\author{Ben Bolker}
\examples{
  m = matrix(1:9,nrow=3)
  m
  makeSymmetric(m)
  makeSymmetric(m,upper=FALSE)
}
\keyword{misc}
