\name{bpmat}
\alias{bpmat}
\alias{bpfun}
\title{Flury back-projection matrix}
\description{
   Constructs and applies the Flury back-projection matrix,
which when multiplied by a data vector removes a particular
directional component
}
\usage{
bpmat(evec)
bpfun(x,f,evec,center=TRUE)
}
\arguments{
 \item{x}{(numeric) A data matrix (observations as rows, variables as
   columns) to back-project, or
   (if f is missing) a list of a data matrix and a grouping variable}
 \item{f}{(factor) A grouping factor to use to calculate CPC1 of x,
   if evec is not specified}
 \item{evec}{(numeric) A direction vector to remove from the data
   by back-projection (typically the first common principal
   component): it should be normalized (sum of squares equals 1)}
 \item{center}{mean-correct data before back-projecting?}
}
\details{
  The back-projection matrix B is I - e*t(e);
  back-projecting a vector x is done by multiplying B * x.
}
\value{
  Either the back-projection matrix or the back-projected data
  matrix
}
\note{mean correction is done so that the grand mean of the
  entire data set is zero.  It doesn't matter for balanced data
  sets, but for unbalanced data sets it should (?) probably be
  done so that the means of the group means are zero instead.}
\references{Flury etc.}
\author{Ben Bolker}
\examples{
  X = matrix(c(1,2,3,4,5,6),ncol=3,byrow=TRUE)
  dir = c(0,1,1)/sqrt(2)
  bpmat(dir)
  bpfun(X,evec=dir)
}
\keyword{array}

