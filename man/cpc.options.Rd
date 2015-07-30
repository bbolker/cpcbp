\name{cpc.options}
\alias{cpc.options}
\title{Set CPC calculation options}
\description{
  Sets or returns global options for computing common principal components
  (currently, just the function to use)
}
\usage{
cpc.options(cpcvecfun, neg.bpvar.tol)
}
\arguments{
  \item{cpcvecfun}{Program to compute CPC}
  \item{neg.bpvar.tol}{tolerance for warnings about negative BP variances}
}
\details{
   \code{\link{phillips.cpc}} uses a call to
    a function written by Phillips; \code{\link{pooled.cpc}} assumes common
    principal components, calculating them by mean-correcting each group
    and computing principal components
}
\value{
  \item{cpcvecfun}{CPC computation function}
  \item{neg.bpvar.tol}{tolerance for warnings about negative BP variances}
}
\author{Ben Bolker}
\keyword{misc}
