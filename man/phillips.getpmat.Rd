\name{phillips.getpmat}
\alias{phillips.getpmat}
\alias{phillips.getvec}
\alias{phillips.getsubmat}
\alias{phillips.getmat}
\alias{phillips.getcolonval}
\alias{phillips.getsection}
\title{utility functions for reading output from Phillips' CPC program}
\description{
  Read various forms of input from the non-machine-friendly output
  format of Phillips' code for calculating common principal components
}
\usage{
phillips.getvec(con, n=3, hdrlines = 1)
phillips.getmat(con, n=3, hdrlines = 1)
phillips.getsubmat(con, n=3)
phillips.getpmat(con, n, hdrlines = 1)
phillips.getcolonval(con)
phillips.getsection(con,hdrlines,ntests,
                                 nvar,
                                 ngrp,
                                 matvals=FALSE,
                                 grpevals=FALSE,
                                 grps=TRUE,
                                 nevecs=1,
verbose=FALSE)
}
\arguments{
  \item{con}{a connection that is open for reading}
  \item{n}{(numeric) number of elements in vector or dimension of square
    matrix to read from file}
  \item{hdrlines}{(numeric) number of header lines to skip}
  \item{nvar}{number of variables}
  \item{ngrp}{number of groups}
  \item{ntests}{number of nested hypothesis tests in a section}
  \item{matvals}{does each matrix also have an associated scalar value
    (e.g. proportionality constants)?}
  \item{grpevals}{are there separate eigenvalues for each group?}
  \item{grps}{are there matrices for each group?}
  \item{nevecs}{how many eigenvectors are there to read?}
  \item{verbose}{print debugging information?}
}
\details{\code{getvec()} reads an n-element vector; \code{getmat()}
  reads an n-by-n matrix, with helper file \code{getsubmat()};
  \code{getpmat()} reads a table representing a series of chi-square
  tests of hypotheses; \code{getconline()} reads a single parameter at
  the end of a line of text, delimited by a colon.
  \code{getsection()} gets a whole section.
  }
\value{
  A numeric vector or matrix or a matrix of chi-square, degrees of
  freedom, and p-values, as appropriate.  \code{getsection}
  a list containing a criterion for equivalence (\code{crit});
  number of parameters (\code{par}); a matrix of hypothesis tests
  (\code{testmat}); a set of eigenvalues \code{evals};
  a set of eigenvector matrices (\code{evecs}); and a set of
  variance-covariance matrices (\code{cov}).
}
\author{Ben Bolker}
\seealso{\code{\link{read.cpc}}}
\keyword{misc}
