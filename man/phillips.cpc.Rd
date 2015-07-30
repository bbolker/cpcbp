\name{phillips.cpc}
\alias{phillips.cpc}
\alias{phillips.cpcvec}
\alias{run.cpc}
\alias{read.cpc}
\title{Run Phillips's CPC program from R}
\description{
  Provides an interface to P. C. Phillips's implementation of Flury's
  algorithm for computing common principal components
}
\usage{
phillips.cpc(x, f, covs, npts,
                     progdir,
                     progname,
                     ansfn,
                     datfn,
                     outfn,
                     unlink.temp=TRUE,
                     use="complete.obs",verbose=FALSE)
phillips.cpcvec(x, f, covs, npts,
                     progdir,
                     progname,
                     ansfn,
                     datfn,
                     outfn,
                     unlink.temp=TRUE,
                     use="complete.obs")
read.cpc(outfn,ngrp,nvar,verbose=FALSE)
run.cpc(covs,npts,progdir,
                    progname,
                    ansfn,
                    datfn,
                    outfn)
}
\arguments{
  \item{x}{a numeric matrix (or data frame with all numeric values,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{a factor describing the group structure of the data}
  \item{covs}{list of pre-computed covariance matrices of all groups}
  \item{npts}{numeric vector, number of observations in each group}
  \item{ansfn}{file name for CPC control file}
  \item{datfn}{file name for CPC temporary data file}
  \item{outfn}{file name for CPC output file}
  \item{unlink.temp}{(logical) delete working files when finished?}
  \item{progdir}{(character) directory containing Phillips's CPC program}
  \item{progname}{(character) name of Phillips's CPC program}
  \item{ngrp}{(numeric) number of groups}
  \item{nvar}{(numeric) number of variables}
  \item{use}{method for missing observations when computing covariances
    (see \code{\link{cov}} for details)}
  \item{verbose}{print debugging information?}
}
\details{
  \code{phillips.cpc()} is a wrapper for \code{run.cpc()}, which writes
  appropriate input files to disk and then calls Phillips's CPC program,
  and \code{read.cpc()}, which reads the resulting output file back into
  R.  At the moment \code{read.cpc()} and \code{phillips.cpc()} only return
  the CPC eigenvector matrix, although they could be extended.
  \code{phillips.cpcvec()} is a wrapper for \code{phillips.cpc()}
  which returns only the common eigenvector matrix (rather than a
  longer list of results).
}
\value{
  \code{run.cpc()} returns nothing (it calls the CPC program which
  creates an output file). 
  \code{read.cpc()} returns a long list (see below).
  \code{phillips.cpcvec()} returns just the eigenvector matrix.
}
\note{
  The list returned by \code{read.cpc} consists of:
  \describe{
    \item{evecs.CPC}{the CPC eigenvector matrix returned by
      the CPC program}
    \item{cpc1.pval}{the p-value for the test of 1 CPC vs unrelated
      matrices}
    \item{datlist}{a list of lists for sub-hypotheses, each
      containing crit (criterion), par (number of parameters for
      the hypothesis), testmat (statistical test matrix),
      evals (list of eigenvalues for each group),
      evecs (eigenvector matrix or matrices),
      cov (list of variance-covariance matrices for each group)}
    \item{pool}{a list for the pooled hypothesis}
    \item{propdat}{a list for the proportional hypothesis}
  }
}
\references{Uses CPC V0.96  (c) 1994-1999 Patrick C. Phillips,
  pphillips@uta.edu}
\seealso{\code{\link{pooled.cpc}}}
\author{Ben Bolker}
\examples{
   set.seed(1001)
   v = c(10,5,2,1)
   corr = c(0.8,0.6,0.5,0.5)
   X = simdata(npts=1000,vars=v,cors=corr)
   eigen(covmat(v,cor=corr))
   P1 = phillips.cpc(X)
   P1$cpc1.pval
   P1$evecs.CPC
   phillips.cpcvec(X)
   pooled.cpc(X)
}
\keyword{misc}



