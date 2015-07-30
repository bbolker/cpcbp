\name{pbperr}

\alias{rbperr}
\alias{pbperr}
\alias{qbperr}
\alias{dbperr}
\title{Distribution of F statistics corrected for back-projection error}
\description{
  Computes (by random simulation) the null F statistic (ratio of mean
  squares) for a back-projection ANOVA with a specified level of
  back-projection error
}
\usage{
rbperr(n, df1, df2, bpms, atab, scaled, dvar=FALSE)
pbperr(x, df1, df2, bpms, atab, nsim = 20000,
lower.tail= TRUE, ...)
qbperr(q, df1, df2, bpms, atab, nsim = 20000, ...)
dbperr(x, df1, df2, bpms, atab, nsim = 20000, ...)
}
\arguments{
  \item{x}{An observed value}
  \item{n}{Number of deviates}
  \item{q}{Quantile}
  \item{df1}{Numerator (treatment) degrees of freedom}
  \item{df2}{Denominator (error) degrees of freedom}
  \item{bpms}{Estimated scaled mean square associated with back-projection
    error: if s2 is the back-projection variance, then
    bpms = s2*(N/df1)/ms.err}
  \item{atab}{an \code{anova} object: if specified, \code{df1}, \code{df2},
    and \code{s2.err} will be extracted from it}
  \item{nsim}{Number of simulations}
  \item{lower.tail}{Return lower tail of cumulative distribution?}
  \item{scaled}{interpret bpms as scaled mean square or as
    back-projection variance?  (default is TRUE if atab is provided,
    otherwise FALSE)}
  \item{dvar}{incorporate variance of MSE in the simulation?
    (I'm not yet sure which way is right)}
  \item{\dots}{additional arguments to rbperr}
}
\details{
  The statistic returned is like an F(df1,df2) statistic,
  but the numerator adds an additional term
  c*bp.ms, where c is a chi-squared deviate
  with 1 df.  The p, q, and d functions work by simulating
  values and finding the empirical cumulative density,
  quantile, or density (all three of these could be quite crude).
}
\value{
  Random deviates, estimated quantiles or cumulative densities
  as appropriate.
}
\examples{
  curve(df(x,5,10),from=0,to=10,lwd=2,col=1,n=200,
    xlab="Probability density")
  bvec = c(0,0.5,1,2)
  for (i in seq(along=bvec)) {
     curve(dbperr(x,5,10,bpms=bvec[i]),col=i+1,add=TRUE,n=200)
  }
  legend("topright",c("F(5,10)",paste("bperr (bpms=",bvec,")",sep="")),
     col=1:(length(bvec)+1),lty=1)
}
\author{Ben Bolker}
\note{A "p value" of 0 should best be interpreted as a p-value
  of (<1/nsim).  WARNING: the functions are not completely vectorized.}
\seealso{\code{\link{bp.anova}}}
\keyword{misc}
