\name{plot_dat.theor}
\alias{plot_dat.theor}
\title{Plot multigroup data along with theoretical predictions}
\description{
  Plots multigroup data along with the major axis line and
  95\% confidence ellipses
}
\usage{
plot_dat.theor(x, f, vars = 1:2, cols = 1:length(levels(f)), theor, lines = TRUE, ellipses = TRUE, ...)
}
\arguments{
  \item{x}{a numeric matrix or data frame,
    or (if f is missing) a list of a data matrix and a grouping variable}
  \item{f}{factor describing group structure}
  \item{vars}{vector of the (two) variables to plot}
  \item{cols}{colors}
  \item{theor}{list of theoretical values, as output by \code{\link{sim.theor}}}
  \item{lines}{(logical) plot major axis lines?}
  \item{ellipses}{(logical) plot 95\% confidence ellipses?}
  \item{\dots}{other options to \code{\link{plot}}}
}
\value{
  Produces a plot.
}
\author{Ben Bolker}
\examples{
X = simdata()
Th = sim.theor()
plot_dat.theor(X,theor=Th)
}
\keyword{misc}
