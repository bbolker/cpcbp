\name{strip.blanks}
\alias{strip.blanks}
\alias{skip.blank.lines}
\alias{skip.pattern.lines}
\title{string utility functions}
\description{
  strip blank characters from the ends of a character string or skip
  blank lines (or lines with a particular pattern) in a file
}
\usage{
strip.blanks(str)
skip.blank.lines(con)
skip.pattern.lines(con,pattern,blank.ok=TRUE)
}
\arguments{
  \item{str}{a character vector}
  \item{con}{an open connection}
  \item{pattern}{a regular expression pattern}
  \item{blank.ok}{skip blank lines as well as lines containing patterns?}
}
\details{
  \code{skip.blank.lines()} and \code{skip.pattern.lines()}
  read lines from a file until the first
  non-blank line or the first line without a specified pattern
  then pushes that line back
}
\value{
  A character vector with blanks removed (for \code{strip.blanks})
  or the number of lines skipped (\code{skip.blank.lines} or
  \code{skip.pattern.lines})
}
\author{Ben Bolker}
\examples{
strip.blanks(c("   a ","b","   c "))
tfn = tempfile("test")
tfile = file(tfn,open="w")
writeLines(c("","  ","abc"),tfile)
close(tfile)
tfile=file(tfn,open="r")
skip.blank.lines(tfile)
readLines(tfile)
close(tfile)
unlink(tfn)
}
\keyword{character}
