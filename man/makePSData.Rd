% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makePSData.R
\name{makePSData}
\alias{makePSData}
\alias{makeData}
\alias{createPSData}
\title{Create a survey data set manually}
\usage{
makePSData(n, count = NULL, type = c("P", "S"), notes = NULL)

makeData(n, count = NULL, type = c("P", "S"), notes = NULL)

createPSData(n, count = NULL, type = c("P", "S"), notes = NULL)
}
\arguments{
\item{n}{Either the number of groups of glass or the size of different groups
of glass, or a \code{vector} of observed groups of glass, or group sizes.
See details for a longer explanation.}

\item{count}{Either the number of people in the survey sample who had
\eqn{n}{n} groups of glass on their clothing, or the number of people who
had a group of glass of size \eqn{n}{n}.}

\item{type}{either \code{"P"} or \code{"S"}}

\item{notes}{a \code{\link[utils]{bibentry}} or a character string which
allows extra information about the data to be stored, such as the source,
or reference. \code{NULL} by default.}
}
\value{
an object of type \code{psData}---see \code{\link{readData}} for more
  details.
}
\description{
Create a survey data set from the command line rather than reading data in
from a file. This function is likely to be only useful where there are a very
small number of group sizes, or sizes of groups of glass.
}
\details{
If \code{count} is \code{NULL}, then it is assumed that \code{n}
  consists of actual observed group sizes or numbers of groups of glass found
  on a survey of N individuals. That is, one could provide \code{n = rep(0:1,
  98, 1)} or \code{n = 0:1, count = c(98, 1)}. The former is more useful when
  performing simulation studies.
}
\examples{
## recreate the data read in the readData example
p1 = makePSData(n = c(0, 1, 2), count = c(98, 1, 1), type = "P")
s1 = makePSData(n = 1:3, count = c(1, 1, 1), type = "S")
p1
s1

}
\seealso{
readData
}
