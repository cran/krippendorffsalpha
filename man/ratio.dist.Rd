% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/krippendorff.R
\name{ratio.dist}
\alias{ratio.dist}
\title{Apply a ratio distance function to two scores.}
\usage{
ratio.dist(x, y)
}
\arguments{
\item{x}{a score.}

\item{y}{a score.}
}
\value{
\eqn{(x-y)^2/(x+y)^2}, or 0 if \code{x} or \code{y} is \code{NA}.
}
\description{
Apply a ratio distance function to two scores.
}
\details{
This function applies a ratio distance function to two scores. This may be an appropriate distance function for the ratio level of measurement. \code{NA}'s are handled gracefully.
}
\seealso{
\code{\link{interval.dist}}, \code{\link{nominal.dist}}
}
