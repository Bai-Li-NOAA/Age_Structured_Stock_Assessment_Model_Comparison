% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_utils.R
\name{logistic}
\alias{logistic}
\title{Compute various logistic patterns}
\usage{
logistic(pattern, x, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL)
}
\arguments{
\item{pattern}{logistic function option.
Pattern 1. simple logistic pattern
Pattern 2. double logisitc pattern}

\item{x}{Vector or scalar of indepenent variable}

\item{a1}{Slope parameter}

\item{b1}{Location parameter}

\item{a2}{Slope parameter for pattern 2}

\item{b2}{Location parameter for pattern 2}
}
\value{
Values of fcn at x: \code{P}.
}
\description{
Compute various logistic patterns
}
\examples{
\dontrun{
logistic(pattern=1, x=c(1:12), a1=3, b1=1.5)
logistic(pattern=2, x=c(1:12), a1=3, b1=1.5, a2=2, b1=10)
}
}
