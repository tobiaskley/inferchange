% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{lambdapath}
\alias{lambdapath}
\title{Find suitable grid of lambda values}
\usage{
lambdapath(lambda_max, n, p, K = 100)
}
\arguments{
\item{lambda_max}{maximum value of lambda values on grid}

\item{n}{value \code{n}}

\item{p}{value \code{p}}

\item{K}{number of lambda values}
}
\value{
vector with lambda values
}
\description{
Generates K lambda values from \code{lambda_max} to 0.01 (if \code{n < p}) or 0.0001
(if \code{n >= p}) in decreasing order that are equaly spaced on the log scale.
This is the grid of lambda values employed as in \code{glmnet}.
In applications lambda_max is chosen as
}
\examples{
lambdapath(10, 100, 200, K = 5)
}
