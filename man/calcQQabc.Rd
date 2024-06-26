% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcQQabc.R
\name{calcQQabc}
\alias{calcQQabc}
\title{Compute the area between curves (abc) in the Q-Q plot}
\usage{
calcQQabc(P, Q, sample1, sample2)
}
\arguments{
\item{P}{A numeric vector of sampled data points from the first distribution.}

\item{Q}{A numeric vector of sampled data points from the second distribution.}

\item{sample1}{A character to denote the first distribution.}

\item{sample2}{A character to denote the second distribution.}
}
\value{
calcQQabc returns a vector of the area between curves (abc) values for each gene using the Q-Q method.
}
\description{
This function computes the area between curves in the Q-Q plot.
             It takes the empirical and theoretical quantiles as input and calculates the area between the 45 degree line and the the quantile of the theoretical distribution against the quantile of observed distribution.
}
