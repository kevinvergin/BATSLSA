% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/acf_seasonality_check.R
\name{s_matrixacf}
\alias{s_matrixacf}
\title{Checks non-normal distribution of residuals from linear model by examining the autocorrelation
function}
\usage{
s_matrixacf(matrix, boolean, fn)
}
\arguments{
\item{inputParameter1}{matrix of NTU abundances in rows and sample dates in columns \code{inputParameter1}}

\item{inputParameter2}{boolean is a matrix of 1's and 0's from the previous linear model code indicating which NTUs are connected (1's) \code{inputParameter2}}

\item{inputParameter3}{fn is a file name for the output  \code{inputParameter3}}
}
\value{
output A boolean matrix of connections that fail (1’s). Results are in the upper triangle. Need to subtract this matrix from the matrix of connections.
}
\description{
Need a dataset with no gaps in the time series. The autocorrelation function is sensitive to the
order of the variables in the linear model so generates data for both cases. Difference
the data (lag=12 for monthly data) to pre-whiten the data.
}
\examples{
Run diagnostic with this command:  s_matrixacf(matrix=, boolean=,fn=)
}
\keyword{keywords}
