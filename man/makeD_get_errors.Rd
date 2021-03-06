% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_errors.R
\name{makeD_get_errors}
\alias{makeD_get_errors}
\title{Filter Clerical Errors In Input Set}
\usage{
makeD_get_errors(phenotype_data, all_data)
}
\arguments{
\item{all_data}{Marker data set as data frame}

\item{phen_data}{Data frame with phenotype data on all input individuals}
}
\description{
Cross references input individual set (phenotype data) and
pool of individuals which have on-file marker data.
Produces a 1-column data frame containing all input individuals not found
in the master marker data set.
}
\examples{
makeD_get_errors()
}
\keyword{data}
\keyword{errors}
\keyword{input}
\keyword{marker}
\keyword{phenotype}
\keyword{set}
