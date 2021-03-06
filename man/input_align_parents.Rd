% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_out.R
\name{input_align_parents}
\alias{input_align_parents}
\title{Align Input Indices and Parents}
\usage{
input_align_parents(phenotype_data, all_data)
}
\arguments{
\item{phenotype_data}{Data frame with phenotype data on all input individuals}

\item{all_data}{Marker data set as data frame}
}
\description{
Given a phenotype and marker data, assembles a data frame with every input individual
in column 1, and parents 1 and 2 in columns 2 and 3. Erroneous individual IDs/indices are
not included.
}
\examples{
input_align_parents()
}
\keyword{.dat}
\keyword{data}
\keyword{formatting}
\keyword{marker}
\keyword{phenotype}
