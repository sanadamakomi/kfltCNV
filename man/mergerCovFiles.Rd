% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalizeCovMatrix.R
\name{mergerCovFiles}
\alias{mergerCovFiles}
\title{Merge coverage files.}
\usage{
mergerCovFiles(x)
}
\arguments{
\item{x}{A character vector contains several coverage files path.}
}
\value{
A data frame, of which columns are chromosome, start position, end
  position and depths in input coverage files.
}
\description{
Align and merge coverage files (<filename>.cov) with chromosome,
  start and end position. Four required fields in a coverage file are
  chromosome name, start and end position, and depth of coverage.
}
\examples{
####### Merge coverage files #######
mergerCovFiles(c(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
  system.file("extdata", 'controlSample.cov', package = "kfltCNV")))
}
\author{
Zhan-Ni Chen
}
