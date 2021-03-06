% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kfltLog2Ratio.R
\name{performRunKflt}
\alias{performRunKflt}
\alias{runKflt}
\alias{runKfltMultiCore}
\title{Perform Kalman filtering to estimate the log2 ratio level.}
\usage{
performRunKflt(file, binSize, outPrefix = NULL, model = list(p = 0, d =
  1, q = 5, bin.size = NULL, sn = 2.5, method = "nlminb"), s = 7,
  alpha = 0.05, thread = 1, tmpDir = NULL)

runKflt(x, model = list(p = 0, d = 1, q = 5, bin.size = 1e+05, sn = 2.5,
  method = "nlminb"), s, alpha)

runKfltMultiCore(x, model = list(p = 0, d = 1, q = 5, bin.size = 1e+05,
  sn = 2.5, method = "nlminb"), s, alpha, thread, tmpDir = NULL)
}
\arguments{
\item{file}{A string value of a log2 ratio file (<filename>.fit) path to do
Kalman filtering.}

\item{binSize}{A non-negative integer giving the width of each bin or window,
such as \strong{1E2}, \strong{1e2} or \strong{100}.}

\item{outPrefix}{A character string, indicates output files' prefix.}

\item{model}{A list object contains elements named \emph{p}, \emph{d},
\emph{q}, \emph{bin.size}, \emph{sn}, and \emph{method}. See details in the
function \code{\link{kflt}}. By default it is for log2 ratio data.}

\item{s}{A non-negative integer indicating the step width to smooth the state
after Kalman filtering.}

\item{alpha}{A numeric value providing the significant level.}

\item{thread}{A non-negative integer providing the number of threads.}

\item{tmpDir}{A character string of the directory path to write to.}

\item{x}{A data frame of a log2 ratio file (<filename>.fit).}
}
\description{
Use the Kalman filter to estimate a log2 ratio level and mark
  regions that significantly different from average level as copy number
  variantions.
}
\details{
The input file is a log2 ratio file (<filename>.fit) which has six
  columns respectively named \emph{chr}, \emph{start}, \emph{end},
  \emph{log2ratio}, \emph{depth} and \emph{baseline}. And \emph{log2ratio} is
  the log2 ratio of test sample's depth to baseline. The columns named
  \emph{depth} and \emph{baseline} are depth of coverage respectively in test
  sample and the baseline.

Results are two files named by prefix.

One is a estimated log2 ratio state file (<filename>.state). It has ten
  columns and the first six columns are from input log2 ratio file
  (<filename>.fit). Then another four columns are estimated log2 ratio state by
  the Kalman filter:

    \emph{state} - the log2 ratio state after kalman filtering.

    \emph{statUp} - the upper limmit of state confidence interval.

    \emph{statDown} - the lower limmit of state confidence interval.

    \emph{report} - the label representing the copy number variation level of
    a region. It has three levels i.e. gain, loss and average.

The other file is a parameter file (<filename>.parameter). It stores a best
  set of parameters found by \code{optim()}.
  \code{performRunKflt()} use the \code{fkf()} function in the package
  \pkg{FKF} to perform Kalman filtering and \code{optim()} to estimate
  parameters. By default it will simulate a \strong{ARMA(0,5)-process}, so a parameter
  file has 6 columns named \emph{par1} to \emph{par6} representing six
  parameters in ARMA process, while the first column is chromosome name because
  it performs Kalman filtering by chromosome.
  See the details in the function \code{\link{kflt}}.

It can perform in multithreading with number of threads you input. If the
  \emph{thread} is over 1, it will perform in multithreading, default by single
  thread.
\code{runKflt()} and \code{runKfltMultiCore()} can run from a data frame and
  return a list of result.
}
\examples{
####### Perform Kalman filtering to estimate the log2 ratio level #######
fitFile <- system.file("extdata", 'testSample.fit', package = "kfltCNV")
performRunKflt(fitFile, 'testSample', binSize = 1E2)

fit.dat <- read.delim(fitFile)
runKflt(fit.dat, s = 7, alpha = 0.05)
}
\author{
Zhan-Ni Chen
}
