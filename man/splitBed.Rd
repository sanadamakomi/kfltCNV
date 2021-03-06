% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splitBed.R
\name{splitBed}
\alias{splitBed}
\title{Spliting a BED file into windows.}
\usage{
splitBed(file, width, genomeSeqinfo = NULL)
}
\arguments{
\item{file}{A character string of a BED file path.}

\item{width}{An integer giving the width of each window, such as
\strong{1E2}, \strong{1e2} or \strong{100}.}

\item{genomeSeqinfo}{A GenomeInfoDb object which contains sequence
information. It can also be extracted from a BAM file by \pkg{Rsamtools}
package \code{seqinfo(BamFile("bam_file_path"))}.}
}
\value{
A GenomicRanges object of splited regions from the BED file input.
}
\description{
Read a BED file to a GenomicRanges object and split it into
  windows with provided width. Three required BED fields are chromosome name,
  start and end position.
}
\examples{
####### Spliting a BED file into windows #######
require(GenomicRanges)
require(Rsamtools)
gr <- GRanges(Rle(c("1", "1", "2")),
  IRanges(start = c(1 ,10000, 100), end = c(1000, 100000, 10000)))
grangetoBed(gr, 'target.bed')
splitBed('target.bed', 1E2)

# genome seqinfo
bamSeqinfo <- seqinfo(BamFile(system.file("extdata", 'testSample.bam', package = "kfltCNV")))
splitBed('target.bed', 1E2, genomeSeqinfo = bamSeqinfo)
}
\author{
Zhan-Ni Chen
}
