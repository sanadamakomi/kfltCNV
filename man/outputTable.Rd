% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callCNV.R
\name{outputTable}
\alias{outputTable}
\alias{outputVcf}
\title{Export CNVs result.}
\usage{
outputTable(x, filepath, info.vec, threshold = 0.5, min.probes = 30)

outputVcf(x, filepath, id, info.vec, threshold = 0.5, min.probes = 30,
  addinfo = NULL, addformat = NULL)
}
\arguments{
\item{x}{A data frame which is result of calling CNVs.}

\item{filepath}{A character string of file path to write to.}

\item{info.vec}{A character vector of column names in x to export.}

\item{threshold}{A non-negative numeric value. A CNV of which log2 ratio is
out of range from inverse \emph{threshold} to \emph{threshold} will not
export.}

\item{min.probes}{A non-negative integer. A CNV of which probe is under
\code{min.probes} will not export.}

\item{id}{A character string of the sample id in VCF file.}

\item{addinfo}{A character string or vecter of column names in x to add to
VCF file header and \emph{INFO} column. Default header can be created by
\code{makeVcfHeader()}.}

\item{addformat}{A character string or vecter of column names in x to add to
VCF file header, \emph{FORMAT} and sample column. Default header can be
created by running \code{makeVcfHeader}.}
}
\description{
Write CNVs result to tabular or VCF format file
}
\examples{
loss <- callLossCNV(system.file("extdata", 'testSample.state', package = "kfltCNV"), 0)
gene <- callCNVGene(system.file("extdata", 'testSample.anno.state', package = "kfltCNV"), 0)

# tabular
outputTable(loss, filepath = 'cnv.result.txt',
  info.vec =  c("svtype", "svlen", "end", "log2ratio", "probe", "depth", "baseline"))

# vcf file
outputVcf(loss, filepath = 'cnv.result.vcf', id = 'testSample',
  info.vec =  c("svtype", "svlen", "end", "log2ratio", "probe", "depth", "baseline"))

# vcf file adding gene
outputVcf(gene, filepath = 'cnv.result.vcf', id = 'testSample',
  info.vec =  c("svtype", "end", 'gene', "log2ratio", "probe", "depth", "baseline"),
  addinfo = '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">')
}
\author{
Zhan-Ni Chen
}
