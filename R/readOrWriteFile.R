#' @title Read BED file and return a GenomicRanges object.
#' @description Read BED file and return a GenomicRanges object. Three required
#'   BED fields are chromosome name, start and end position.
#' @param file A character string of the BED file path.
#' @param genomeSeqinfo A GenomeInfoDb object which contain the sequence
#'   information. It can also be extracted from a BAM file by \pkg{Rsamtools}
#'   package \code{seqinfo(BamFile("bam_file_path"))}.
#' @return A GenomicRanges object of input BED file.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @author Zhan-Ni Chen
#' @examples
#' ####### Read BED file and return a GenomicRanges object. #######
#' require(GenomicRanges)
#' require(Rsamtools)
#' gr <- GRanges(Rle(c("1", "1", "2")),
#'   IRanges(start = c(1 ,10000, 100), end = c(1000, 100000, 10000)))
#' grangetoBed(gr, 'target.bed')
#' bedtoGRange('target.bed')
#'
#' # genome seqinfo
#' bamSeqinfo <- seqinfo(
#'   BamFile(system.file("extdata", 'testSample.bam', package = "kfltCNV")))
#' bedtoGRange('target.bed', genomeSeqinfo = bamSeqinfo)
bedtoGRange <- function(file, genomeSeqinfo = NULL) {
    indat <- read.table(file, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
    fill = TRUE, stringsAsFactors = FALSE)
    chr <- as.character(as.vector(indat[,1]))
    chr <- gsub('chr', '', chr)
    idx <- which(chr %in% c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    if (! is.null(genomeSeqinfo)) {
        if(grepl('chr', seqnames(genomeSeqinfo)[1])) chr <- paste0('chr', chr)
    }
    GRanges(Rle(chr[idx]), IRanges(start = as.numeric(as.vector(indat[idx, 2])),
        end = as.numeric(as.vector(indat[idx, 3]))))
}

#' @title Write a GenomicRanges object into a BED file.
#' @description Write a GenomicRanges object into a BED file with three required
#'   BED fields, i.e. chromosome name, start and end position.
#' @param gr A GenomicRanges object.
#' @param path A character string of the BED file path to write to.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @author Zhan-Ni Chen
#' @examples
#' ####### Write a GenomicRanges object into a BED file #######
#' require(GenomicRanges)
#' gr <- GRanges(Rle(c("1", "1", "2")),
#'   IRanges(start = c(1 ,10000, 100), end = c(1000, 100000, 10000)))
#' grangetoBed(gr, 'target.bed')
grangetoBed <- function(gr, path) {
    gr <- sort(gr)
    cat(paste(as.character(as.vector(seqnames(gr))), start(gr), end(gr), sep = '\t'),
        sep = '\n', file = path)
}

#' @title Read a coverage file.
#' @description Read a coverage file (<filename>.cov) and change it to a sorted
#'   GenomicRanges object. Four required a coverage file fields are chromosome
#'   name, start, end and coverage of depth.
#' @param file A character string of coverage file path.
#' @return A GenomicRanges object of genome region with a column named depth.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @author Zhan-Ni Chen
#' @examples
#' ####### Read a coverage file #######
#' readCovFile(system.file("extdata", 'testSample.cov', package = "kfltCNV"))
readCovFile <- function(file) {
    indat <- read.table(file, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
        fill = TRUE, stringsAsFactors = FALSE)
    sort( GRanges(Rle(as.character(as.vector(indat[,1]))),
        IRanges(start = as.numeric(as.vector(indat[,2])), end = as.numeric(as.vector(indat[,3]))),
        depth = as.numeric(as.vector(indat$V4))) )
}

#' @title Write a GenomicRanges object into a coverage file.
#' @description Write a GenomicRanges object with a column named depth into a
#'   coverage file (<filename>.cov). The result file does not have a header
#'   line, and its four columns respectively indicate chromosome name, start
#'   position, end position, and depth.
#' @param gr A GenomicRanges object with a column named depth.
#' @param path A character string of a coverage file path to write to.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Write a GenomicRanges object into a coverage file #######
#' gr <- readCovFile(system.file("extdata", 'testSample.cov', package = "kfltCNV"))
#' writeCovFile(gr, 'sample.cov')
writeCovFile <- function(gr, path) {
    df <- as.data.frame(gr)
    cat(paste(as.character(as.vector(df[, "seqnames"])),
        as.numeric(as.vector(df[, "start"])), as.numeric(as.vector(df[, "end"])),
        as.numeric(as.vector(df[, "depth"])), sep = '\t'), sep = '\n', file = path)
}
