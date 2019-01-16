#' @title Spliting a BED file into windows.
#' @description Read a BED file to a GenomicRanges object and split it into
#'   windows with provided width. Three required BED fields are chromosome name,
#'   start and end position.
#' @param file A character string of a BED file path.
#' @param width An integer giving the width of each window, such as
#'   \strong{1E2}, \strong{1e2} or \strong{100}.
#' @param genomeSeqinfo A GenomeInfoDb object which contains sequence
#'   information. It can also be extracted from a BAM file by \pkg{Rsamtools}
#'   package \code{seqinfo(BamFile("bam_file_path"))}.
#' @return A GenomicRanges object of splited regions from the BED file input.
#' @export
#' @importFrom S4Vectors width
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom GenomicRanges slidingWindows
#' @importFrom GenomicRanges findOverlaps
#' @author Zhan-Ni Chen
#' @examples
#' ####### Spliting a BED file into windows #######
#' require(GenomicRanges)
#' require(Rsamtools)
#' gr <- GRanges(Rle(c("1", "1", "2")),
#'   IRanges(start = c(1 ,10000, 100), end = c(1000, 100000, 10000)))
#' grangetoBed(gr, 'target.bed')
#' splitBed('target.bed', 1E2)
#'
#' # genome seqinfo
#' bamSeqinfo <- seqinfo(BamFile(system.file("extdata", 'testSample.bam', package = "kfltCNV")))
#' splitBed('target.bed', 1E2, genomeSeqinfo = bamSeqinfo)
splitBed <- function(file, width, genomeSeqinfo = NULL) {
    if (is.null(genomeSeqinfo)) {
        gr <- bedtoGRange(file)
    } else {
        gr <- bedtoGRange(file, genomeSeqinfo)
    }
    grlist <- slidingWindows(gr, width = width, step = width)
    idx0 <- which(lengths(grlist) == 1)
    grlist2 <- unlist(grlist)
    idx <- which(width(grlist2) <= width / 2)
    if (length(idx0) > 0) {
        gridx0 <- grlist[idx0]
        gridx0 <- unlist(gridx0)
        newidx0 <- subjectHits(findOverlaps(gridx0, grlist2, type = "equal"))
        idx <- setdiff(idx, newidx0)
    }
    idx2 <- setdiff(1:length(grlist2), idx)
    idx3 <- setdiff(1:length(grlist2), idx - 1)
    GRanges(Rle(seqnames(grlist2)[idx2]), IRanges(start = start(grlist2)[idx2], end = end(grlist2)[idx3]))
}

#' @title Spliting genome into windows.
#' @description Split genome into windows with provided width.
#' @param genomeSeqinfo A GenomeInfoDb object which contain the sequence
#'   information. It can also be extracted from a BAM file by \pkg{Rsamtools}
#'   package \code{seqinfo(BamFile("bam_file_path"))}.
#' @param width An integer giving the width of each window, such as
#'   \strong{1E2}, \strong{1e2} or \strong{100}.
#' @return A GenomicRanges object of splited regions from genome.
#' @export
#' @import S4Vectors
#' @import IRanges
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @author Zhan-Ni Chen
#' @examples
#' ####### Spliting genome into windows #######
#' require(Rsamtools)
#' bamSeqinfo <- seqinfo(BamFile(system.file("extdata", 'testSample.bam', package = "kfltCNV")))
#' splitGenome <- splitGenome(bamSeqinfo, width = 1E7)
splitGenome <- function(genomeSeqinfo, width) {
    chr <- as.character(as.vector(seqnames(genomeSeqinfo)))
    chr <- gsub('chr', '', chr)
    idx <- which(chr %in% c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    gr <- GRanges(Rle(seqnames(genomeSeqinfo)[idx]), IRanges(start = 1, width = seqlengths(genomeSeqinfo)[idx]),
        seqinfo = genomeSeqinfo)
    grlist <- slidingWindows(gr, width = width, step = width)
    unlist(grlist)
}
