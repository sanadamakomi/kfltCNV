#' @title Call gender from BAM file.
#' @description Call gender by read depth of chromosome X and Y. Compare the
#'   chi-squared values obtained to infer whether the male or female assumption
#'   fits read depth better.
#' @param bam A character string of BAM file path.
#' @param lower A non-negative integer. Position of which coverage is lower than
#'   the integer will not be count.
#' @return A character string, \emph{Unknow}, \emph{Female} or \emph{Male}.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import GenomeInfoDb
#' @import S4Vectors
#' @importFrom stats chisq.test
#' @importFrom Rsamtools idxstatsBam
#' @importFrom Rsamtools ScanBamParam
#' @author Zhan-Ni Chen
#' @examples
#' ####### Call gender from BAM file #######
#' callGenderByBam(system.file("extdata", 'testSample.bam', package ="kfltCNV"))
callGenderByBam <- function(bam, lower = 20) {
    df <- idxstatsBam(bam)
    df <- df[1:24,]
    auto_chr <- df[which(df[,'mapped'] == max(df[,'mapped'])), 'seqnames']
    auto_chr <- auto_chr[1]
    bamSeqinfo <- seqinfo(BamFile(bam))
    seqn <- seqnames(bamSeqinfo)
    sexn <- seqn[sapply(c('X', 'Y'), grep, seqn)]
    which <- GRanges(Rle(c(auto_chr,sexn)), ranges = IRanges(start = c(1, 1, 1), end = seqlengths(bamSeqinfo)[c(auto_chr, sexn)]), seqinfo = bamSeqinfo)
    param <- ScanBamParam(which = which)
    aln <- readGAlignments(file = bam, index = bam, param = param)
    cov <- coverage(bam, param = param)
    df <- sapply(cov[c(auto_chr, sexn)], function(c, lower) {
        v <- Views(c, c >= lower)
        return(c(sum(width(v)), sum(sum(v))))
    }, lower = lower )
    ratio <- df[2,] / df[1,]
    ratio_a <- ratio[grep(auto_chr, names(ratio))]
    ratio_x <- ratio[grep('X', names(ratio))]
    ratio_y <- ratio[grep('Y', names(ratio))]
    if (is.na(ratio_x) & is.na(ratio_y)) return('Unknow')
    f_a_x <- chisq.test(c(ratio_a, ratio_x), p = c(2/4, 2/4))$p.value
    m_a_y <- chisq.test(c(ratio_a, ratio_y), p = c(2/3, 1/3))$p.value
    if (ratio_x == 0 & ratio_y == 0) return('Unknow')
    if (ratio_x == 0) {
        if (m_a_y > 5E-2) return('Male')
        return('Female')
    } else if (ratio_y == 0) {
        if (f_a_x > 5E-2) return('Female')
        return('Male')
    } else {
        if (f_a_x > m_a_y & f_a_x > 5E-2) return('Female')
        if (f_a_x < m_a_y & m_a_y > 5E-2) return('Male')
        return('Unknow')
    }
}

#' @title Call gender from coverage file.
#' @description Call gender by read depth of chromosome X and Y. Compare the
#'   chi-squared values obtained to infer whether the male or female assumption
#'   fits read depth better.
#' @param x A character string of coverage file (<fileName>.cov) path.
#' @return A character string, \emph{Unknow}, \emph{Female} or \emph{Male}.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom stats chisq.test
#' @author Zhan-Ni Chen
#' @examples
#' ####### Call gender from coverage file #######
#' callGenderByCov(system.file("extdata", 'testSample.cov', package = "kfltCNV"))
callGenderByCov <- function (x) {
    gr <- readCovFile(x)
    seqn <- as.character(as.vector(runValue(seqnames(gr))))
    sexn <- lapply(c('X', 'Y'), function(x) {seqn[grep(x, seqn)]})
    sexn <- unlist(sexn)
    if (length(sexn) == 0) return('Unknow')
    sexgr <- split(gr, seqnames(gr))
    sexgr <- sexgr[sexn]
    ratio <- sapply(sexgr, function(x) {mean(x$depth, na.rm = TRUE)})
    ratio_a <- mean(gr$depth, na.rm = TRUE )
    ratio_x <- 0
    ratio_y <- 0
    if (grepl('X', names(ratio))) ratio_x <- ratio[grep('X', names(ratio))]
    if (grepl('Y', names(ratio))) ratio_y <- ratio[grep('Y', names(ratio))]
    f_a_x <- chisq.test(c(ratio_a, ratio_x), p = c(2/4, 2/4))$p.value
    m_a_y <- chisq.test(c(ratio_a, ratio_y), p = c(2/3, 1/3))$p.value
    if (ratio_x == 0 & ratio_y == 0) return('Unknow')
    if (ratio_x == 0) {
        if (m_a_y > 5E-2) return('Male')
        return('Female')
    } else if (ratio_y == 0) {
        if (f_a_x > 5E-2) return('Female')
        return('Male')
    } else {
        if (f_a_x > m_a_y & f_a_x > 5E-2) return('Female')
        if (f_a_x < m_a_y & m_a_y > 5E-2) return('Male')
        return('Unknow')
    }
}

#' @title Check BAM file.
#' @description It will stop if BAM file is illegal.
#' @param x A character string or vector of BAM File path.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import S4Vectors
#' @importFrom Rsamtools BamFile
#' @author Zhan-Ni Chen
#' @examples
#' ####### Check BAM file #######
#' checkBam(system.file("extdata", 'testSample.bam', package = "kfltCNV"))
checkBam <- function(x) {
    a <- sapply(x, function(bam) {
        if (! file.exists(bam)) stop(paste0(bam, ' file is missing.'))
        bai <- BamFile(bam)$index
        if (is.na(bai)) stop(paste0(bam, '.bai index file is missing.'))
        bam_info <- file.info(bam)
        bai_info <- file.info(bai)
        dt <- difftime(bai_info$mtime, bam_info$mtime, units = 'secs')
        dt <- as.numeric(dt)
        if (dt < 0) stop(paste0(bam, ' index file is older than BAM file.'))
    })
}
