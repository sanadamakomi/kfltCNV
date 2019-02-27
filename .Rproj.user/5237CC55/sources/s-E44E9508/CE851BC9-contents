#' @title Test two BED files for equal.
#' @description Read two BED files and change them to sorted GenomicRanges
#'   objects, then use \code{identical()} to test for exactly equal. It returns
#'   \code{TRUE} in this case, \code{FALSE} in every other case. Three required
#'   BED fields are chromosome name, start and end position.
#' @param x,y Character strings, a BED file or a file with first three columns
#'   the same format as a BED file, such as a coverage file (<filename>.cov).
#' @return It returns \code{TRUE} if two BED file contain the same regions,
#'   \code{FALSE} in every other case.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Test two BED files for equal #######
#' isSameBedFile(x = system.file("extdata", 'chr10_exome.bed', package = "kfltCNV"),
#' y = system.file("extdata", 'testSample.cov', package = "kfltCNV"))
isSameBedFile <- function(x, y) {
    grX <- sort(bedtoGRange(x))
    grY <- sort(bedtoGRange(y))
    identical(grX, grY)
}

#' @title Align a coverage file.
#' @description Align a coverage file (<filename>.cov) by a BED-like file and
#'   fit \code{NA} data point by linear interpolation.
#' @details Four required fields in a coverage file are chromosome name, start
#'   and end position, and depth of coverage. Three required BED-like fields are
#'   chromosome name, start and end position.
#' @param covFile A character string of a coverage file path.
#' @param axisFile A character string of a BED or BED-like file with first three
#'   columns of chromosome name, start and end position.
#' @return A GenomicRanges object of which region is consistent with
#'   \emph{axisFile} and has an updated depth column.
#' @export
#' @importFrom FKF fkf
#' @importFrom xts rbind.xts
#' @importFrom xts xts
#' @importFrom zoo index
#' @importFrom zoo na.locf
#' @importFrom zoo as.Date
#' @importFrom zoo na.approx
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @author Zhan-Ni Chen
#' @examples
#' ####### Align a coverage file #######
#' require(GenomicRanges)
#' gr <- readCovFile(system.file("extdata", 'testSample.cov', package = "kfltCNV"))
#' writeCovFile(shift(gr, shift = 10L), 'shift.cov')
#'
#' linearInterpolationCov(covFile = 'shift.cov',
#' axisFile = system.file("extdata", 'testSample.cov', package = "kfltCNV"))
linearInterpolationCov <- function(covFile, axisFile) {
    covGR <- readCovFile(covFile)
    axisGR <- bedtoGRange(axisFile)
    covGRlst <- split(covGR, seqnames(covGR))
    axisGRlst <- split(axisGR, seqnames(axisGR))
    newCovGRlst <- axisGRlst
    allchr <- names(newCovGRlst)[which(length(newCovGRlst) > 0)]
    # do each chromosome in covFile
    do <- sapply(allchr, function(chr, ...) {
        gr1 <- axisGRlst[[chr]]
        gr2 <- covGRlst[[chr]]
        # compute cumsum base along position
        orgcov <- width(gr2) * gr2$depth
        cumcov <- cumsum(ifelse(is.na(orgcov), 0, orgcov)) + orgcov * 0
        x1 <- xts(cumcov, as.Date(end(gr2)))
        x2 <- xts(rep(NA, length(end(gr1))), as.Date(end(gr1)))
        x3 <- rbind.xts(x1, x2)
        # use linear interpolation of time series
        x3 <- na.approx(x3, na.rm = FALSE)
        x3 <- na.locf(x3)
        x3 <- na.locf(x3, fromLast = TRUE)
        x4 <- cbind(V1 = as.matrix(x3), V2 = index(x3))
        newcumcov <- sapply(1 : length(end(gr1)), function(k, ...) {
            gatch <- which(x4[,2] %in% end(gr1)[k])
            if(length(gatch) > 0) return(x4[gatch[1],1])
        })

        newcov <- c(newcumcov[1], diff(newcumcov, lag = 1))
        newCovGRlst[[chr]]$depth <- newcov/width(gr2)
    })

    unlist(newCovGRlst)
}

#' @title Merge coverage files.
#' @description Align and merge coverage files (<filename>.cov) with chromosome,
#'   start and end position. Four required fields in a coverage file are
#'   chromosome name, start and end position, and depth of coverage.
#' @param x A character vector contains several coverage files path.
#' @return A data frame, of which columns are chromosome, start position, end
#'   position and depths in input coverage files.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Merge coverage files #######
#' mergerCovFiles(c(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample.cov', package = "kfltCNV")))
mergerCovFiles <- function(x) {
    coverageMatrix <- NULL
    for (i in x) {
        id <- rev(unlist(strsplit(i, "/")))[1]
        id <- unlist(strsplit(id, ".cov"))[1]
        indat <- read.table(i, header = FALSE, sep = "\t", quote = "", comment.char = "#",
            na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
        indat <- indat[, 1:4]
        names(indat)[1:4] <- c("chr", "start", "end", id)
        if (is.null(coverageMatrix)) {
            coverageMatrix <- indat
        } else {
            coverageMatrix <- merge(coverageMatrix, indat, by = c("chr", "start", "end"))
        }
    }
    coverageMatrix
}

#' @title Normalize a coverage matrix.
#' @description Normalize every sample's coverage to the same total bases, i.e.
#'   the minimum total coverage of input samples. Input the result of function
#'   \code{\link{mergerCovFiles}}, which is a data frame with first three
#'   columns of chromosome, start, end, and following by columns of several
#'   samples' depths, and each row means different regions.
#' @param x A data frame or a matrix, the result of function
#'   \code{mergerCovFiles()}.
#' @param total.base A non-negative integer, indicates a total base to scale to.
#' @return A data frame with the same format as input.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Normalize a coverage matrix #######
#' covmatrix <- mergerCovFiles(c(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample.cov', package = "kfltCNV")))
#' normalizeCovMatrix(covmatrix)
normalizeCovMatrix <- function(x, total.base = NULL) {
    w <- as.numeric(as.vector(x[,3])) - as.numeric(as.vector(x[,2])) + 1
    allTotalBase <- apply(x[, 4:ncol(x), drop = FALSE], 2, function(k, ...) { sum(as.numeric(k) * w) }, na.rm = TRUE)
    if (is.null(total.base)) {
        coef <- min(allTotalBase) / allTotalBase
    } else {
        coef <- total.base / allTotalBase
    }
    for (i in colnames(x)[4:ncol(x)]) {
        x[, i] <- x[, i] * coef[i]
    }
    x
}

#' @title Scale a coverage file to a provided total base.
#' @description Scale a coverage file (<filename>.cov) by a provided total base.
#'   Four required fields in a coverage file are chromosome name, start and end
#'   position, and depth of coverage.
#' @param x A character string of a coverage file path.
#' @param path A character string of a file to write to.
#' @param total.base A non-negative integer, indicates a total base to scale to.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Scale a coverage file to a provided total base #######
#' scaleCovFile(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   path = 'testSample.nor.cov', total.base = 1E6)
scaleCovFile <- function(x, path, total.base) {
    indat <- read.table(x, header = FALSE, sep = "\t", quote = "", comment.char = "#",
            na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
    w <- as.numeric(as.vector(indat[,3])) - as.numeric(as.vector(indat[,2])) + 1
    allTotalBase <- sum(as.numeric(as.vector(indat[, 4])) * w, na.rm = TRUE)
    indat[, 4] <- as.numeric(as.vector(indat[, 4])) * total.base / allTotalBase
    cat(paste(as.character(as.vector(indat[, 1])),
        as.numeric(as.vector(indat[, 2])), as.numeric(as.vector(indat[, 3])),
        as.numeric(as.vector(indat[, 4])), sep = '\t'), sep = '\n', file = path)
}

#' @title Compute minimum total bases of several coverage files.
#' @description Calculate minimum total bases of several coverage files
#'   (<filename>.cov). Four required fields in a coverage file are chromosome
#'   name, start, end, and depth of coverage.
#' @param x A data frame or a matrix, the result of function 'mergerCovFiles'.
#' @param method Default \emph{min}, also can input \emph{max}, \emph{median},
#'   or \emph{mean}.
#' @param use.method \code{TRUE} by default, if \code{FALSE}, return a numeric
#'   vector of total bases.
#' @return A numeric value of total bases.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Compute minimum total bases of several coverage files #######
#' computeMinTotalBase(c(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'  system.file("extdata", 'controlSample.cov', package = "kfltCNV")))#'
#'
#' # output total base
#' computeMinTotalBase(c(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'  system.file("extdata", 'controlSample.cov', package = "kfltCNV")), use.method = FALSE)
computeMinTotalBase <- function(x, method = "min", use.method = TRUE) {
    totalbase <- sapply(x, function(i) {
        indat <- read.table(i, header = FALSE, sep = "\t", quote = "", comment.char = "#",
            na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
        w <- as.numeric(as.vector(indat[,3])) - as.numeric(as.vector(indat[,2])) + 1
        sum(as.numeric(as.vector(indat[,4])) * w, na.rm = TRUE)
    })
    if (use.method) {
        return(eval(parse(text = paste0(method, "(c(", paste(totalbase, collapse = ","), "), na.rm = TRUE)"))))
    } else {
        return(totalbase)
    }
}
