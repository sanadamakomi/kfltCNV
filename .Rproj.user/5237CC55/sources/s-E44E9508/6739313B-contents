#' @title Create a baseline to call CNV.
#' @description Creating a baseline by \emph{control}, \emph{shuffle} or
#'   \emph{movingAverage}.
#' @details The \emph{control} mode (by default) will compute median depth of
#'   control samples in each region as the baseline. More than one coverage file
#'   (<filename>.cov) is need to pass and if total.base is null, it will scale
#'   to min total base.
#'
#'   The \emph{shuffle} mode will shuffle the depth of each chromosome or the
#'   whole sample (when \emph{isByChr} set to \code{FALSE}) to create a
#'   baseline. If several files input, the first one will be useing.
#'
#'   The \emph{movingAverage} mode will compute moving average depth of each
#'   chromosome or the whole sample to create a baseline. If several files
#'   input, the first one will be useing.
#' @param x A string value or vector of coverage file path.
#' @param path A string value or vector of coverage file path.
#' @param mode A character string, \emph{control}, \emph{shuffle} or \emph{movingAverage}.
#' @param total.base A non-negative numeric value used in mode \emph{control}
#'   indicates a total base to scale to.
#' @param ifByChr A logical object used in mode \emph{shuffle}; if \code{TRUE},
#'   it will shuffle depth in each chromosome. If \code{FALSE} of the whole
#'   sample.
#' @param step A non-negative integer used in mode \emph{movingAverage}; it
#'   specifies the step by performing moving average.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Create a baseline to call CNV #######
#' # control
#' createBaseline(c(system.file("extdata", 'controlSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample2.cov', package = "kfltCNV")),
#'   path = 'baseline.control.cov')
#'
#' # shuffle
#' createBaseline(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   mode = "shuffle", path = 'baseline.shuffle.cov')
#'
#' # movingAverage
#' createBaseline(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   mode = "movingAverage", path = 'baseline.movingAverage.cov')
createBaseline <- function(x, path, mode = "control", total.base = NULL, step = 10, ifByChr = FALSE) {
    if (! mode %in% c("control", "shuffle", "movingAverage")) stop("Mode must be 'control', 'shuffle' or 'movingAverage'.")
    if (mode %in% "control" & length(x) < 2) stop("Create control baseline need at least 2 samples.")
    if (mode %in% c("shuffle", "movingAverage") & length(x) > 1) warning("Only ", paste0(x[1], " will be used to create baseline."))
    for (file in x) {
        if (! file.exists(file)) stop(paste0(file, " no exists."))
    }

    if (mode == "control") {

        for (file in x[2:length(x)]) {
            if (! isSameBedFile(x[1], file)) stop("*cov files do not share the same region.")
        }
        write(paste0('Start creating baseline\nMode:\t', mode, '\nFiles:\n',
            paste(x, collapse='\n'), '\nScale all sample to total base:\t', total.base), stdout())
        coverageMatrix <- mergerCovFiles(x)
        norcoverageMatrix <- normalizeCovMatrix(coverageMatrix, total.base = total.base)
        baseline <- createControlBaseline(norcoverageMatrix)

    } else if (mode == "shuffle") {

        write(paste0('Start creating baseline\nMode:\t', mode, '\nFiles:\n',
            x[1], '\nIf shuffle depth by chr:\t', ifByChr), stdout())
        indat <- read.table(x[1], header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
        fill = TRUE, stringsAsFactors = FALSE)
        baseline <- createShuffleBaseline(indat, ifByChr = ifByChr)

    } else {

        write(paste0('Start creating baseline\nMode:\t', mode, '\nFiles:\n',
            x[1], '\nStep by moving average:\t', step), stdout())
        indat <- read.table(x[1], header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
        fill = TRUE, stringsAsFactors = FALSE)
        baseline <- createMoveAveBaseline(indat, step = step)

    }
    cat(paste(as.character(as.vector(baseline[, 1])),
        as.numeric(as.vector(baseline[, 2])), as.numeric(as.vector(baseline[, 3])),
        as.numeric(as.vector(baseline[, 4])), sep = '\t'), sep = '\n', file = path)
}

#' @title Create a baseline with control samples to call CNV.
#' @description Compute median depth of control samples in each region as the
#'   baseline for calling CNV.
#' @param x A data frame or a matrix with the same format as the result of
#'   function \code{\link{mergerCovFiles}}, of which first three columns are
#'   chromosome, start and end, and following with several samples' depth. It
#'   returns a four-column data frame which has the same first three columns as
#'   input but the forth column is median depth of samples named \emph{baseline}.
#' @return A data frame with four columns named \emph{chr}, \emph{start},
#'   \emph{end}, \emph{baseline}.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Create baseline to call CNV #######
#' covmatrix <- mergerCovFiles(c(
#'   system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample2.cov', package = "kfltCNV")))
#'
#' normal.coveragematrix <- normalizeCovMatrix(covmatrix)
#' baseline.control <- createControlBaseline(normal.coveragematrix[, -4])
createControlBaseline <- function(x) {
    baseline <- x[, 1:4]
    colnames(baseline)[4] <- "baseline"
    baseline[, 4] <- apply(x[, 4:ncol(x), drop = TRUE], 1, function(k, ...) { median(as.numeric(k), na.rm = TRUE) })
    baseline
}

#' @title Creat a baseline by shuffling.
#' @description Shuffle the depth of each chromosome or the whole sample to
#'   create a baseline.
#' @param x A data frame or a matrix as the same format as the result of
#'   function \code{\link{mergerCovFiles}}, with first three columns of
#'   chromosome, start and end, and forth column is depth to shuffle. It returns
#'   a data frame of which first three columns are the same as input but the
#'   forth one is a shuffled depth named \emph{baseline}.
#' @param ifByChr Logical; if \code{TRUE}, it will shuffle depth of each
#'   chromosome. If \code{FALSE} of the whole sample.
#' @return A data frame with four columns named \emph{chr}, \emph{start},
#'   \emph{end}, \emph{baseline}.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Creat a baseline by shuffling #######
#' tumor.coveragematrix <- read.table(
#'   system.file("extdata", 'testSample.cov', package = "kfltCNV"), header = FALSE)
#' colnames(tumor.coveragematrix) <- c('chr', 'start', 'end', 'testSample')
#' baseline.shuffle <- createShuffleBaseline(tumor.coveragematrix)
createShuffleBaseline <- function(x, ifByChr = FALSE) {
    baseline <- x[, 1:4]
    colnames(baseline)[4] <- "baseline"
    if (ifByChr) {
        chr <- as.character(as.vector(baseline[, 1]))
        uniqueChr <- unique(sort(chr))
        for (i in uniqueChr) {
            idx <- which(chr == i)
            baseline[idx, 4] <- shuffling(baseline[idx, 4])
        }
    } else {
        baseline[, 4] <- shuffling(baseline[, 4])
    }
    baseline
}

#' @title Creat a baseline by moving average.
#' @description Compute moving average depth of each chromosome to create a
#'   baseline depth.
#' @param x a data frame or a matrix as the same format as the result of
#'   function \code{\link{mergerCovFiles}}, with first three columns of
#'   chromosome, start and end position, and forth column is depth to compute
#'   moving average. It returns a data frame of which first three columns are
#'   the same as input but the forth one is a moving average depth named
#'   \emph{baseline}.
#' @param step A non-negative integer; it specifies the step by performing
#'   moving average.
#' @return A data frame with four columns.
#' @export
#' @importFrom stats filter
#' @author Zhan-Ni Chen
#' @examples
#' ####### Creat a baseline by moving average #######
#' tumor.coveragematrix <- read.table(
#'   system.file("extdata", 'testSample.cov', package = "kfltCNV"), header = FALSE)
#' colnames(tumor.coveragematrix) <- c('chr', 'start', 'end', 'testSample')
#' baseline.moveaverage <- createMoveAveBaseline(tumor.coveragematrix, step=10)
createMoveAveBaseline <- function(x, step) {
    baseline <- x[, 1:4]
    colnames(baseline)[4] <- "baseline"
    chr <- as.character(as.vector(baseline[, 1]))
    uniqueChr <- unique(sort(chr))
    for (i in uniqueChr) {
        idx <- which(chr == i)
        smoothCov <- filter(baseline[idx, 4], filter = rep(1 / step, step), method = "convolution", sides = 2, circular = TRUE)
        baseline[idx, 4] <- smoothCov
    }
    baseline
}

#' @title Compute the log2 ratio between sample depth and baseline.
#' @description Input vectors of sample depth and baseline to compute the log2
#'   ratio.
#' @param x A numeric vector of sample depth.
#' @param baseline A numeric vector contains baseline depth of region with the
#'   same order in sample depth vector.
#' @param badDepth A numeric value; a baseline depth low than badDepth will be
#'   set to \code{NA}.
#' @return A numeric vector of log2 ratio.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Compute log2ratio #######
#' covmatrix <- mergerCovFiles(c(
#'   system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample2.cov', package = "kfltCNV")))
#'
#' normal.coveragematrix <- normalizeCovMatrix(covmatrix)
#' baseline.control <- createControlBaseline(normal.coveragematrix[, -4])
#' log2ratio <- calculateLog2ratio(normal.coveragematrix[, 4],
#'   baseline.control[, 4], badDepth = 10)
calculateLog2ratio <- function(x, baseline, badDepth) {
    x <- as.numeric(as.vector(x))
    baseline <- as.numeric(as.vector(baseline))
    x[which(x == 0)] <- 0.01
    baseline[which(baseline < badDepth)] <- NA
    log2ratio <- log2(x / baseline)
    log2ratio[which(is.na(log2ratio))] <- NA
    log2ratio
}

shuffling <- function(x) {
    set.seed(123)
    x <- as.numeric(x)
    theta.boot <- replicate(1000, expr = {
        sample(x, size = length(x), replace = FALSE)
    })
    apply(theta.boot, 1, median, na.rm = TRUE)
}
