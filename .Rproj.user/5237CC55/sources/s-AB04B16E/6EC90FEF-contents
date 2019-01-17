#' @title Batch pipeline of calling CNV.
#' @description A pipeline for using BAM files or baseline coverage files to
#'   call CNV by the Kalman Filter.
#' @param testBamFile A character string of test BAM file path.
#' @param controlBamFile A character string or vector of control BAM files path.
#' @param baselineFile A character string of baseline coverage file path.
#' @param mode A character string, \emph{control}, \emph{shuffle} or \emph{movingAverage}.
#' @param bedFile A character string of BED file path.
#' @param binSize A non-negative integer giving the width of each bin or window,
#'   such as \strong{1E2}, \strong{1e2} or \strong{100}.
#' @param outDir A character string of directory path to write to.
#' @param annote.database A character string of \emph{refGene.txt} file path. It
#'   can be download from \href{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database}{UCSC database ftp}.
#'   And \emph{refGene.txt} from \emph{Annovar} database is supported, too.
#' @param gapwidth A non-negative integer indicates the maximum gap width
#'   between CNV region to ignore.
#' @param threshold A non-negative numeric value. A CNV of which log2 ratio is
#'   out of range from inverse \emph{threshold} to \emph{threshold} will not
#'   export.
#' @param min.probes A non-negative integer. A CNV of which probe is under
#'   \emph{min.probes} will not export.
#' @param thread A non-negative integer providing the number of threads.
#' @param is.plot Logical. Whether to plot results.
#' @param is.run.by.gender Logial. Whether to call CNV based on sample's gender.
#' @param shuffle.ifByChr A logical object used in mode \emph{shuffle}; if
#'   \code{TRUE}, it will shuffle depth of each chromosome. If \code{FALSE} of
#'   the whole sample.
#' @param moveAverage.step A non-negative integer used in mode
#'   \emph{movingAverage}; it specifies the step by performing moving average.
#' @export
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @author Zhan-Ni Chen
#' @examples
#' ####### Batch pipeline of calling CNV #######
#' testBamFile <- system.file("extdata", 'testSample.bam', package ="kfltCNV")
#' controlBamFile <- c(system.file("extdata", 'controlSample.bam', package ="kfltCNV"),
#'   system.file("extdata", 'controlSample2.bam', package ="kfltCNV"))
#' bedFile <- system.file("extdata", 'chr10_exome.bed', package ="kfltCNV")
#' binSize <- 1E2
#' kfltBatch(testBamFile, controlBamFile = controlBamFile, bedFile = bedFile, binSize = 1E2)
kfltBatch <- function(testBamFile, controlBamFile = NULL, baselineFile = NULL, mode = 'control', bedFile = NULL, binSize = NULL, outDir = NULL, annote.database = NULL, thread = 1, gapwidth = 0, threshold = 0, min.probes = 1, is.plot = TRUE, is.run.by.gender = FALSE, shuffle.ifByChr = FALSE, moveAverage.step = 10 ) {

    if (is.null(outDir)) outDir <- normalizePath('.')
    if (! dir.exists(outDir)) dir.create(outDir)

    # mode: control
    if (mode %in% 'control') {
        if (is.null(baselineFile) & is.null(controlBamFile) ) stop("Mode 'control': baselineFile and controlFile must not be null simultaneously.")
        if (! is.null(baselineFile)) {
            if (! file.exists(baselineFile)) stop(paste0(baselineFile, ' not exsits.'))
        }
        performCreateCovFile(bamFiles = testBamFile, bedFile = bedFile, outDir = outDir, width = binSize, thread = thread)
        testID <- getID(testBamFile)
        testFile <- list.files(path = outDir, pattern = paste0(testID, ".cov$"), all.files = TRUE, full.names = TRUE, recursive = FALSE)
        controlFile <- NULL
        if (! is.null(controlBamFile)) {

            if (length(controlBamFile) == 1 ) {
                if (grepl(',', controlBamFile)){
                    controlBamFile <- unlist(strsplit(controlBamFile, ','))
                }
            }

            # if call cnv by gender
            if (is.run.by.gender) {
                testGender <- callGenderByBam(testBamFile)
                controlGender <- sapply(controlBamFile, callGenderByBam)
                idx <- which(controlGender %in% testGender)
                if (length(idx) > 0) {
                    controlBamFile <- controlBamFile[idx]
                    write(paste0('Use the following control samples to Call CNVs:\n', paste(controlBamFile, collapse = '\n')), stdout())
                } else {
                    write(paste0('Test and control sample share different gender. \nStill use all control samples to Call CNVs.', paste(controlBamFile, collapse = '\n') ), stdout())
                }
            }
            performCreateCovFile(bamFiles = controlBamFile, bedFile = bedFile, outDir = outDir, width = binSize, thread = thread)
            controlID <- sapply(controlBamFile, getID)
            controlFile <- lapply(controlID, function(id) {list.files(path=outDir, pattern=paste0(id, ".cov$"), all.files=TRUE, full.names=TRUE, recursive=FALSE)})
            controlFile <- unlist(controlFile)

        }
        performFitCovFile(testFile = testFile, path = paste0(outDir, '/', testID, '.fit'), baselineFile = baselineFile, controlFile = controlFile)

    } else if (mode %in% 'shuffle' | mode %in% 'movingAverage') {

        # mode: shuffle | movingAverage
        controlFile <- NULL
        performCreateCovFile(bamFiles = testBamFile, bedFile = bedFile, outDir = outDir, width = binSize, thread = thread)
        testID <- getID(testBamFile)
        testFile <- list.files(path = outDir, pattern = paste0(testID, ".cov$"), all.files = TRUE, full.names = TRUE, recursive = FALSE)
        createBaseline(x = paste0(outDir, '/', testID, ".cov"), path = paste0(outDir, "/baseline.cov"), mode = mode, ifByChr = shuffle.ifByChr, step = moveAverage.step)
        performFitCovFile(testFile = testFile, path = paste0(outDir, '/', testID, '.fit'), baselineFile = paste0(outDir, "/baseline.cov"))

    } else {
        stop("Mode must be 'control', 'shuffle' or 'movingAverage'.")
    }

    performRunKflt(file = paste0(outDir, '/', testID, '.fit'), binSize = binSize, outPrefix = paste0(outDir, '/', testID), thread = thread, tmpDir = outDir)
    performCallCNV(file = paste0(outDir, '/', testID, '.state'), outPrefix = paste0(outDir, '/', testID), id = testID, gapwidth = gapwidth, annote.database = annote.database, threshold = threshold, min.probes = min.probes)

    # plot
    if (is.plot) {
        pdf(paste0(outDir, '/', testID, "_state.pdf"), width = 8, height = 3)
        plotKfltResult(file = paste0(outDir, '/', testID, '.state'))
        dev.off()
    }
}

getID <- function(file, ptn=NULL) {
    sampleid  <- rev(unlist(strsplit(file, "/")))[1]
    if(is.null(ptn)) {
        sampleid  <- unlist(strsplit(sampleid, '\\.'))[1]
    }else{
        sampleid  <- unlist(strsplit(sampleid, ptn))[1]
    }
    sampleid
}
