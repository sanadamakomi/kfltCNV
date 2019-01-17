#' @title Prepare a log2 ratio file.
#' @description Prepare a log2 ratio file (<filename>.fit) to do Kalman filtering.
#' @details A log2 ratio file has six columns respectively named \emph{chr},
#'   \emph{start}, \emph{end}, \emph{log2ratio}, \emph{depth} and
#'   \emph{baseline}.
#'
#'   And \emph{log2ratio} is the log2 ratio of test sample's depth to the
#'   baseline. The columns named \emph{depth} and \emph{baseline} are depth of
#'   coverage respectively in test sample and the baseline.
#'
#'   If \emph{baselineFile} isn't \code{NULL}, it will use baseline coverage
#'   file to create the log2 ratio file. When you input control samples'
#'   coverage files, it will use control samples to create a baseline at first.
#' @param testFile A character string of a coverage file path to call CNV.
#' @param path A character string of the log2 ratio file to output.
#' @param baselineFile A character string of a baseline coverage file.
#' @param controlFile A character string or vector of control samples' coverage
#'   files.
#' @param badDepth A non-negative numeric value. The region of which depth in
#'   baseline low than badDepth will be set to \code{NA}.
#' @export
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom stats dgamma
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats dbeta
#' @importFrom stats optim
#' @importFrom stats nlminb
#' @importFrom stats filter
#' @author Zhan-Ni Chen
#' @examples
#' ####### Prepare a log2 ratio file #######
#' testFile <- system.file("extdata", 'testSample.cov', package = "kfltCNV")
#' controlFile <- c(system.file("extdata", 'controlSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample2.cov', package = "kfltCNV"))
#'
#' performFitCovFile(testFile, path = 'testSample.fit', controlFile = controlFile)
performFitCovFile <- function(testFile, path, baselineFile = NULL, controlFile = NULL, badDepth = 10) {
    if (is.null(baselineFile) & is.null(controlFile)) stop('BaselineFile and controlFile must not be null simultaneously.')
    allFile <- c(testFile, baselineFile, controlFile)
    allFile <- allFile[which(! sapply(allFile, is.null))]
    for (file in allFile) {
        if (! file.exists(file)) stop(paste0(file, " no exists."))
    }
    for (file in setdiff(allFile, testFile)) {
        if (! isSameBedFile(testFile, file)) stop("Coverage files do not share the same region.")
    }
    write(paste0('Start to create log2 ratio file\ntestFile:\n',
        testFile, '\nbaselineFile:\n', baselineFile, '\ncontrolFile:\n',
        paste(controlFile, collapse = '\n')), stdout())
    if (! is.null(baselineFile) ) {
        coverageMatrix <- mergerCovFiles(c(testFile, baselineFile))
        log2ratio <- calculateLog2ratio(x = coverageMatrix[,4],
            baseline = coverageMatrix[, 5], badDepth = badDepth)
        log2ratioDF <- data.frame(chr = coverageMatrix[,1],
                          start = coverageMatrix[,2],
                          end = coverageMatrix[,3],
                          log2ratio = log2ratio,
                          depth = coverageMatrix[,4],
                          baseline = coverageMatrix[,5])
        write.table(log2ratioDF, file = path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    } else {
        coverageMatrix <- mergerCovFiles(c(testFile, controlFile))
        norcoverageMatrix <- normalizeCovMatrix(coverageMatrix)
        if (length(controlFile) > 1) {
            baseline <- createControlBaseline(norcoverageMatrix[, setdiff(1:ncol(norcoverageMatrix), 4)])
        } else {
            baseline <- norcoverageMatrix[, setdiff(1:ncol(norcoverageMatrix), 4)]
        }
        log2ratio <- calculateLog2ratio(x = norcoverageMatrix[,4],
            baseline = baseline[, 4], badDepth = badDepth)
        log2ratioDF <- data.frame(chr = norcoverageMatrix[,1],
                          start = norcoverageMatrix[,2],
                          end = norcoverageMatrix[,3],
                          log2ratio = log2ratio,
                          depth = coverageMatrix[,4],
                          baseline = baseline[,4])
        write.table(log2ratioDF, file = path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }
}

#' @title Perform Kalman filtering to estimate the log2 ratio level.
#' @description Use the Kalman filter to estimate a log2 ratio level and mark
#'   regions that significantly different from average level as copy number
#'   variantions.
#' @details The input file is a log2 ratio file (<filename>.fit) which has six
#'   columns respectively named \emph{chr}, \emph{start}, \emph{end},
#'   \emph{log2ratio}, \emph{depth} and \emph{baseline}. And \emph{log2ratio} is
#'   the log2 ratio of test sample's depth to baseline. The columns named
#'   \emph{depth} and \emph{baseline} are depth of coverage respectively in test
#'   sample and the baseline.
#'
#' Results are two files named by prefix.
#'
#' One is a estimated log2 ratio state file (<filename>.state). It has ten
#'   columns and the first six columns are from input log2 ratio file
#'   (<filename>.fit). Then another four columns are estimated log2 ratio state by
#'   the Kalman filter:
#'
#'     \emph{state} - the log2 ratio state after kalman filtering.
#'
#'     \emph{statUp} - the upper limmit of state confidence interval.
#'
#'     \emph{statDown} - the lower limmit of state confidence interval.
#'
#'     \emph{report} - the label representing the copy number variation level of
#'     a region. It has three levels i.e. gain, loss and average.
#'
#' The other file is a parameter file (<filename>.parameter). It stores a best
#'   set of parameters found by \code{optim()}.
#'   \code{performRunKflt()} use the \code{fkf()} function in the package
#'   \pkg{FKF} to perform Kalman filtering and \code{optim()} to estimate
#'   parameters. By default it will simulate a \strong{ARMA(0,5)-process}, so a parameter
#'   file has 6 columns named \emph{par1} to \emph{par6} representing six
#'   parameters in ARMA process, while the first column is chromosome name because
#'   it performs Kalman filtering by chromosome.
#'   See the details in the function \code{\link{kflt}}.
#'
#' It can perform in multithreading with number of threads you input. If the
#'   \emph{thread} is over 1, it will perform in multithreading, default by single
#'   thread.
#' \code{runKflt()} and \code{runKfltMultiCore()} can run from a data frame and
#'   return a list of result.
#' @param file A string value of a log2 ratio file (<filename>.fit) path to do
#'   Kalman filtering.
#' @param x A data frame of a log2 ratio file (<filename>.fit).
#' @param outPrefix A character string, indicates output files' prefix.
#' @param binSize A non-negative integer giving the width of each bin or window,
#'   such as \strong{1E2}, \strong{1e2} or \strong{100}.
#' @param model A list object contains elements named \emph{p}, \emph{d},
#'   \emph{q}, \emph{bin.size}, \emph{sn}, and \emph{method}. See details in the
#'   function \code{\link{kflt}}. By default it is for log2 ratio data.
#' @param s A non-negative integer indicating the step width to smooth the state
#'   after Kalman filtering.
#' @param alpha A numeric value providing the significant level.
#' @param thread A non-negative integer providing the number of threads.
#' @param tmpDir A character string of the directory path to write to.
#' @export
#' @import parallel
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom methods as
#' @author Zhan-Ni Chen
#' @examples
#' ####### Perform Kalman filtering to estimate the log2 ratio level #######
#' fitFile <- system.file("extdata", 'testSample.fit', package = "kfltCNV")
#' performRunKflt(fitFile, 'testSample', binSize = 1E2)
#'
#' fit.dat <- read.delim(fitFile)
#' runKflt(fit.dat, s = 7, alpha = 0.05)
performRunKflt <- function(file, binSize, outPrefix = NULL, model = list(p = 0, d = 1, q = 5, bin.size = NULL, sn = 2.5, method = "nlminb"), s = 7, alpha = 0.05, thread = 1, tmpDir = NULL) {
    if (is.null(outPrefix)) {
        id <- rev(unlist(strsplit(file, "/")))[1]
        id <- unlist(strsplit(id, "\\."))[1]
        outPrefix <- paste0(normalizePath('.'), '/', id)
    }
    if (is.null(tmpDir)) tmpDir <- '.'
    tmpDir <- normalizePath(tmpDir)
    model$bin.size <- binSize
    log2ratioDF <- read.table(file, header = TRUE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
    fill = TRUE, stringsAsFactors = FALSE)
    log2ratioDF[, 'chr'] <- as.character(as.vector(log2ratioDF[, 'chr']))
    if (is.null(model$bin.size)) {
        medianBin <- median(log2ratioDF$end - log2ratioDF$start, na.rm = TRUE)
        model$bin.size <- medianBin
    }
    write(paste0('Start Kalman filtering\nModel:\n',
        paste(paste(names(model), "=", model[names(model)]), collapse = ', '),
        '\nFiles:\n', file), stdout())
    if (thread > 1) {
        result <- runKfltMultiCore(log2ratioDF, model = model, s = s, alpha = alpha, thread = thread, tmpDir = tmpDir)
    } else {
        result <- runKflt(log2ratioDF, model = model, s = s, alpha = alpha)
    }
    state.df <- result[[1]]
    state.df <- do.call("rbind", state.df)
    pars <- result[[2]]
    par.df <- sapply(names(pars), function(n) { c(n, pars[[n]])} )
    par.df <- t(par.df)
    par.df <- as.data.frame(par.df)
    colnames(par.df) <- c("chr", paste0('par', seq(1, ncol(par.df)-1, 1)))
    write.table(state.df, file = paste0(outPrefix, '.state'), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    write.table(par.df, file = paste0(outPrefix, '.parameter'), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

#' @rdname performRunKflt
#' @export
runKflt <- function(x, model = list(p = 0, d = 1, q = 5, bin.size = 1E5, sn = 2.5, method = "nlminb"), s, alpha) {
    finalDF <- data.frame(A=x[, "log2ratio"],
                      M=rep(NA, nrow(x)),
                      position=x[, "end"],
                      chr=x[, "chr"],
                      start=x[, "start"],
                      end=x[, "end"],
                      depth=x[, "depth"],
                      baseline=x[, "baseline"])
    finalDF.list <- split(finalDF, as.factor(finalDF$chr))
    state.list <- list()
    par.list <- list()
    for (i in 1:length(finalDF.list)) {
        fitData <- finalDF.list[[i]]
        result <- runKfltPer(fitData, model = model, s = s, alpha =  alpha)
        chr <- as.character(result$chr)
        state.list[[chr]] <- result$state
        par.list[[chr]] <- result$par
    }
    return(list(state.list, par.list))
}

#' @rdname performRunKflt
#' @export
runKfltMultiCore <- function(x, model = list(p = 0, d = 1, q = 5, bin.size = 1E5, sn = 2.5, method = "nlminb"), s, alpha, thread, tmpDir = NULL) {
    timestamp <- paste0('kflt_', as.character(as.integer(Sys.time())))
    if (is.null(tmpDir)) tmpDir <- '.'
    tmp <- normalizePath(tmpDir)
    tmp <- paste0(tmp, '/tmp_', timestamp)
    if (! file.exists(tmp) ) dir.create(path = tmp)
    finalDF <- data.frame(A=x[, "log2ratio"],
                      M=rep(NA, nrow(x)),
                      position=x[, "end"],
                      chr=x[, "chr"],
                      start=x[, "start"],
                      end=x[, "end"],
                      depth=x[, "depth"],
                      baseline=x[, "baseline"])
    finalDF.list <- split(finalDF, as.factor(finalDF$chr))

    clsp <- makeCluster(thread, type = "FORK", outfile = paste0(tmp, "/log.txt"))
    a <- parLapply(clsp, finalDF.list, function(x, ...) {
        result <- runKfltPer(x, model = model, s = s, alpha =  alpha)
        save(result, file = paste0(tmp, "/", timestamp, "_", x[1, "chr"],".Rdata"))
    })
    stopCluster(clsp)
    allrdata <- list.files(path = tmp, pattern = ".Rdata", all.files = TRUE, full.names = TRUE, recursive = FALSE, include.dirs = FALSE)
    state.list <- list()
    par.list <- list()
    for (i in 1:length(allrdata)) {
        result <- get(load(allrdata[i]))
        chr <- as.character(result$chr)
        state.list[[chr]] <- result$state
        par.list[[chr]] <- result$par
    }
    unlink(tmp, recursive = TRUE)
    return(list(state.list, par.list))
}

runKfltPer <- function(fitData, model = list(p = 0, d = 1, q = 5, bin.size = 1E5, sn = 2.5, method = "nlminb"), s, alpha) {
    omitData <- NULL
    if (is.na(fitData[1, "A"])) {
        nonNaIndex <- which(!is.na(fitData[, "A"]))
        if (length(nonNaIndex) == 0) {
            stateDF <- data.frame(chr=fitData[, "chr"],
                           start=fitData[, "start"],
                           end=fitData[, "end"],
                           log2ratio=fitData[, "A"],
                           depth=fitData[, "depth"],
                           baseline=fitData[, "baseline"],
                           state=rep(NA, nrow(fitData)),
                           statUp=rep(NA, nrow(fitData)),
                           statDown=rep(NA, nrow(fitData)),
                           report=rep("average", nrow(fitData)))
            parameters <- rep(NA, (model$p + model$q + 1))
            return(list(chr = fitData[1, "chr"], state = stateDF, par = parameters))
        }
        omit <- min(which(!is.na(fitData[, "A"]))) - 1
        omitData <- fitData[1:omit, , drop = FALSE]
        fitData <- fitData[-c(1:omit), ]
    }
    yt <- matrix(fitData[, "A"], nrow = 1)
    ct <- matrix(0, nrow = 1)
    result <- kflt(yt = yt, ct = ct, model = model, smooth = s)
    if (any(is.na(yt[1, ]))) yt <- matrix(na.locf(yt[1, ]), nrow = 1)
    copyNum <- 2 ^ (1 + as.numeric(result$state))
    conf <- 1 - 0.5 * alpha/nrow(fitData)
    statUp <- as.numeric(result$state) + sqrt(as.numeric(result$variance)) * qnorm(conf)
    statDown <- as.numeric(result$state) - sqrt(as.numeric(result$variance)) * qnorm(conf)
    gainIDX <- which(statDown > 0)
    lossIDX <- which(statUp < 0)
    report <- rep("average", nrow(fitData))
    report[gainIDX] <- "gain"
    report[lossIDX] <- "loss"
    stateDF <- data.frame(chr=fitData[, "chr"],
                           start=fitData[, "start"],
                           end=fitData[, "end"],
                           log2ratio=fitData[, "A"],
                           depth=fitData[, "depth"],
                           baseline=fitData[, "baseline"],
                           state=as.numeric(result$state),
                           statUp=statUp,
                           statDown=statDown,
                           report=report)
    if (! is.null(omitData)) stateDF <- rbind(data.frame(
        chr=omitData[, "chr"],
        start=omitData[, "start"],
        end=omitData[, "end"],
        log2ratio=omitData[, "A"],
        depth=omitData[, "depth"],
        baseline=omitData[, "baseline"],
        state=rep(NA, nrow(omitData)),
        statUp=rep(NA, nrow(omitData)),
        statDown=rep(NA, nrow(omitData)),
        report=rep("average", nrow(omitData))), stateDF)
    parameters <- result$par
    return(list(chr = fitData[1, "chr"], state = stateDF, par = parameters))
}

#' @title Kalman filter.
#' @description This function allows for Kalman filtering a numeric vector use
#'   the \code{fkf()} function in the package \pkg{FKF}. See the details in the
#'   funtion \code{\link[FKF]{fkf}}.
#' @details The Kalman filter is a set of recursive linear estimation steps,
#'   which computes an estimation value at time \emph{'t+1'} based on the
#'   estimation value at time \emph{'t'} and the observation value at time
#'   \emph{'t+1'}. Its estimates are linear combinations of Gaussian data when
#'   the noise disturbances are Gaussian.
#'
#' Read depth data from Next-Generation Sequencing are proved to be a
#' combinations of Gaussian data. Perform log-ratio transformation on read depth
#' to reduce system error from sequencing process, and smooth by moving-average.
#' So the log2 ratio data is Gaussian and can be simulated by the Kalman filter.
#'
#' Based on the Box-Jenkins Approach, the Log2 ratio data fits a model of
#' \emph{ARMA(0,5)-process} with zero (p) autoregressive term and five (q)
#' moving-average terms.
#'
#' Ande parameters in detail are:
#'
#'     \emph{yt} - a matrix containing the observation and \code{NA} values are
#'     allowed. It is a \emph{1 * n} matrix with log2 ratio data of n region.
#'
#'     \emph{ct} - a matrix giving the intercept of the measurement equation. It
#'     is a \emph{1 * 1} matrix and value is 0 indicating that log2 ratio is
#'     symmetrical around 0.
#'
#'     \emph{model} - a list object containing six elements:
#'
#'         \emph{p} - a non-negative integer providing the number of
#'         autoregressive terms.
#'
#'         \emph{q} - a non-negative integer providing the number of
#'         moving-average terms.
#'
#'         \emph{d} - a non-negative integer. If d is over 0, it will difference
#'         the data to reduce the changes in the level.
#'
#'         \emph{bin.size} - a non-negative integer indicates the width of each
#'         window.
#'
#'         \emph{sn} - a numberic value providing an average signal-to-noise
#'         ratio. If input data is unstable it can be turn down.
#'
#'         \emph{method} - the method to optimize the model. By default it is
#'         \emph{nlminb} which use the function \code{\link[stats]{nlminb}}.
#'         Other methods are \emph{Nelder-Mead}, \emph{BFGS}, \emph{CG},
#'         \emph{L-BFGS-B}, \emph{SANN} and \emph{Brent} which use the function
#'         \code{\link[stats]{optim}}. See details in description of the
#'         function \code{\link[stats]{optim}}.
#'
#'     \emph{nafun} - the function name in the package \pkg{zoo} to fill a
#'     \code{NA} data point. By default is \code{\link[zoo]{na.locf}}.
#'
#'     \emph{smooth} - an non-negative integer indicates the step when do
#'     moving-average. If smooth is 0, it return the original state of Kalman
#'     filtering.
#'
#' @param yt A matrix containing the observations. \code{NA} values are allowed.
#' @param ct A matrix giving the intercept of the measurement equation.
#' @param model A list object contains elements named \emph{p}, \emph{d},
#'   \emph{q}, \emph{bin.size}, \emph{sn} and \emph{method}. See details in the
#'   function \code{\link{kflt}}. By default it is for log2 ratio data.
#' @param nafun A function's name to fill \code{NA} value. Default is
#'   \code{\link[zoo]{na.locf}}.
#' @param smooth A non-negative integer indicating the step width to smooth the
#'   state after Kalman filtering.
#' @export
#' @importFrom FKF fkf
#' @importFrom xts xts
#' @importFrom zoo na.locf
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom methods as
#' @importFrom stats dgamma
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom stats dbeta
#' @importFrom stats optim
#' @importFrom stats nlminb
#' @importFrom stats filter
#' @author Qi-Yuan Li
#' @examples
#' ####### Kalman filter #######
#' kflt(
#' yt = matrix(data = log2(rnorm(1000, mean = 2, sd = 0.1)/2), nrow = 1),
#' ct = matrix(data = 1, nrow = 1))
kflt <- function(yt, ct, model = list(p = 0, d = 1, q = 5, bin.size = 1E5, sn = 2.5, method = "nlminb"), nafun = na.locf, smooth = 0) {
    par0 <- rep(.5, model$p + model$q + 1)
    lower.bound <- c(rep(1E-6, model$p + model$q + 1))
    upper.bound <- c(rep(1, model$p + model$q), 9999)
    if (any(is.na(yt[1, ]))) yt <- matrix(nafun(yt[1, ]), nrow = 1)
    if (any(is.na(ct[1, ]))) ct <- matrix(nafun(ct[1, ]), nrow = 1)
    myt <- mean(yt[1, ], na.rm = T, trim = 0.03)
    if (model$d > 0) {
        dyt <- myDiff(yt[1, ], d = model$d)
        yt <- matrix(dyt$D, nrow = 1)
        yt0 <- matrix(dyt$X0, nrow = 1)
    }
    set.seed(sample(100000, 1))
    if (model$method == "L-BFGS-B") {

        opt <- try(optim(par = par0, fn = kfilt.pp, yt = yt, ct = ct,
                p = model$p, q = model$q, sn = model$sn,
                bin.size = model$bin.size,
                lower = lower.bound,
                upper = upper.bound,
                method = "L-BFGS-B"), silent = T)

    } else if (model$method == "nlminb") {

        opt <- nlminb(start = par0, objective = kfilt.pp, yt = yt, ct = ct,
            p = model$p, q = model$q, sn = model$sn,
            bin.size = model$bin.size,
            lower = lower.bound,
            upper = upper.bound)

    } else {
        opt <- try(optim(par = par0, fn = kfilt.pp, yt = yt, ct = ct,
            p = model$p, q = model$q, sn = model$sn,
            bin.size = model$bin.size,
            fixed = model$fixed,
            method = model$method), silent = T)
    }

    if (class(opt) == "try-error") {
        kfilts <- "No convergence"
        return(list(state = NA, var = NA, error = NA, logLik = NA))
    }
    if (opt$convergence == 0) {
        param <- opt$par
        r <- max(model$p, model$q + 1)
        if (model$p == 0) a <- rep(0, r)
        if (model$p < r && model$p > 0) a <- c(param[1:model$p], rep(0, r - model$p))
        if (model$p == r) a <- param[1:r]
        if (model$q == 0) b <- rep(0, r)
        if (model$q < r && model$q > 0) b <- c(1, param[(model$p + 1):(model$p + model$q)], rep(0, r - model$q - 1))
        sigma <- param[model$p + model$q + 1]
        omega <- sigma / model$sn
        a0 <- c(1, rep(0, r - 1))
        P0 <- diag(999, r)
        dt <- matrix(0, nrow = r)
        Tt <- cbind(a, rbind(diag(1, r - 1), rep(0, r - 1)))
        Zt <- matrix(c(1, rep(0, r - 1)), nrow = 1)
        HHt <- matrix(b, ncol = 1) %*% matrix(b, nrow = 1) * sigma ^ 2
        GGt <- matrix(omega ^ 2)
        kfilts <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = yt)
        state <- kfilts$att[1, ]
        if (model$d > 0) state <- myCSum(list(D = state, X0 = yt0))
        if (smooth > 0) state <- filter(state, filter = rep(1 / smooth, smooth), method = "convolution", sides = 2)
    } else {
        kfilts <- "No convergence"
        return(list(state = NA, variance = NA, error = NA, logLik = NA))
    }
    list(state = state,
        variance = c(NA, kfilts$Ptt[1, 1, ]),
        error = c(NA, kfilts$vt[1, ]),
        logLik = kfilts$logLik,
        par = param)
}

kfilt.pp <- function(par, yt, ct, p, q, sn, bin.size) {
    r <- max(p, q + 1)
    if (length(par) < p + q + 1) stop("Missing parameters!")
    delta <- bin.size * (1 : p) / 1E6
    lambda <- exp(-2 * delta / 100) / 2 + 0.5
    if (p == 0) a <- rep(0, r)
    if (p < r && p > 0) a <- c(par[1:p], rep(0, r - p))
    if (p == r) a <- par[1:r]
    if (q == 0) b <- rep(0, r)
    if (q < r && q > 0) b <- c(1, par[(p + 1):(p + q)], rep(0, r - q - 1))
    sigma <- par[p + q + 1]
    omega <- sigma / sn
    if (p > 0) theta <- 1E-3
    a0 <- c(1, rep(0, r - 1))
    P0 <- diag(999, r)
    dt <- matrix(0, nrow = r)
    Tt <- cbind(a, rbind(diag(1, r - 1), rep(0, r - 1)))
    Zt <- matrix(c(1, rep(0, r - 1)), nrow = 1)
    HHt <- matrix(b, ncol = 1) %*% matrix(b, nrow = 1) * sigma ^ 2
    GGt <- matrix(omega ^ 2)
    kfilt <- fkf(a0 = a0, P0 = P0, dt = dt, ct = ct, Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = yt)
    pr_a <- 0
    pr_b <- 0
    if (p > 0) {
        pr_0 <- sapply(1:p, function(i) dnorm(a[i], mean = lambda[i], sd = sigma / theta, log = T))
        pr_a <- sum(pr_0) + dgamma(1 / sigma ^ 2, shape = 2, rate = 1, log = T)
    }
    if (q > 0) pr_b <- sum(dbeta(b[-1], shape1 = 1, shape2 = 5, log = T))
    -kfilt$logLik - pr_a - pr_b
}

myDiff <- function(x, lag = 1, d = 1) {
    x0 <- rep(0, d)
    for (i in 1:d) {
        x0[i] <- x[1]
        x <- diff(x, lag = lag, differences = 1)
    }
    list(D = x, X0 = x0)
}

myCSum <- function(D) {
    x0 <- D$X0
    x <- D$D
    d <- length(x0)
    for (i in d:1) {
        x <- cumsum(c(x0[i], x))
    }
    x
}
