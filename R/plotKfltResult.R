#' @title Plot the unifomity of coverage among samples.
#' @description Compute correlation between columns of a numeric matrix or data
#'   frame and plot the distribution of correlations.
#' @param x A numeric matrix or a data frame.
#' @return A boxplot R object.
#' @export
#' @importFrom S4Vectors cor
#' @importFrom grDevices heat.colors
#' @importFrom graphics boxplot
#' @importFrom graphics axis
#' @importFrom graphics abline
#' @author Zhan-Ni Chen
#' @examples
#' ####### Plot the coverage unifomity of samples input #######
#' covmatrix <- mergerCovFiles(c(system.file("extdata", 'testSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample.cov', package = "kfltCNV"),
#'   system.file("extdata", 'controlSample2.cov', package = "kfltCNV")))
#' plotCoverageUniformity(covmatrix[,4:6])
plotCoverageUniformity <- function(x) {
    cordata <- cor(x)
    cordata2 <- data.frame(matrix(data = NA, ncol = ncol(x), nrow = ncol(x) - 1))
    colnames(cordata2) <- colnames(x)
    i <- 1
    while (i <= ncol(x)) {
        step1 <- cordata[, i]
        step2 <- step1[-i]
        cordata2[, i] <- step2
        i <- i + 1
    }
    boxplot(cordata2, xaxt="n", na.rm="T",
        ylim = c(0,1), cex.lab = 1.2, cex.axis = 1, cex.main = 1.5,
        col = heat.colors(ncol(cordata2)),
        main = "Inter-sample Pearson Correlation Coefficient of Bin Coverage", xlab = "Sample ID",
        ylab = "Pearson Correlation Coefficient")
    axis(side = 1,at = seq(1, ncol(cordata2), 1), labels = colnames(cordata2), cex.axis = 0.7)
    abline(h = 0.5, lwd = 2, col = "blue")
    abline(h = 1.0, lwd = 2, col = "red")
}

#' @title Plot kflt results.
#' @description Plot the log2 ratio state after Kalman filtering.
#' @param file A character string of Kalman filter result file
#'   (<fileName>.state).
#' @export
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @author Zhan-Ni Chen
#' @examples
#' ####### Plot kflt results #######
#' plotKfltResult(system.file("extdata", 'testSample.state', package = "kfltCNV"))
plotKfltResult <- function(file) {
    write(paste0('Start to plot\nfile:\n', file), stdout())
    indat <- read.table(file, header = TRUE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA",
        fill = TRUE, stringsAsFactors = FALSE)
    allchr <- as.character(as.vector(indat[,1]))
    unique_chr <- unique(allchr)
    if (grepl('chr', unique_chr[1])) {
        unique_chr <- factor(unique_chr, levels = c(paste0('chr', seq(1, 22, 1)), 'chrX', 'chrY'))
    } else {
        unique_chr <- factor(unique_chr, levels = c(as.character(seq(1, 22, 1)), 'X', 'Y'))
    }
    unique_chr <- sort(unique_chr)
    unique_chr <- as.character(as.vector(unique_chr))
    for (chr in unique_chr) {
        idx <- which(allchr %in% chr)
        df <- indat[idx, ]
        if (nrow(na.omit(df)) == 0) next
        df <- df[with(df, order(start, end)), ]
        gr <- slidingWindows(GRanges(Rle(chr), IRanges(start = 1, end = max(df[, 'end']))), width = 1E7, step =  1E7)
        gr <- unlist(gr)
        pos_label <- as.character(end(gr)/1E6)
        pos_label[length(pos_label)] <- ''
        pos_label <- c('0', pos_label)
        df$report <- as.factor(df$report)
        cols <- c("average" = "gray50", "gain" = "red", "loss" = "green1")
        cols <- cols[unique(as.character(as.vector(df$report)))]
        p <- ggplot(df, aes_string(x = 'end', y = 'log2ratio', color = 'report')) +
            geom_ribbon(aes_string(ymin = 'statUp', ymax = 'statDown'),
                        fill = "gray", color = NA, alpha = 0.5) +
            geom_point(aes_string(y = "log2ratio"), na.rm = TRUE, size = 0.4) +
            scale_colour_manual(values = cols, guide = FALSE) +
            scale_y_continuous(limits = c(-4, 4)) +
            scale_x_continuous(breaks = c(0, end(gr)), labels = pos_label )+
            labs(title = paste0('chromosome ', chr),
                x = 'Chromosome Position (Mb)',
                y = 'Log2 Ratio') +
            theme_bw()+
            theme(legend.title = element_text( size = rel(2.5)), legend.position = 'bottom',
                  axis.text.x = element_text(size = 10, angle = 90))
        print(p)
    }
}
