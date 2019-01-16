#' @title Perform call CNVs.
#' @description Perform call CNVs from log2 ratio state file (<filename>.state).
#'   If \emph{annote.database} isn't \code{NULL}, it will create the gene CNV
#'   result.
#' @details The results are tabular or
#'   \href{https://samtools.github.io/hts-specs/VCFv4.2.pdf}{VCFv4.2} format
#'   file.
#'
#' The tabular file has nine columns while a gene CNV result will also have a column named \emph{gene}.
#'
#'     \emph{chr} - chromosome name.
#'
#'     \emph{start} - start position of CNV segment.
#'
#'     \emph{end} - end position of CNV segment.
#'
#'     \emph{svtype} - \emph{DUP} or \emph{DEL} representing copy number gains or losses in CNV segment.
#'
#'     \emph{log2ratio} - log2 ratio of CNV segment.
#'
#'     \emph{probe} - the number of bins covered by CNV segment.
#'
#'     \emph{depth} - median depth of CNV segment in test sample.
#'
#'     \emph{baseline} - median depth of CNV segment in baseline.
#'
#'     \emph{cn} - copy number of CNV segment.
#'
#'     \emph{gene} - gene symbol.
#' @param file A character string of a log2 ratio state file (<filename>.state)
#'   path.
#' @param outPrefix A character string, indicates output files' prefix. By
#'   default it is the prefix of input file.
#' @param id A character string of sample id in VCF file. By default it is the
#'   prefix of input file.
#' @param gapwidth A non-negative integer indicates the maximum gap width
#'   between CNV region to ignore.
#' @param annote.database A character string of \emph{refGene.txt} file path. It
#'   can be download from \href{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database}{UCSC database ftp}.
#'   And \emph{refGene.txt} from \emph{Annovar} database is supported, too.
#' @param out.type A character string, \emph{vcf}, \emph{table} or \emph{all}.
#'   By default it will export both VCF file and tabular file.
#' @param threshold A non-negative numeric value. A CNV of which log2 ratio is
#'   out of range from inverse \emph{threshold} to \emph{threshold} will not
#'   export.
#' @param min.probes A non-negative integer. A CNV of which probe is under
#'   \emph{min.probes} will not export.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Perform call CNVs #######
#' performCallCNV(
#'     file =  system.file("extdata", 'testSample.state', package = "kfltCNV"),
#'     annote.database = system.file("extdata", 'hg19_refGene_chr10.txt', package = "kfltCNV"))
performCallCNV <- function(file, outPrefix = NULL, id = NULL, gapwidth = 500, annote.database = NULL, out.type = 'all', threshold = 0, min.probes = 1) {

    if (is.null(id)) {
        id <- rev(unlist(strsplit(file, "/")))[1]
        id <- unlist(strsplit(file, "\\."))[1]
    }
    if (is.null(outPrefix)) {
        outPrefix <- id
    }
    if (length(which(! out.type %in% c('vcf', 'table', 'all'))) == length(out.type)) stop('Error in out.type.')
    if (! is.null(annote.database)) {
        if (! file.exists(annote.database)) stop(paste0(annote.database, ' file no exists.'))
    }
    if ('all' %in% out.type) out.type <- c('vcf', 'table')

    write(paste0('Start calling CNVs\nFile:\n', file,
        '\nMerge CNV gap small than gapwidth:\t', gapwidth,
        ' bp\nLog2 ratio threshold:\t', threshold,
        '\nMinimun probe counts supporting CNV:\t', min.probes, '\nOutput file type:\t',
        paste(out.type, collapse = ', ')), stdout())

    gainDF <- callGainCNV(file, gapwidth = gapwidth)
    lossDF <- callLossCNV(file, gapwidth = gapwidth)

    if (is.null(gainDF) & is.null(lossDF)) return()
    if (is.null(gainDF)) {
        cnvDF <- lossDF
    } else if (is.null(lossDF)) {
        cnvDF <- gainDF
    } else {
        cnvDF <- rbind(gainDF, lossDF)
    }
    cnvDF <- cnvDF[with(cnvDF, order(cnvDF$chr, cnvDF$start, cnvDF$end)), ]
    info.vec <- c("svtype", "svlen", "end", "log2ratio", "probe", "depth", "baseline")
    if ('table' %in% out.type) outputTable(cnvDF, filepath = paste0(outPrefix, '.cnv.txt'), info.vec = info.vec, threshold = threshold, min.probes = min.probes)
    if ('vcf' %in% out.type) outputVcf(cnvDF, filepath = paste0(outPrefix, '.cnv.vcf'), id = id, info.vec = info.vec, threshold = threshold, min.probes = min.probes )
    if (! is.null(annote.database) ){
        annote.file <- annovateBedFormat(x = file, y = paste0(outPrefix, '.anno.state'), annote.database = annote.database)
        geneDF <- callCNVGene(file = paste0(outPrefix, '.anno.state'), gapwidth = gapwidth)
        info.vec2 <- c("svtype", "end", 'gene', "log2ratio", "probe", "depth", "baseline")
        if ('table' %in% out.type) outputTable(geneDF, filepath = paste0(outPrefix, '.cnv.gene.txt'), info.vec = info.vec2, threshold = threshold, min.probes = min.probes)
        if ('vcf' %in% out.type) outputVcf(geneDF, filepath = paste0(outPrefix, '.cnv.gene.vcf'), id = id, info.vec = info.vec2, threshold = threshold, min.probes = min.probes, addinfo = '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">')
    }
}

#' @title Call CNVs.
#' @description Read a log2 ratio state file (<filename>.state) and call CNVs.
#'
#' \code{callGainCNV()} will merge and output copy number gains.
#'
#' \code{callLossCNV()} will merge and output copy number losses.
#'
#' \code{callCNVGene()} will compute gene's copy number in the CNV segment. It need an input file with a column named \emph{gene}. It can be created by a log2 ratio state file (<filename>.state) using the function \code{\link{annovateBedFormat}}.
#' @details The result is a data frame with nine columns and a \emph{gene} column if \code{callCNVGene()} is used:
#'
#'     \emph{chr} - chromosome name.
#'
#'     \emph{start} - start position of CNV segment.
#'
#'     \emph{end} - end position of CNV segment.
#'
#'     \emph{svtype} - \emph{DUP} or \emph{DEL} representing copy number gains or losses in CNV segment.
#'
#'     \emph{log2ratio} - log2 ratio of CNV segment.
#'
#'     \emph{probe} - the number of bins covered by CNV segment.
#'
#'     \emph{depth} - median depth of CNV segment in test sample.
#'
#'     \emph{baseline} - median depth of CNV segment in baseline.
#'
#'     \emph{cn} - copy number of CNV segment.
#'
#'     \emph{gene} - gene symbol.
#' @param file A character string of a log2 ratio state file (<filename>.state)
#'   path.
#' @param gapwidth A non-negative integer indicates the maximum gap width
#'   between CNV region to ignore.
#' @return If find no CNVs it will retrun NULL. Or else return a data frame with
#'   nine columns named \emph{chr}, \emph{start}, \emph{end}, \emph{svtype},
#'   \emph{log2ratio}, \emph{probe}, \emph{depth}, \emph{baseline} and
#'   \emph{cn}.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import S4Vectors
#' @importFrom GenomicAlignments findOverlaps
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom methods as
#' @importFrom stats setNames
#' @author Zhan-Ni Chen
#' @examples
#' ####### Call CNVs #######
#' callGainCNV(system.file("extdata", 'testSample.state', package = "kfltCNV"), gapwidth = 0)
callGainCNV <- function(file, gapwidth) {
    stateDF <- read.table(file, header = TRUE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
    if (length(which(stateDF[, 'report'] %in% "gain")) == 0) return(NULL)
    gainDF <- stateDF[which(stateDF[, 'report'] %in% "gain"), ]
    gainDF.list <- split(gainDF, as.factor(gainDF[, "chr"]))
    stat.vec <- c("log2ratio", "depth", "baseline", "state", "statUp", "statDown")
    method <- "median"
    gain.merge.list <- list()
    for (i in 1:length(gainDF.list)) {
        gain.merge.list[[i]] <- mergeDF(gainDF.list[[i]], stat.vec = stat.vec, method = method, gapwidth = gapwidth)
    }
    gain.merge.list <- do.call("rbind", gain.merge.list)
    if (nrow(gain.merge.list) > 0) {
        gain.merge.list$svtype <- rep("DUP", nrow(gain.merge.list))
        gain.merge.list$svlen <- as.numeric(as.vector(gain.merge.list[,"end"])) - as.numeric(as.vector(gain.merge.list[,"start"])) + 1
        gain.merge.list$cn <- round(2*2^as.numeric(as.vector(gain.merge.list[,"log2ratio"])), 0)
        gain.merge.list$depth <- round(as.numeric(as.vector(gain.merge.list[,"depth"])), 0)
        gain.merge.list$baseline <- round(as.numeric(as.vector(gain.merge.list[,"baseline"])), 0)
        gain.merge.list$log2ratio <- round(as.numeric(as.vector(gain.merge.list[,"log2ratio"])), 6)
        return(gain.merge.list)
    } else {
        return(NULL)
    }
}

#' @rdname callGainCNV
#' @export
callLossCNV <- function(file, gapwidth) {
    stateDF <- read.table(file, header = TRUE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
    if (length(which(stateDF[, 'report'] %in% "loss")) == 0) return(NULL)
    lossDF <- stateDF[which(stateDF[, 'report'] %in% "loss"), ]
    lossDF.list <- split(lossDF, as.factor(lossDF[, "chr"]))
    stat.vec <- c("log2ratio", "depth", "baseline", "state", "statUp", "statDown")
    method <- "median"
    loss.merge.list <- list()
    for (i in 1:length(lossDF.list)) {
        loss.merge.list[[i]] <- mergeDF(lossDF.list[[i]], stat.vec = stat.vec, method = method, gapwidth = gapwidth)
    }
    loss.merge.list <- do.call("rbind", loss.merge.list)
    if (nrow(loss.merge.list) > 0) {
        loss.merge.list$svtype <- rep("DEL", nrow(loss.merge.list))
        loss.merge.list$svlen <- - (as.numeric(as.vector(loss.merge.list[, "end"])) - as.numeric(as.vector(loss.merge.list[, "start"])) + 1)
        loss.merge.list$cn <- round(2*2^as.numeric(as.vector(loss.merge.list[, "log2ratio"])), 0)
        loss.merge.list$depth <- round(as.numeric(as.vector(loss.merge.list[,"depth"])), 0)
        loss.merge.list$baseline <- round(as.numeric(as.vector(loss.merge.list[,"baseline"])), 0)
        loss.merge.list$log2ratio <- round(as.numeric(as.vector(loss.merge.list[,"log2ratio"])), 6)
        return(loss.merge.list)
    } else {
        return(NULL)
    }
}

#' @rdname callGainCNV
#' @export
callCNVGene <- function(file, gapwidth) {
    stateDF <- read.table(file, header = TRUE, sep = "\t", quote = "", comment.char = "#",
        na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
    cnvDF <- stateDF[which(stateDF[, 'report'] %in% c("gain", "loss")), ]
    if (length(cnvDF) > 0) {
        cnvgr <- GRanges(Rle(as.character(as.vector(cnvDF[, 'chr']))),
            IRanges(start = as.numeric(as.vector(cnvDF[,'start'])), end = as.numeric(as.vector(cnvDF[,'end']))))
        gainDF <- callGainCNV(file, gapwidth = gapwidth)
        lossDF <- callLossCNV(file, gapwidth = gapwidth)
        cnvSegmentDF <- rbind(gainDF, lossDF)
        cnvSegment.gr <- GRanges(Rle(as.character(as.vector(cnvSegmentDF[, 'chr']))),
            IRanges(start = as.numeric(as.vector(cnvSegmentDF[,'start'])),
                end = as.numeric(as.vector(cnvSegmentDF[,'end']))),
            probe=as.numeric(as.vector(cnvSegmentDF[, 'probe'])))
        probes <- rep(NA, nrow(cnvDF))
        hit <- findOverlaps(cnvgr, cnvSegment.gr, type = "within")
        probes[queryHits(hit)] <- cnvSegment.gr$probe[subjectHits(hit)]
        cnvDF$cnvprobe <- probes
        cnvDF.list <- split(cnvDF, as.factor(cnvDF[, "gene"]))
        depth <- sapply(cnvDF.list, function(x) {median(x[,'depth'], na.rm = TRUE)} )
        baseline <- sapply(cnvDF.list, function(x) {median(x[,'baseline'], na.rm = TRUE)} )
        depth[which(depth <= 1)] <- 1
        baseline[which(baseline <= 1)] <- 1
        log2ratio <- log2(depth/baseline)
        svtype <- rep("DUP", length(cnvDF.list))
        svtype[which(log2ratio<0)] <- "DEL"
        df <- data.frame(chr=sapply(cnvDF.list, function(x) {x[1, 'chr']} ),
            start=sapply(cnvDF.list, function(x) {x[1, 'start']} ),
            end=sapply(cnvDF.list, function(x) {x[nrow(x), 'end']} ),
            svtype=svtype,
            log2ratio=round(log2ratio, 6),
            bins=lengths(cnvDF.list),
            probe=sapply(cnvDF.list, function(x) {max(x[, 'cnvprobe'], na.rm = TRUE)} ),
            depth=round(sapply(cnvDF.list, function(x) {median(x[,'depth'], na.rm = TRUE)} ), 0),
            baseline=round(sapply(cnvDF.list, function(x) {median(x[,'baseline'], na.rm = TRUE)} ), 0),
            cn=round(2 * 2 ^ log2ratio, 0),
            gene=names(cnvDF.list))
        df[with(df, order(df$chr, df$start, df$end)), ]
    }
}

#' @title Annovating BED format file with gene symbol.
#' @description Annovating BED format file with gene symbol. It will create a
#'   new file with a column named \emph{gene} adding to input file.
#' @param x A character string of file path.
#' @param y A character string of a file to write to.
#' @param annote.database A character string of \emph{refGene.txt} file path. It
#'   can be download from \href{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database}{UCSC database ftp}.
#'   And \emph{refGene.txt} from \emph{Annovar} database is supported, too.
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import S4Vectors
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom methods as
#' @importFrom GenomicAlignments findOverlaps
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' ####### Annovating BED format file with gene symbol #######
#' annovateBedFormat(x = system.file("extdata", 'testSample.state', package = "kfltCNV"),
#'   y = 'test.anno.state',
#'   annote.database = system.file("extdata", 'hg19_refGene_chr10.txt', package = "kfltCNV"))
annovateBedFormat <- function(x, y, annote.database) {
    indat <- read.table(x, header = TRUE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
    chr.flag <- grepl("chr", indat[1,1])
    indat.gr <- GRanges(Rle(as.character(as.vector(indat[,1]))),
                IRanges(start = as.numeric(as.vector(indat[,2])), end = as.numeric(as.vector(indat[,3]))))
    mcols(indat.gr) <- indat
    dbdat <- read.table(annote.database, header = FALSE, sep = "\t", quote = "", comment.char = "#", na.strings = "NA", fill = TRUE, stringsAsFactors = FALSE)
    chr <- as.character(as.vector(dbdat[,3]))
    idx <- which(chr %in% c(paste0("chr", seq(1, 22, 1)), "chrX", "chrY"))
    dbdat <- dbdat[idx,]
    chr <- as.character(as.vector(dbdat[,3]))
    if (! chr.flag) chr <- gsub("chr", "", chr)
    db.gr <- GRanges(Rle(chr),
                IRanges(start = as.numeric(as.vector(dbdat[,5])), end = as.numeric(as.vector(dbdat[,6]))),
                gene=as.character(as.vector(dbdat[,13])))
    db.gr1 <- db.gr[unique(subjectHits(findOverlaps(indat.gr, db.gr)))]
    indat.gr.list <- GenomicRanges::split(as(indat.gr, "GAlignments"), seq(1:length(indat.gr)))
    names(indat.gr.list) <- paste0(as.character(seqnames(indat.gr)), ":", start(indat.gr), "-", end(indat.gr))
    db.gr.list <- splitGrByGRlist(gr = db.gr1, grlist = indat.gr.list, split.type = "any")
    gene.list <- split(x = unlist(db.gr.list)$gene, f = names(unlist(db.gr.list)))
    genes <- lapply(gene.list, function(x, ...) { paste(unique(x), collapse = ",")})
    allgene <- rep(NA, length(indat.gr))
    names(allgene) <- names(indat.gr.list)
    allgene[names(genes)] <- genes[names(genes)]
    df <- mcols(indat.gr)
    df$gene <- unlist(allgene)
    write.table(df, file=y, row.names=FALSE, col.names = TRUE, quote=FALSE, sep="\t")
}

#' @title Export CNVs result.
#' @description Write CNVs result to tabular or VCF format file
#' @param x A data frame which is result of calling CNVs.
#' @param filepath A character string of file path to write to.
#' @param info.vec A character vector of column names in x to export.
#' @param threshold A non-negative numeric value. A CNV of which log2 ratio is
#'   out of range from inverse \emph{threshold} to \emph{threshold} will not
#'   export.
#' @param min.probes A non-negative integer. A CNV of which probe is under
#'   \code{min.probes} will not export.
#' @param id A character string of the sample id in VCF file.
#' @param addinfo A character string or vecter of column names in x to add to
#'   VCF file header and \emph{INFO} column. Default header can be created by
#'   \code{makeVcfHeader()}.
#' @param addformat A character string or vecter of column names in x to add to
#'   VCF file header, \emph{FORMAT} and sample column. Default header can be
#'   created by running \code{makeVcfHeader}.
#' @export
#' @author Zhan-Ni Chen
#' @examples
#' loss <- callLossCNV(system.file("extdata", 'testSample.state', package = "kfltCNV"), 0)
#' gene <- callCNVGene(system.file("extdata", 'testSample.anno.state', package = "kfltCNV"), 0)
#'
#' # tabular
#' outputTable(loss, filepath = 'cnv.result.txt',
#'   info.vec =  c("svtype", "svlen", "end", "log2ratio", "probe", "depth", "baseline"))
#'
#' # vcf file
#' outputVcf(loss, filepath = 'cnv.result.vcf', id = 'testSample',
#'   info.vec =  c("svtype", "svlen", "end", "log2ratio", "probe", "depth", "baseline"))
#'
#' # vcf file adding gene
#' outputVcf(gene, filepath = 'cnv.result.vcf', id = 'testSample',
#'   info.vec =  c("svtype", "end", 'gene', "log2ratio", "probe", "depth", "baseline"),
#'   addinfo = '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">')
outputTable <- function(x, filepath, info.vec, threshold = 0.5, min.probes = 30) {
    if (! is.null(x)) {
        x <- na.omit(x)
        idx <- which((as.numeric(as.vector(x[, "log2ratio"])) >= threshold |
            as.numeric(as.vector(x[, "log2ratio"])) <= -1 * threshold) &
            as.numeric(as.vector(x[, "probe"])) >= min.probes)
        if ( length(idx) > 0) {
            write.table(x[idx, c("chr", "start", "end", setdiff(info.vec, c("chr", "start", "end", "cn")), "cn")], file = filepath, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        }
    }
}

#' @rdname outputTable
#' @export
outputVcf <- function(x, filepath, id, info.vec, threshold = 0.5, min.probes = 30, addinfo = NULL, addformat = NULL) {
    if (! is.null(x)) {
        x <- na.omit(x)
        idx <- which((as.numeric(as.vector(x[, "log2ratio"])) >= threshold |
            as.numeric(as.vector(x[, "log2ratio"])) <= -1 * threshold) &
            as.numeric(as.vector(x[, "probe"])) >= min.probes)
        if ( length(idx) > 0) {
            x <- x[idx, , drop = FALSE]
            header <- makeVcfHeader(addinfo = addinfo, addformat = addformat)
            vcf.matrix <- makeVcfMatrix(x, info.vec)
            cat(header, file=filepath, sep="\n")
            cat(paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", id), collapse="\t"),
             file=filepath, sep="\n", append = TRUE)
            cat(apply(vcf.matrix, 1, paste, collapse="\t"), file=filepath, sep="\n", append = TRUE)
        }
    }
}

mergeDF <- function(x, stat.vec, method, gapwidth) {
    if ( length( which( ! stat.vec %in% colnames(x))) > 0) stop("Values to stat do not in the data frame.")
    gr <- GRanges(Rle(as.character(as.vector(x[, "chr"]))),
            IRanges(start = as.numeric(as.vector(x[, "start"])), end = as.numeric(as.vector(x[, "end"]))))
    mcols(gr) <- x[, stat.vec]
    gr <- sort(gr)
    gr.merge <- reduce(gr, min.gapwidth = gapwidth, ignore.strand = TRUE)
    gr.merge.list <- GenomicRanges::split(gr.merge, seq(1, length(gr.merge), 1))
    gr.each.list <- splitGrByGRlist(gr = gr, grlist = gr.merge.list, split.type = "within")
    probes <- lengths(gr.each.list)
    merge.dat <- mergeGR(gr.each.list, stat.vec = stat.vec, method = method)
    merge.dat <- cbind(data.frame(chr=as.character(as.vector(seqnames(gr.merge))),
        start=start(gr.merge),
        end=end(gr.merge)), merge.dat)
    merge.dat$probe <- probes
    merge.dat
}

splitGrByGRlist <- function(gr, grlist, split.type) {
    hits <- findOverlaps(gr, grlist, ignore.strand = TRUE, type = split.type)
    gr2grlist <- setNames(as(t(hits), "List"), names(grlist))
    gr_by_grlist <- extractList(gr, gr2grlist)
    gr_by_grlist
}

mergeGR <- function(gr.list, stat.vec, method) {
    domerge <- function(gr, stat.vec = NULL, method = NULL, pst.vec = NULL, ...){
        if (length(gr) > 0) {
            stat.value <- apply(mcols(gr)[, stat.vec, drop = FALSE], 2, method, na.rm = TRUE)
        } else {
            stat.value <- rep(NA, length(stat.vec))
        }
      return(stat.value)
    }
    result <- lapply(gr.list, domerge, stat.vec = stat.vec, method = method)
    result <- do.call("rbind", result)
    result <- as.data.frame(result)
    colnames(result) <- stat.vec
    result
}

makeVcfMatrix <- function(x, info.vec) {
    CHROM <- as.character(as.vector(x[, "chr"]))
    POS <- as.numeric(as.vector(x[, "start"]))
    if ("id" %in% names(x)) {
        ID <- as.character(as.vector(x[, "id"]))
    } else {
        ID <- rep(".", nrow(x))
    }
    REF <- rep("N", nrow(x))
    ALT <- paste0("<", as.character(as.vector(x[, "svtype"])), ">")
    if ("qual" %in% names(x)) {
        QUAL <- as.character(as.vector(x[, "qual"]))
    } else {
        QUAL <- rep(".", nrow(x))
    }
    if ("filter" %in% names(x)) {
        FILTER <- as.character(as.vector(x[, "filter"]))
    } else {
        FILTER <- rep(".", nrow(x))
    }
    idx <- which(info.vec %in% colnames(x))
    info.vec <- info.vec[idx]
    if (length(info.vec) == 0) stop("Error info.vec: column names do not contain any string in info.vec.")
    INFO <- makeVcfInfo(x[, info.vec])
    x$gt <- ifelse(as.integer(as.vector(x[,'cn'])) == 0, "1/1", "0/1")
    if(! "gq" %in% names(x)) {
        x$gq <- rep('.', nrow(x))
    }
    FORMAT <- makeVcfFormat(x[, c("gt", "gq", "cn")])
    cbind(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, rep("GT:GQ:CN", nrow(x)), FORMAT)
}

makeVcfInfo <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    col_n <- colnames(x)
    x <- apply(x, 2, function(y) { as.character(as.vector(y)) })
    x <- data.frame(matrix(data = x, ncol = n_col, nrow = n_row))
    colnames(x) <- col_n
    sapply(1:nrow(x), function(i) {
        paste(paste0(toupper(colnames(x)), "=", as.character(as.vector(x[i,]))), collapse=";")
    })
}

makeVcfFormat <- function(x) {
    n_col <- ncol(x)
    n_row <- nrow(x)
    col_n <- colnames(x)
    x <- apply(x, 2, function(y) { as.character(as.vector(y)) })
    x <- data.frame(matrix(data = x, ncol = n_col, nrow = n_row))
    colnames(x) <- col_n
    sapply(1:nrow(x), function(i) {
        paste(as.character(as.vector(x[i,])), collapse=":")
    })
}

makeVcfHeader <- function(addinfo = NULL, addformat = NULL) {
    fileformat <- "##fileformat=VCFv4.2"
    fileDate <- paste0("##fileDate=", format(Sys.time(), "%Y%m%d"))
    alt <- c('##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=DUP,Description="Duplication">')
    info <- c('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=DEPTH,Number=1,Type=Float,Description="Weighted mean of normalized read depths across all this gene bins.">',
    '##INFO=<ID=BASELINE,Number=1,Type=Float,Description="Weighted mean of normalized read depths across control samples.">',
    '##INFO=<ID=LOG2RATIO,Number=1,Type=Float,Description="Weighted mean of log2 ratios of all the gene bins, including any off-target intronic bins.">',
    '##INFO=<ID=PROBE,Number=1,Type=Integer,Description="Number of probes in CNV">')
    if (! is.null(addinfo)) info <- unique(c(info, addinfo))
    format <- c('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">')
    if (! is.null(format)) format <- unique(c(format, addformat))
    c(fileformat, fileDate, alt, info, format)
}
