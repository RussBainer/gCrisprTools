##' @title Remove low-abundance elements from an ExpressionSet object
##' @description This function removes gRNAs only present in very low abundance across all samples of a pooled Crispr
##' screening experiment. In most cases very low-abundance guides are the
##' result of low-level contamination from other libraries, and often distort standard normalization approaches. This
##' function trims gRNAs in a largely heuristic way, assuming that the majority of 'real' gRNAs within the library are
##' comparably abundant in at least some of the samples (such as unexpanded controls), and that contaminants are
##' present at negligible levels. Specifically, the function trims the \code{trim}
##' most abundant guides from the upper tail of each log-transformed sample distribution, and then omits gRNAs whose
##' abundances are always less than 1/(2^\code{log2.ratio}) of this value.
##' @param eset An unnormalized \code{ExpressionSet} object containing, at minimum, a matrix of gRNA counts accessible with \code{exprs()}.
##' @param trim The number of gRNAs to be trimmed from the top of the distribution before estimating the abundance range. Empirically, this usually should be equal to about 2 to 5 percent of the guides in the library.
##' @param log2.ratio Maximum abundance of contaminant gRNAs, expressed on the log2 scale from the top of the trimmed range 
##' of each sample. That is, \code{log2.ratio = 4} means to discard all gRNAs whose abundance is (1/2)^4 of the trimmed maximum. 
##' @param sampleKey An (optional) sample key, supplied as an ordered factor linking the samples to experimental
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control set is assumed to be
##' the first \code{level}.
##' @param plot.it Logical value indicating whether to plot the adjusted gRNA densities on the default device.
##' @param read.floor Optionally, the minimum number of reads required for each gRNA.
##' @return An \code{ExpressionSet} object, with trace-abundance gRNAs omitted.
##' @author Russell Bainer
##' @examples data('es')
##' ct.filterReads(es)
##' @export
ct.filterReads <- function(eset, trim = 1000, log2.ratio = 4, sampleKey = NULL, plot.it = TRUE, read.floor = NULL) {

    if (!methods::is(eset, "ExpressionSet")) {
        stop(paste(deparse(substitute(eset)), "is not an ExpressionSet."))
    }
    if (!is.numeric(trim)) {
        stop("trim is not a numeric value.")
    }
    if (!is.numeric(log2.ratio)) {
        stop("log2.ratio is not a numeric value.")
    }
    if (!is.null(read.floor) & (!is.numeric(read.floor) | (length(read.floor) > 1))) {
        stop("If provided, read.floor must be a numeric value of length 1.")
    }

    e <- log2(exprs(eset) + 1)
    trim < nrow(e) || stop("'trim' must be less than the total number of features in the 'eset': ", "trim:", trim, ", total number of features: ", nrow(e), ".")
    if (!is.null(sampleKey)) {
        sampleKey <- ct.keyCheck(sampleKey, eset)
        control.samples <- names(sampleKey)[sampleKey == levels(sampleKey)[1]]
        e <- e[, control.samples]
    }
    

    # Trim and discard the elements that never cross the minimum threshold
    e.cuts <- apply(e, 2, sort, decreasing = TRUE)[trim, ] - log2.ratio
    if (!is.null(read.floor)) {
        message(paste("Using the supplied minimum threshold of", read.floor, "reads for each guide."))
        read.floor <- log2(read.floor)
        newcuts <- vapply(e.cuts, function(x) {
            ifelse(x < read.floor, read.floor, x)
        }, numeric(1))
        names(newcuts) <- names(e.cuts)
        e.cuts <- newcuts
    }

    whitelist <- row.names(e)[colSums(apply((t(e) - e.cuts), 2, sign)) != -ncol(e)]

    new.es <- eset[whitelist, ]

    # Raw
    if (plot.it) {
        par(mfrow = c(2, 1))
        ds <- apply(log2(exprs(eset) + 1), 2, density)
        ymax <- max(unlist(lapply(ds, function(x) {
            x$y
        })))
        xr <- range(unlist(lapply(ds, function(d) {
            d$x
        })))
        plot(ds[[1]], main = "Untrimmed gRNA Density", ylim = c(0, ymax), xlim = xr, xlab = "Raw Log2 gRNA Count", ylab = "Density")
        invisible(lapply(ds, lines))
        if (!is.null(sampleKey)) {
            ds <- apply(log2(exprs(eset))[, control.samples], 2, density)
            invisible(lapply(ds, lines, col = "red"))
            legend("topright", "Control", fill = "red")
        }

        # corrected
        ds <- apply(log2(exprs(new.es) + 1), 2, density)
        ymax <- max(unlist(lapply(ds, function(x) {
            x$y
        })))
        xr <- range(unlist(lapply(ds, function(d) {
            d$x
        })))
        plot(ds[[1]], main = "Trimmed gRNA Density", ylim = c(0, ymax), xlim = xr, xlab = "Trimmed Log2 gRNA Count", ylab = "Density")
        invisible(lapply(ds, lines))
        if (!is.null(sampleKey)) {
            ds <- apply(log2(exprs(new.es))[, control.samples], 2, density)
            invisible(lapply(ds, lines, col = "red"))
            legend("topright", "Control", fill = "red")
        }
    }

    return(new.es)
}





