##' @title Visualization of gRNA GC Content Trends
##' @description This function visualizes relationships between gRNA GC content and their measured abundance or 
##' various differential expression model estimates. 
##' @param data.obj A \code{SummarizedExperiment}, \code{ExpressionSet}, or fit (\code{MArrayLM}) object to be 
##' analyzed for the presence of GC content bias. 
##' @param ann An annotation \code{data.frame}, used to estimate GC content for each guide. Guides are annotated by 
##' row, and the object must minimally contain a \code{target} column containing a character 
##' vector that indicates the corresponding nucleotide sequences. Ignored if `data.obj` is a `SummarizedExperiment`, 
##' in which case annotation is extracted from the `rowData`.  
##' @param sampleKey An optional sample key, supplied as a factor linking the samples to experimental 
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the 
##' control set is assumed to be the first \code{level}. Ignored in the analysis of model fit objects. 
##' @param lib.size An optional vector of voom-appropriate library size adjustment factors, usually calculated with \code{\link[edgeR]{calcNormFactors}} 
##' and transformed to reflect the appropriate library size. These adjustment factors are interpreted as the total library sizes for each sample, 
##' and if absent will be extrapolated from the columnwise count sums of the \code{exprs} slot of the \code{eset}.
##' @return An image relating GC content to experimental observations on the default device. If the 
##' provided \code{data.obj} is an \code{ExpressionSet}, this takes the form of a scatter plot where the 
##' GC% and abundance estimates of each guide are represented as points overlaid
##' with a smoothed trendline within each sample. If \code{data.obj} is a fit object describing a linear model 
##' contrast, then four panels are returned describing the GC content distribution and its relationship 
##' to guide-level fold change, standard deviation, and P-value estimates.  
##' @author Russell Bainer
##' @examples 
##' data('es')
##' data('ann')
##' data('fit')
##' 
##' ct.GCbias(es, ann)
##' ct.GCbias(fit, ann)
##' @export

ct.GCbias <- function(data.obj, ann, sampleKey = NULL, lib.size = NULL) {

    # check inputs
    if (methods::is(data.obj, "SummarizedExperiment")) {
      data.obj <- ct.extractSE('es', data.obj)
      ann <- ct.extractSE('ann', data.obj)
    }
    if (methods::is(data.obj, "ExpressionSet")) {
        is.fit <- FALSE
        d <- exprs(data.obj)

        if (is.null(lib.size)) {
            lib.size <- colSums(d)
        } else if (!is.numeric(lib.size) | (length(lib.size) != ncol(d))) {
            stop("If specified, lib.size must be a numeric vector of the same 
             length as the number of samples in the eset.")
        }

        d <- t(log2(t(d + 0.5)/(lib.size + 1) * 1e+06))

        if (is.null(sampleKey)) {
            sampleKey <- factor(colnames(d))
            names(sampleKey) <- colnames(d)
        }

        sampleKey <- ct.keyCheck(sampleKey, data.obj)
        ann <- ct.prepareAnnotation(ann, object = data.obj, throw.error = FALSE)

    } else if (is(data.obj, "MArrayLM")) {

        is.fit <- TRUE
        ann <- ct.prepareAnnotation(ann, object = data.obj, throw.error = FALSE)

    } else {
        stop("Please provide an Expressionset or MArrayLM object for GC analysis.")
    }

    if (!("target" %in% names(ann))) {
        stop("The provided annotation object does not contain a \"target\" column annotating the gRNA sequences.")
    } else if (!is.character(ann$target)) {
        message("The \"target\" column of the provided annotation object is not a character vector. Coercing.")
    }

    # Calculate GC
    seq <- toupper(ann$target)
    seqlen <- vapply(seq, nchar, numeric(1))
    gc <- vapply(gregexpr("[CG]", seq, perl = TRUE), length, integer(1))/seqlen

    if (any(seqlen == 0)) {
        gc[seqlen == 0] <- NA
        warning("Could not calculate GC% for ", sum(seqlen == 0), " guides. Omitting.")
    }

    # if a fit object, plot relationship between GC and P, coefficients, variance
    if (is.fit) {
        ltp <- -log10(data.obj$p.value)
        ylims <- range(c(ltp, data.obj$coefficients[, 1]))

        par(mfrow = c(2, 2))
        hist(gc, xlim = c(0, 1), xlab = "gRNA %GC", main = "GC Content by gRNA", breaks = max(seqlen))

        plot(gc, data.obj$coefficients[, 1], xlim = c(0, 1), xlab = "gRNA %GC", ylab = "", main = "Log2 Fold Change", pch = 16, col = rgb(0, 0, 1, 0.2))
        lines(spline(gc, data.obj$coefficients[, 1]), col = "black", lwd = 2)

        plot(gc, data.obj$stdev.unscaled[, 1], xlim = c(0, 1), xlab = "gRNA %GC", ylab = "", main = "Standard Deviation", pch = 16, col = rgb(1, 0, 0, 0.2))
        lines(spline(gc, data.obj$stdev.unscaled[, 1]), col = "black", lwd = 2)

        plot(gc, ltp, xlim = c(0, 1), xlab = "gRNA %GC", ylab = "", main = "-Log10 P-value", pch = 16, col = rgb(0, 1, 0, 0.2))
        lines(spline(gc, ltp), col = "black", lwd = 2)
        return(invisible())
    } else {
        colors <- colorRampPalette(c(rgb(1, 0, 0, 0.2), rgb(0, 0, 1, 0.2)), alpha = TRUE)(ncol(d))
        spline.color <- substr(colors, 1, 7)
        plot(NA, NA, xlab = "gRNA %GC", ylab = "Log2 Counts Per Million", main = "GC Content vs. Measured Abundance", ylim = range(d), xlim = c(0, 1))
        invisible(lapply(seq_len(ncol(d)), function(x) {
            points(gc, d[, x], pch = 16, col = colors[x])
            lines(spline(gc, d[, x]), col = spline.color[x], lwd = 4)
        }))
        legend("topleft", names(sampleKey), fill = spline.color)
    }
}
