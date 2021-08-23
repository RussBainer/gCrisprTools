##' @title Generate a Receiver-Operator Characteristic (ROC) Curve from a CRISPR screen
##' @description Given a set of targets of interest, this function generates a ROC curve and associated statistics from the results of
##' a CRISPR screen. Specifically, it orders the elements targeted in the screen in the specified direction, and then plots the cumulative
##' proportion of positive hits on the y-axis. The corresponding vectors and Area Under the Curve (AUC) statistic are returned as a list.
##'
##' Note that ranking statistics in CRISPR screens are (usually) permutation-based, and so some granularity is expected. This
##' function does a little extra work to ensure that hits are counted as soon as the requisite value of the ranking statistic is reached
##' regardless of where the gene is located within the block of equally-significant genes. Functionally, this means that the drawn curve is
##' somewhat anticonservative in cases where the gene ranks are not well differentiated.
##'
##' @param summaryDF A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}.
##' @param target.list A character vector containing the names of the targets to be tested. Only targets contained in the \code{geneSymbol}
##' column of the provided \code{summaryDF} are considered.
##' @param direction Direction by which to order target signals (`enrich` or `deplete`).  
##' @param condense Logical indicating whether the returned x and y coordinates should be 'condensed', returning only the points at which
##' the detected proportion of \code{target.list} changes. If set to \code{FALSE}, the returned \code{x} and \code{y} vectors will explicitly
##' indicate the curve value at every position (useful for performing curve arithmetic downstream).
##' @param plot.it Logical value indicating whether to plot the curves.
##' @param ... Additional parameters for `ct.simpleResult()`
##' @return A list containing the the x and y coordinates of the curve, and the AUC statistic (invisibly).
##' @author Russell Bainer
##' @examples data('resultsDF')
##' data('essential.genes') #Note that this is an artificial example.
##' roc <- ct.ROC(resultsDF, essential.genes, direction = 'deplete')
##' str(roc)
##' @export
ct.ROC <- function(summaryDF, target.list, direction = c("enrich", "deplete"), condense = TRUE, plot.it = TRUE, ...) {

    direction <- match.arg(direction)
    stopifnot(is(plot.it, "logical"), is(condense, "logical"))

    if (!is.character(target.list)) {
        warning("Supplied target.list is not a character vector. Coercing.")
        target.list <- as.character(target.list)
    }

    # Infer whether Gsdb is ID or feature centric
    gids <- sum(target.list %in% summaryDF$geneID)
    gsids <- sum(target.list %in% summaryDF$geneSymbol)

    if (all(c(gsids, gids) == 0)) {
        stop("None of the features in the GeneSetDb are present in either the geneID or geneSymbol slots of the first provided result.")
    }

    collapse <- ifelse(gids > gsids, "geneID", "geneSymbol")
    simpleDF <- ct.simpleResult(summaryDF, collapse)

    present <- intersect(target.list, row.names(simpleDF))

    if (length(present) != length(target.list)) {
        if (length(present) < 1) {
            stop(paste0("None of the genes in the input list are present in the ", collapse, " column of the input data.frame."))
        }
        warning(paste(length(present), "of", length(target.list), "genes are present in the supplied results data.frame. Ignoring the remainder of the target.list."))
    }

    # Subset signals
    values <- simpleDF[simpleDF$direction == direction, ]
    values <- c(values$best.p[order(values$best.p, decreasing = FALSE)], rep(1, times = sum(!(simpleDF$direction %in% direction))))

    targvals <- vapply(target.list, function(x) {
        ifelse(simpleDF[x, "direction"] %in% direction, simpleDF[x, "best.p"], 1)
    }, numeric(1))

    out <- list()
    values.unique <- sort(unique(values), decreasing = FALSE)
    out$sensitivity <- unlist(lapply(values.unique, function(x) {
        sum(targvals <= x, na.rm = TRUE)/length(targvals)
    }))
    out$specificity <- unlist(lapply(values.unique, function(x) {
        (sum(values > x) - sum(targvals > x, na.rm = TRUE))/(nrow(simpleDF) - length(present))
    }))

    # Calculate the AUC/Enrichment

    out$AUC <- sum(unlist(lapply(2:length(out$sensitivity), function(x) {
        (out$sensitivity[(x)]) * (values.unique[x] - values.unique[(x - 1)])
    })))


    enrich <- switch(direction, enrich = ct.targetSetEnrichment(simpleDF, target.list, enrich = TRUE, collapse = collapse), deplete = ct.targetSetEnrichment(simpleDF, 
        target.list, enrich = FALSE, collapse = collapse))
    out <- c(out, enrich)

    # Plot it?
    if (plot.it) {
        plot(c(0, (1 - out$specificity), 1), c(0, out$sensitivity, 1), xlim = c(0, 1), ylim = c(0, 1), type = "l", ylab = "Sensitivity", xlab = "1-Specificity", main = paste("AUC:", 
            round(out$AUC, 3)), col = "blue", lwd = 3)
        abline(0, 1, lty = "dashed", col = "red")
    }

    if (!condense) {
        out <- .rocXY(out)
    }

    return(invisible(out))
}


.rocXY <- function(roc) {
    elements <- 0:max(roc$specificity)
    y <- lapply(elements, function(value) {
        pos <- length(roc$specificity[roc$specificity <= value])
        return(roc$sensitivity[pos])
    })

    return(list(x = elements, y = unlist(y), AUC = roc$AUC, targets = roc$targets, P.values = roc$P.values, Q.values = roc$Q.values))
}
