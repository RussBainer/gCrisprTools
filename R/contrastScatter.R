##' @title Compare Two CRISPR Screen Contrasts via a Scatter Plot
##' @description This is a function for comparing the results of two screening experiments. Given two \code{summaryDF}, 
##' the function places them in register with one another, generates a simplified scatter plot where enrichment or depletion 
##' in each contrast is represented by the associated "signed" log10 (*P*/*Q*)-value (where enriched signals are represented 
##' in the positive direction and depleted signals are shown in the negative direction), and returns an invisible `data.frame`
##' containing the target X-axis and Y-axis coordinates and corresponding quadrant. 
##' 
##' This is a target-level analysis, and some minor simplifications are introduced to screen signals for the sake of clarity. 
##' Principal among these is the decision to collapse gene signals to a single directional enrichment statistic. Target-level
##' signals are typically aggregates of many guide-level signals, it is formally possible for targets to be both significantly 
##' enriched and significantly depleted within a single screen contrast as a result of substantially divergent reagent activity. 
##' This behavior is uncommon, however, and so targets are represented by selecting the direction of enrichment or depletion 
##' associated with the most significant (*P*/*Q*)-value. This directionality is then encoded into the X-axis and Y-axis 
##' position of the target as the sign of the signal as described above.
##'
##' @param dflist A (named) list of results dataframes, of length 2. See \code{\link{ct.generateResults}}. 
##' @param targets Column of the provided \code{summaryDF} to consider. Must be \code{geneID} or \code{geneSymbol}.
##' @param statistic Statistic to plot on each axis (after -log10 transformation). Must be 'p', 'q', or 'rho'.
##' @param cutoff significance cutoff used to define the significance quadrants (cannot be exactly zero).
##' @param plot.it Logical indicating whether to compose the plot on the default device. 
##' @return Invisibly, a list of length 4 containing the genes passing significance for the respective quadrants.
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' scat <- ct.scatter(list('FirstResult' = resultsDF[100:2100,], 'SecondResult' = resultsDF[1:2000,]))
##' head(scat)
##' @export
ct.scatter <- function(dflist, targets = c("geneSymbol", "geneID"), statistic = c("best.p", "best.q"), cutoff = 0.05, plot.it = TRUE) {

    # Check the input:
    targets <- match.arg(targets)
    statistic <- match.arg(statistic)
    stopifnot(is.numeric(cutoff), is.logical(plot.it), (cutoff > 0), is.list(dflist), (length(dflist) == 2))
    cutoff <- -log10(cutoff)

    dfs <- ct.regularizeContrasts(dflist, collapse = targets)
    dfs <- lapply(dfs, function(x) {
        x[, c("Rho_enrich", "Rho_deplete", "best.p", "best.q")] <- apply(x[, c("Rho_enrich", "Rho_deplete", "best.p", "best.q")], 2, ct.softLog)
        return(x)
    })
    
    # Form output & divide into quadrants:
    out <- cbind(dfs[[1]][, seq_len(6)], dfs[[2]][, 3:6])
    colnames(out) <- c("geneID", "geneSymbol", paste0(c("Rho.enrich.", "Rho.deplete.", "p.", "q."), names(dfs)[1]), paste0(c("Rho.enrich.", "Rho.deplete.", "p.", "q."), 
        names(dfs)[2]))


    c1.sig <- (dfs[[1]][, statistic] >= cutoff)  # & (dfs[[1]][,'direction'] == dfs[[2]][,'direction'])
    c2.sig <- (dfs[[2]][, statistic] >= cutoff)  # & (dfs[[1]][,'direction'] == dfs[[2]][,'direction'])

    maxval <- max(c(dfs[[1]][, statistic], dfs[[2]][, statistic]))

    out$quadrant <- 5
    out$quadrant[(c1.sig & c2.sig & (dfs[[1]]$direction == "deplete") & (dfs[[2]]$direction == "enrich"))] <- 1
    out$quadrant[(!c1.sig & c2.sig & (dfs[[2]]$direction == "enrich"))] <- 2
    out$quadrant[(c1.sig & c2.sig & (dfs[[1]]$direction == "enrich") & (dfs[[2]]$direction == "enrich"))] <- 3
    out$quadrant[(c1.sig & !c2.sig & (dfs[[1]]$direction == "deplete"))] <- 4
    out$quadrant[(c1.sig & !c2.sig & (dfs[[1]]$direction == "enrich"))] <- 6
    out$quadrant[(c1.sig & c2.sig & (dfs[[1]]$direction == "deplete") & (dfs[[2]]$direction == "deplete"))] <- 7
    out$quadrant[(!c1.sig & c2.sig & (dfs[[2]]$direction == "deplete"))] <- 8
    out$quadrant[(c1.sig & c2.sig & (dfs[[1]]$direction == "enrich") & (dfs[[2]]$direction == "deplete"))] <- 9

    out$x <- dfs[[1]][, statistic] * vapply(dfs[[1]]$direction, function(x) {
        ifelse(x == "enrich", 1, -1)
    }, numeric(1))
    out$y <- dfs[[2]][, statistic] * vapply(dfs[[2]]$direction, function(x) {
        ifelse(x == "enrich", 1, -1)
    }, numeric(1))

    # Plot it
    if (plot.it) {

        dims <- c(-maxval, maxval)

        plot(out$x, out$y, main = paste0(names(dfs)[1], " vs ", names(dfs)[2], " (", statistic, ")"), xlab = paste0(names(dfs)[1], " Signed log10(", statistic, ")"), ylab = paste0(names(dfs)[2], 
            " Signed log10(", statistic, ")"), xlim = dims, ylim = dims, pch = 19, col = rgb(0, 0, 0, 0.1), cex = 0.7)
        abline(v = c(-cutoff, cutoff), h = c(-cutoff, cutoff), lty = 2)

        points(out$x[out$quadrant %in% c(1, 3, 7, 9)], out$y[out$quadrant %in% c(1, 3, 7, 9)], pch = 19, col = rgb(0.5, 0, 0, 0.4), cex = 0.7)
        points(out$x[out$quadrant %in% c(2, 4, 6, 8)], out$y[out$quadrant %in% c(2, 4, 6, 8)], pch = 19, col = rgb(0, 0, 0.6, 0.4), cex = 0.7)

        text(-maxval, maxval, sum(out$quadrant == 1), cex = 0.5, col = rgb(0, 0, 0, 0.4))
        text(maxval, maxval, sum(out$quadrant == 3), cex = 0.5, col = rgb(0, 0, 0, 0.4))
        text(-maxval, -maxval, sum(out$quadrant == 7), cex = 0.5, col = rgb(0, 0, 0, 0.4))
        text(maxval, -maxval, sum(out$quadrant == 9), cex = 0.5, col = rgb(0, 0, 0, 0.4))
        text(-maxval, maxval, sum(out$quadrant == 1), cex = 0.5, col = rgb(0, 0, 0, 0.4))
    }
    return(invisible(out))
}













