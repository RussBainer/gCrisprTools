##' @title Identify Replicated Signals in Pooled Screens Using Conditional Scoring
##' @description This function identifies signals that are present in one or more screening experiment contrasts using a conditional 
##' strategy. Specifically, this function identifies all significant signals (according to user definitions) in a set of provided 
##' results DF and returns a `simplifiedResult` dataframe derived from the first provided contrast with an appended logical column 
##' indicating whether there is evidence for signal replication in the other provided resultsDFs.
##' 
##' Signals are considered replicated if they cross the specified stringent threshold (default: 10% FDR) in one or more of the provided 
##' contrasts, and are similarly enriched or depleted at the relaxed threshold (default: P = 0.1) in all of the remaining contrasts. If 
##' a single contrast is provided, all signals crossing the stringent threshold are conisered replicated.
##' 
##' Signals are compared across screens on the basis of \code{\link{ct.regularizeContrasts}}, so users must provide an identifier 
##' with which to standardize targets (`geneID` by default). 
##' 
##' @param dflist A list of (possibly simplified) results data.frames produced by \code{\link{ct.generateResults}}. 
##' @param statistics Statistics to use to define congruence; may be a single value, but internally coerced to a vector of length 2 where the first 
##' value corresponds to the stringent cutoff annd the second value is used for the relaxed cutoff. Must be 'best.p' or 'best.q'. 
##' @param cutoffs Numeric value(s) corresponding to the significance cutoff(s) used to define stringent and relaxed values of `statistics`. 
##' Internally coerced to a vector of length 2.
##' @param same.dir Logical vector of the same length as `dflist` indicating whether replicating signals are expected to go in the same direction 
##' (e.g., enrich/deplete in their respective screens). For example, a `dflist` of length 3 could be specified as c(TRUE, TRUE, FALSE), indicating 
##' that replicating signals should be enriched in both of the first two contrasts and depleted in the third to be considered replicated (or 
##' vise-versa). Default is `rep(TRUE, length(dflist))`.
##' @param return.stats When TRUE, return the significance of overlap instead of the logical vector (by permutation).
##' @param nperm numeric indicating number of permutations when `return.stats` is true (default 10000).  
##' @param ... Other arguments to `ct.simpleResult()`, especially `collapse`.
##' @return If `return.stats` is `FALSE`, returns the first contrast as a `simplifiedResult` data.frame, with a `replicated` logical column 
##' indicating whether each signal replicates in all of the provided screens according to the specified logic. 
##' 
##' If `return.stats` is `TRUE`, returns a dataframe indicating the permutation-based test statistics summarizing the evidence for significantly 
##' enriched signal replication across the provided contrasts (enrich, deplete, and all together).   
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' summary(ct.compareContrasts(list(resultsDF, resultsDF[1:5000,]))$replicated)
##' ct.compareContrasts(list(resultsDF, resultsDF[1:5000,]), return.stats = TRUE)
##' @export
ct.compareContrasts <- function(dflist, statistics = c("best.q", "best.p"), cutoffs = c(0.1, 0.1), same.dir = rep(TRUE, length(dflist)), return.stats = FALSE, nperm = 10000, 
    ...) {

    # Check the input:
    dflist <- ct.regularizeContrasts(dflist, ...)
    stopifnot(all(statistics %in% c("best.p", "best.q")), (length(statistics) <= 2), (length(statistics) > 0), is(return.stats, "logical"), is(cutoffs, "numeric"), 
        (length(cutoffs) <= 2), (length(cutoffs) > 0), is.logical(same.dir), (length(same.dir) == length(dflist)))
    if (length(statistics) == 1) {
        statistics <- rep(statistics, 2)
    }
    if (length(cutoffs) == 1) {
        cutoffs <- rep(cutoffs, 2)
    }

    if (return.stats) {
        stopifnot(is(nperm, "numeric"))
    }

    # Find validated signals.
    stringent <- vapply(dflist, function(x) {
        x[, statistics[1]] <= cutoffs[1]
    }, logical(nrow(dflist[[1]])))
    lax <- vapply(dflist, function(x) {
        x[, statistics[2]] <= cutoffs[2]
    }, logical(nrow(dflist[[1]])))
    dirs <- vapply(dflist, function(x) {
        x$direction %in% "enrich"
    }, logical(nrow(dflist[[1]])))
    dirs[, !same.dir] <- !dirs[, !same.dir]
    valid <- (rowSums(dirs) %in% c(ncol(dirs), 0)) & (rowSums(stringent) > 0) & (rowSums(lax) == ncol(lax))

    mainresult <- dflist[[1]]
    mainresult$replicated <- valid

    if (return.stats) {
        obs <- c(sum((rowSums(dirs) == ncol(dirs)) & valid), sum((rowSums(dirs) == 0) & valid), sum(valid))

        if (length(dflist) == 1) {
            out <- data.frame(expected = obs, observed = obs, p = c(1, 1, 1))
            row.names(out) <- c("enrich", "deplete", "all")
            return(out)
        }

        perm <- t(replicate(nperm, expr = {
            pr <- lapply(seq_len(ncol(stringent)), function(x) {
                sample(seq_len(nrow(stringent)))
            })
            nr <- nrow(stringent)

            p.stringent <- vapply(seq_len(length(pr)), function(x) {
                stringent[pr[[x]], x]
            }, logical(nr))
            p.lax <- vapply(seq_len(length(pr)), function(x) {
                lax[pr[[x]], x]
            }, logical(nr))
            p.dirs <- vapply(seq_len(length(pr)), function(x) {
                dirs[pr[[x]], x]
            }, logical(nr))

            dn <- sum((rowSums(p.dirs) == 0) & (rowSums(p.stringent) > 0) & (rowSums(p.lax) == ncol(p.lax)))
            up <- sum((rowSums(p.dirs) == ncol(p.dirs)) & (rowSums(p.stringent) > 0) & (rowSums(p.lax) == ncol(p.lax)))
            return(c(up, dn, up + dn))
        }, simplify = TRUE))
        out <- data.frame(expected = colMeans(perm), observed = obs, p = c(sum(perm[, 1] > obs[1])/nperm, sum(perm[, 2] > obs[2])/nperm, sum(perm[, 3] > obs[3])/nperm))

        row.names(out) <- c("enrich", "deplete", "all")
        return(out)
    }
    return(mainresult)
}


##' @title Consolidate shared signals across many contrasts in an UpSet Plot
##' @description This function takes in a named list of `results` dataframes produced by `ct.generateResults()` or similar, 
##' harmonizes them, and identifies overlaps between them using the logic implemented in `ct.compareContrasts()`. It then uses the
##' overlaps of these sets to compose an UpSet plot summarizing shared overlaps of the provided contrasts. These overlaps can be 
##' specified with some detail via arguments passed to the `ct.compareContrasts()` function; see documentation for more details.
##' 
##' Note that the UpSet plot is constructed to respect signal directionality, and by default constructs overlaps conditionally, 
##' but in a *bidirectional* manner. That is, a signal is considered observed in two (or more) contrasts regardless of the 
##' contrast from which the stringent signal is observed, so a signal replicated in three contrasts is interpreted as a target 
##' for which the evidence crosses the stringent threshold in one or more of the contrasts and passes the lax contrast in the others. 
##' 
##' Note that multiple important parameters are passed directly to `ct.compareContrasts()` if not specified in the command. Users 
##' are advised to study the corresponding manual page to better understand their options regarding contrast thresholding, 
##' orientation, etc. 
##' 
##' @param dflist a named list of (possibly simplified) `resultsDf`s. 
##' @param add.stats Logical indicating whether the significance of set overlaps should be included in the visualization. 
##' @param nperm Number of permutations for P-value generation. Ignored if `add.stats` is `FALSE`.
##' @param ... Other named arguments to `ComplexHeatmap::UpSet()`, `ct.compareContrasts`, or `ct.simpleResult()`. 
##' @return An UpSet plot on the current device. Silently, a combination matrix appropriate for plotting that plot, 
##' containing useful information about the observed intersections.  
##' @author Russell Bainer
##' @importFrom ComplexHeatmap %v% 
##' @examples 
##' data('resultsDF')
##' sets <- ct.upSet(list('first' = resultsDF, 'second' = resultsDF[1:5000,]))
##' @export
ct.upSet <- function(dflist, add.stats = TRUE, nperm = 10000, ...) {
    stopifnot(length(dflist) > 1, is.numeric(nperm))

    if (is.null(names(dflist))) {
        stop("The names() attribute must be set on the provided dflist for this to make any sense.")
    }

    # Subfunction args: dflist <- ct.regularizeContrasts(dflist, ...)
    dots <- list(...)
    dirs <- rep(TRUE, length(dflist))
    if ("same.dir" %in% names(dots)) {
        dirs <- dots$same.dir
    }
    dots <- dots[!(names(dots) %in% c("dflist", "return.stats", "same.dir"))]

    dflist <- do.call("ct.regularizeContrasts", args = c(list(dflist = dflist), dots[names(dots) %in% names(formals("ct.regularizeContrasts"))]))

    # Generate the relevant overlap counts
    combos <- unlist(lapply(seq_len(length(dflist)), function(x) {
        combn(seq_len(length(dflist)), x, simplify = FALSE)
    }), recursive = FALSE)

    # Calculate overlap counts
    overlaps <- lapply(combos, function(x) {
        return(do.call("ct.compareContrasts", args = c(list(dflist = dflist[x], same.dir = dirs[x], return.stats = FALSE), dots[names(dots) %in% names(formals("ct.compareContrasts"))])))
    })
    # overlaps <- lapply(combos, function(x){ if(exists('same.dir')){ dirarg <- same.dir[x] } else {dirarg <- rep(TRUE, length(x))}
    # return(ct.compareContrasts(dflist[x], same.dir = dirarg, return.stats = FALSE)) })

    overlapct <- vapply(overlaps, function(x) {
        sum(x$replicated, na.rm = TRUE)
    }, numeric(1))

    # Create the comb mat object.
    n <- length(dflist)
    comb_mat <- matrix(FALSE, nrow = n, ncol = sum(choose(n, seq_len(n))))
    for (x in seq_len(length(combos))) {
        comb_mat[combos[[x]], x] <- TRUE
    }
    rownames(comb_mat) <- names(dflist)
    comb_mat <- comb_mat + 0

    attributes(overlapct) <- NULL
    attr(comb_mat, "set_size") <- rep(nrow(dflist[[1]]), length(dflist))
    attr(comb_mat, "comb_size") <- overlapct
    attr(comb_mat, "data") <- lapply(overlaps, function(x) {
        x$geneID[x$replicated]
    })
    param <- list(mode = "conditional", value_fun = "gCrisprTools", universal_set = NULL, set_on_rows = TRUE)
    attr(comb_mat, "param") <- param
    class(comb_mat) <- c("comb_mat", "matrix")
    cmorder <- ComplexHeatmap::order.comb_mat(comb_mat)
    comb_mat <- comb_mat[cmorder]
    comb_mat <- comb_mat[overlapct[cmorder] != 0]

    # If appropriate, generate stats:
    if (add.stats) {
        combo_stats <- lapply(combos, function(x) {
            return(do.call("ct.compareContrasts", args = c(list(dflist = dflist[x], same.dir = dirs[x], return.stats = TRUE), dots[names(dots) %in% names(formals("ct.compareContrasts"))])))
        })
        # combo_stats <- lapply(combos, function(x){ if(exists('same.dir')){ dirarg <- same.dir[x] } else {dirarg <- rep(TRUE, length(x))}
        # return(ct.compareContrasts(dflist[x], same.dir = dirarg, return.stats = TRUE, nperm = nperm)) })
        lfc <- log2(vapply(combo_stats, function(x) {
            log2(x[3, 2]/x[3, 1])
        }, numeric(1))[cmorder])
        lfc[!is.finite(lfc)] <- 0
        pv <- ct.softLog(vapply(combo_stats, function(x) {
            x[3, 3]
        }, numeric(1))[cmorder])

        # Make the UpSet Plot us <- UpSet(comb_mat) %v%
        us <- do.call(getExportedValue("ComplexHeatmap", "UpSet"), args = c(list(m = comb_mat), dots[names(dots) %in% names(formals(getExportedValue("ComplexHeatmap", 
            "UpSet")))])) %v% ComplexHeatmap::HeatmapAnnotation(Log2FC = anno_points(lfc, ylim = (rep(max(abs(range(lfc))), 2) * c(-1, 1))), annotation_name_side = "left", 
            annotation_name_rot = 0) %v% ComplexHeatmap::HeatmapAnnotation(`-log10(P)` = anno_barplot(pv, ylim = c(0, (max(pv) + 1))), annotation_name_side = "left", 
            annotation_name_rot = 0)
        show(us)
        ComplexHeatmap::decorate_annotation("Log2FC", {
            pushViewport(viewport(xscale = c(0.5, 10.5), yscale = (rep(max(abs(range(lfc))), 2) * c(-1, 1))))
            grid.lines(c(0.5, 10.5), c(0, 0), gp = gpar(lty = 2), default.units = "native")
            popViewport()
        })
        ComplexHeatmap::decorate_annotation("-log10(P)", {
            pushViewport(viewport(xscale = c(0.5, 10.5), yscale = c(0, (max(pv) + 1))))
            grid.lines(c(0.5, 10.5), c(2, 2), gp = gpar(lty = 2), default.units = "native")
            popViewport()
        })
    } else {
        # show(UpSet(comb_mat))
        do.call(getExportedValue("ComplexHeatmap", "UpSet"), args = c(list(m = comb_mat), dots[names(dots) %in% names(formals(getExportedValue("ComplexHeatmap", "UpSet")))]))
    }

    return(invisible(comb_mat))

}








