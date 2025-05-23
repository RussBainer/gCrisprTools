
## ----------------------------------------------------------------------------- File Utilities


##' Initializes the output directory
##'
##' If outdir is NULL, then no directory is checked/created. This also implies
##' that creating plots is not possible.
##'
##' @param outdir character vector pointing to a directory to check/create
##' @return TRUE if the output directory was created, otherwise FALSE (it might
##' already exist).
##' @keywords internal
##' @author Steve Lianoglou, Russell Bainer
initOutDir <- function(outdir) {
    if (is.null(outdir)) {
        return(FALSE)
    }
    if (!is.character(outdir) && length(outdir) != 1) {
        stop("character required for `outdir`")
    }
    outdir.created <- FALSE
    if (dir.exists(outdir)) {
        if (!dir.writable(outdir)) {
            stop("Can't write to output directory: ", outdir)
        }
    } else {
        pdir <- dirname(outdir)
        if (!dir.exists(pdir)) {
            stop("Path to outdir does not exist: ", pdir)
        }
        if (!dir.writable(pdir)) {
            stop("Can't create output directory in: ", pdir)
        }
        dir.create(outdir)
        outdir.created <- TRUE
    }
    outdir.created
}


##' Checks that the directory provided is writable by the current user
##'
##' This works by testing to put a temporary file into an already existing
##' directory
##'
##' @param path The path to a directory to check.
##'
##' @return \code{logical}, \code{TRUE} if \code{path} is writable by the
##' current user, otherise \code{FALSE}
##' @keywords internal
##' @author Steve Lianoglou
dir.writable <- function(path) {
    if (!dir.exists(path)) {
        stop("The directory provided does not exist: ", path)
    }
    tmp.fn <- tempfile(tmpdir = path, fileext = ".test.tmp")
    on.exit({
        if (file.exists(tmp.fn)) unlink(tmp.fn)
    })

    tryCatch({
        suppressWarnings(writeLines("test", tmp.fn))
        TRUE
    }, error = function(e) FALSE)
}


##' @title Check compatibility of a sample key with a supplied object
##' @description For many gCrisprTools functions, a sample key must be provided that specifies 
##' sample mapping to experimental groups. The sample key should be provided as a single, named factor whose  
##' names exactly correspond to the `colnames()` of the `ExpressionSet` containing the count data to be 
##' processed (or coercible as such). By convention, the first level corresponds to the control sample group.
##'  
##' This function checks whether the specified sample key is of the proper format and has 
##' properties consistent with the specified object. 
##' @param sampleKey A named factor, where the \code{levels} indicate the experimental replicate 
##' groups and the \code{names} match the \code{colnames} of the expression matrix contained in \code{object}. 
##' The first \code{level} should correspond to the control samples, but obviously there is no 
##' way to algorithmically control this. 
##' @param object An \code{ExpressionSet}, \code{EList}, or matrix.  
##' @return A logical indicating whether the objects are compatible.
##' @import limma
##' @author Russell Bainer
ct.inputCheck <- function(sampleKey, object) {

    .Deprecated('ct.keyCheck', package = 'gCrisprTools', msg = 'Please use ct.keyCheck() instead of ct.inputCheck().')
    
    if (!(is.factor(sampleKey))) {
        sampleKey <- as.factor(sampleKey)
        warning("Coercing the provided sample key to a factor. Control group set to: ", levels(sampleKey)[1])
    }
    
    # Check input formats
    if (!any(is(object, "ExpressionSet"), is(object, "EList"), is(object, "matrix"))) {
        stop(deparse(substitute(object)), " is not an ExpressionSet, Elist, or matrix. Class is: ", class(object))
    }

    if (is.null(names(sampleKey))) {
        stop(deparse(substitute(sampleKey)), " must have a names attribute, specifying the sample assignments in ", deparse(substitute(object)), ".")
    }

    # Check to see if the names match properly
    if (methods::is(object, "EList")) {
        dat <- object$E
    } else if (methods::is(object, "ExpressionSet")) {
        dat <- exprs(object)
    } else {
        dat <- object
    }

    if (!setequal(colnames(dat), names(sampleKey))) {
        stop("The names of ", deparse(substitute(sampleKey)), " must exactly match the colnames of the data contained in ", deparse(substitute(object)), ".")
    }

    return(TRUE)
}

##' @title Check compatibility of a sample key with a supplied ExpressionSet or similar object
##' @description For many gCrisprTools functions, a sample key must be provided that specifies 
##' sample mapping to experimental groups. The sample key should be provided as a single, named factor whose  
##' names exactly correspond to the `colnames()` of the `ExpressionSet` containing the count data to be 
##' processed (or coercible as such). By convention, the first level corresponds to the control sample group.
##'  
##' This function checks whether the specified sample key is of the proper format and has 
##' properties consistent with the specified object. 
##' @param sampleKey A named factor, where the \code{levels} indicate the experimental replicate 
##' groups and the \code{names} match the \code{colnames} of the expression matrix contained in \code{object}. 
##' The first \code{level} should correspond to the control samples, but obviously there is no 
##' way to algorithmically control this. 
##' @param object An \code{ExpressionSet}, \code{EList}, or other matrix-like object with defined `colnames()`.  
##' @return Invisibly, a properly formatted `sampleKey`.
##' @import limma
##' @author Russell Bainer
##' @examples data('es')
##' library(limma)
##' library(Biobase)
##' 
##' #Build the sample key
##' sk <- relevel(as.factor(pData(es)$TREATMENT_NAME), 'ControlReference')
##' names(sk) <- row.names(pData(es))
##' ct.keyCheck(sk, es)
##' @export
ct.keyCheck <- function(sampleKey, object) {
    
    #check dimensions
    stopifnot(length(sampleKey) == ncol(object))
    
    if (!(is.factor(sampleKey))) {
        sampleKey <- as.factor(sampleKey)
        warning("Coercing the provided sample key to a factor. Control group set to: ", levels(sampleKey)[1])
    }
    
    #Check that names exist, are valid, and are equal
    stopifnot(!is.null(names(sampleKey)), !is.null(colnames(object)))
    stopifnot(setequal(names(sampleKey), colnames(object)))
    
    if(any(make.names(names(sampleKey)) != names(sampleKey))) {
        stop('sampleKey names are not syntactically valid!')
    }
    
    return(invisible(droplevels(sampleKey)))
}


##' @title Determine whether a supplied object contains the results of a Pooled Screen
##' @description Many gCrisprTools functions operate on a \code{data.frame} of results
##' generated by a CRISPR screen. This function takes in a supplied object and returns 
##' a logical indicating whether the object can be treated as one of these data.frames 
##' for the purposes of downstream analyses. This is largely used internally, but can 
##' be useful if a user needs to build a result object for some reason.  
##' 
##' @param summaryDF A \code{data.frame}, usually returned by \code{ct.generateResults}. 
##' if you need to generate one of these by hand for some reason, see the example 
##' \code{resultsDF} object loaded in the example below. 
##' @return A logical indicating whether the object is of the appropriate format.
##' @author Russell Bainer
##' @examples data('resultsDF')
##' ct.resultCheck(resultsDF)
##' @export
ct.resultCheck <- function(summaryDF) {

    # Check input formats
    if (!is.data.frame(summaryDF)) {
        summaryDF <- as.data.frame(summaryDF, stringsAsFactors = FALSE)
        message("The supplied screen results are not a data.frame.")
        return(FALSE)
    }

    expectedNames <- c("geneID", "geneSymbol", "gRNA Log2 Fold Change", "gRNA Depletion P", "gRNA Depletion Q", "gRNA Enrichment P", "gRNA Enrichment Q", "Target-level Enrichment P", 
        "Target-level Enrichment Q", "Target-level Depletion P", "Target-level Depletion Q", "Median log2 Fold Change", "Rho_enrich", "Rho_deplete")

    if (!all(expectedNames %in% names(summaryDF))) {
        missing <- setdiff(expectedNames, names(summaryDF))
        warning("The supplied result object seems to have some incorrect columns. I was expecting: ")
        print(missing)
        stop("Please supply a summaryDF object generated from ct.generateResults() in the gCrisprTools package.")
        return(FALSE)
    }

    if (!setequal(vapply(summaryDF, class, character(1)), rep(c("character", "numeric"), times = c(4, 12)))) {
        stop("Some of the columns in the supplied result object seem to be of the wrong type. Please supply a summaryDF object generated from ct.generateResults() in the gCrisprTools package.")
        return(FALSE)
    }

    return(TRUE)
}

##' @title Package Screen Data into a `SummarizedExperiment` Object 
##' @description Convenience function to package major components of a screen into a `SummarizedExperiment` container
##' for downstream visualization and analysis. All arguments are optional except for `es`. 
##' @param es An `ExpressionSet` of screen data. Required. 
##' @param sampleKey a gCrisprTools `sampleKey` object, to be added to the `colData`. 
##' @param ann Annotation object to be packaged into the `rowData`
##' @param vm A `voom`-derived normalized object
##' @param fit a `MArrayLM` object containing the contrast information and model results
##' @param summaryList A named list of \code{data.frame}s, returned by \code{ct.generateResults}. 
##' if you need to generate one of these by hand for some reason, see the example 
##' \code{resultsDF} object loaded in the example below. 
##' @return A `SummarizedExperiment` object. 
##' @importFrom SummarizedExperiment SummarizedExperiment
##' @author Russell Bainer
##' @examples 
##' data('ann', 'es', 'fit', 'resultsDF')
##' ct.buildSE(es, ann = ann, fit = 'fit', summaryList = list('resA' = resultsDF, 'resB' = resultsDF))
##' @export
ct.buildSE <- function(es, sampleKey = NULL, ann = NULL, vm = NULL, fit = NULL, summaryList = NULL) {
    stopifnot(methods::is(es, "ExpressionSet"))

    asy <- list(counts = exprs(es))
    met <- list()
    rd <- fData(es)
    cd <- pData(es)

    if (!is.null(sampleKey)) {
        cd$sampleKey <- sampleKey[row.names(cd)]
    }

    if (!is.null(vm)) {
        stopifnot(setequal(colnames(vm), colnames(es)), is(vm, "EList"), setequal(row.names(vm), row.names(es)))
        asy[["voom"]] <- vm$E
        asy[["weights"]] <- vm$weights
        met[["design"]] <- vm$design

        newCols <- (ncol(cd) + 1):(ncol(cd) + ncol(vm$targets))
        cd <- cbind(cd, vm$targets[row.names(cd), ])
        colnames(cd)[newCols] <- colnames(vm$targets)
    }

    if (!is.null(fit)) {
        met[["fit"]] <- fit
    }

    if (!is.null(summaryList)) {
        if ((!is(summaryList, "list")) | (is.null(names(summaryList)))) {
            stop("When supplied, results dataframes must be provided as a named list.")
        }
        invisible(lapply(summaryList, ct.resultCheck))
        met[["results"]] <- summaryList
    }

    if (!is.null(ann)) {
        ann <- ct.prepareAnnotation(ann, es)
        rd <- cbind(rd, ann[row.names(rd), ])
    }

    se <- SummarizedExperiment::SummarizedExperiment(assays = asy, rowData = rd, colData = cd, metadata = met)

    return(se)
}

##' @title Convert a verbose results DF object to a gene-level result object
##' 
##' @description Convenience function to reduce a full results object to a gene-level object 
##' that retains minimal statistics (or alternatively, check that a provided simple result object is valid). 
##' 
##' @param summaryDF A \code{data.frame}, usually returned by \code{ct.generateResults}. 
##' if you need to generate one of these by hand for some reason, see the example 
##' \code{resultsDF} object loaded in the example below. 
##' @param collapse Column of the provided resultsDF on which to collapse values; in most cases this should be 
##' `geneSymbol` or `geneID`.
##' @return A gene-level `data.frame`, with guide-level information omitted 
##' @author Russell Bainer
##' @examples data('resultsDF')
##' ct.simpleResult(resultsDF)
##' @export
ct.simpleResult <- function(summaryDF, collapse = c("geneSymbol", "geneID")) {

    # Check inputs
    collapse <- match.arg(collapse)
    stopifnot(is(summaryDF, "data.frame"), (collapse %in% names(summaryDF)))

    if (length(setdiff(c("geneID", "geneSymbol", "Rho_enrich", "Rho_deplete", "best.p", "best.q", "direction"), (names(summaryDF)))) == 0) {
        # already simplified
        stopifnot(all(vapply(names(summaryDF)[seq_len(7)], function(x) {
            class(summaryDF[, x])
        }, character(1)) == rep(c("character", "numeric", "character"), times = c(2, 4, 1))))
        out <- summaryDF

    } else {
        stopifnot(ct.resultCheck(summaryDF))
        out <- summaryDF

        out$direction <- vapply(seq_len(nrow(out)), function(x) {
            ifelse(out[x, "Target-level Enrichment P"] < out[x, "Target-level Depletion P"], "enrich", "deplete")
        }, character(1))

        out$best.p <- vapply(seq_len(nrow(out)), function(x) {
            ifelse(out$direction[x] %in% "enrich", out[x, "Target-level Enrichment P"], out[x, "Target-level Depletion P"])
        }, numeric(1))

        out$best.q <- vapply(seq_len(nrow(out)), function(x) {
            ifelse(out$direction[x] %in% "enrich", out[x, "Target-level Enrichment Q"], out[x, "Target-level Depletion Q"])
        }, numeric(1))
        out <- out[, c("geneID", "geneSymbol", "Rho_enrich", "Rho_deplete", "best.p", "best.q", "direction")]


    }

    # Cleanup duplicates - always keep strongest signals.
    out <- out[order(out$best.p, decreasing = FALSE), ]
    out <- out[!duplicated(out[, collapse]), ]
    if (any(is.na(out[, collapse]))) {
        warning(paste0("NA detected in column ", collapse, "! Omitting the associated entries."))
        out <- out[!is.na(out[, collapse]), ]
    }
    row.names(out) <- out[, collapse]

    return(out)
}


##' @title Regularize Two Screening Results Objects 
##' @description This function prepares multiple `gCrisprTools` results dataframes for comparison. Specifically, 
##' it checks that all provided data frames are valid result objects, converts each to the target-wise `simpleResult` format, 
##' removes signals that are not shared by all objects, places their rows in identical order, and then returns the simplified dataframes as a list. 
##' 
##' This function is largely meant to be used by other gCrisprtools functions, although there are occasions when an analyst may want to call it directly. 
##' Often, it is useful to pass the `collapse` argument to `ct.simpleresult()` in cases where libraries and technologies differ between screens. 
##' @param dflist A list of results dataframes. Names will be preserved.
##' @param collapse Column of the provided resultsDFs on which to collapse values; should be `geneSymbol` or `geneID`.
##' @return A list of the in-register `simpleResult` objects, with length and names identical to `dflist`.
##' @examples 
##' data('resultsDF')
##' lapply(ct.regularizeContrasts(list('df1' = resultsDF[1:300,], 'df2' = resultsDF[200:400,])), nrow)
##' @export
ct.regularizeContrasts <- function(dflist, collapse = c("geneSymbol", "geneID")) {

    # input check
    stopifnot(is.list(dflist))

    # convert to simple results
    dflist <- lapply(dflist, ct.simpleResult, collapse = collapse)

    # find common rows
    rowcounts <- table(unlist(lapply(dflist, row.names)))
    samerows <- names(rowcounts)[rowcounts == length(dflist)]

    if (length(samerows) == 0) {
        stop("The supplied DFs have no targets in common! Consider specifying `collapse` argument for ct.simpleResult().")
    } else if (!all(rowcounts == length(dflist))) {
        message(length(samerows), " targets in common. Omitting others.")
    }

    return(lapply(dflist, function(x) {
        x[samerows, ]
    }))
}

##' @title Rank Signals in a Simplified Pooled Screen Result Object 
##' @description This function takes in a supplied results data.frame, optionally transforms it into a `simplifiedResult`, 
##' and returns the ranks of the target-level signals. 
##' 
##' @param df A results data.frame, in either raw or simplified form. Will be converted to simplified form if necessary.
##' @param top Determines the directionality of the ranking. `enrich` defines ranks from the most enriched to the most 
##' depleted target; `deplete` does the opposite 
##' @return A numeric vector of ranks, with length equal to the number of rows in the simplified data.frame. 
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' df.simple <- ct.simpleResult(resultsDF) 
##' sr <- ct.rankSimple(resultsDF)
##' all((df.simple$best.p[sr == 1] == 0), (df.simple$direction[sr == 1] == 'enrich'))
##' 
##' sr <- ct.rankSimple(resultsDF, 'deplete')
##' all((df.simple$best.p[sr == 1] == 0), (df.simple$direction[sr == 1] == 'deplete'))
##' @export
ct.rankSimple <- function(df, top = c("enrich", "deplete")) {
    df.simple <- ct.simpleResult(df)
    top <- match.arg(top)

    ranktype <- switch(top, enrich = "min", deplete = "max")

    # rank the enriched:
    enr <- rank(df.simple$best.p[df.simple$direction %in% "enrich"], na.last = TRUE, ties.method = ranktype)

    # Depleted is trickier:
    dep <- rank(df.simple$best.p[df.simple$direction %in% "deplete"], na.last = TRUE, ties.method = ranktype)
    dep <- ((max(dep) + 1) - dep) + max(enr)

    # put them back in order:
    out <- rep(0, nrow(df.simple))
    out[df.simple$direction %in% "enrich"] <- enr
    out[df.simple$direction %in% "deplete"] <- dep

    # Switch the values if depletion
    if (top %in% "deplete") {
        out <- (max(out) + 1) - out
    }

    return(out)
}

##' @title Log10 transform empirical P-values with a pseudocount
##' @description This function -log10 transforms empirical P-values by adding a pseudocount of 1/2 the minimum nonzero value. 
##' @param x numeric vector. 
##' @return -log10-transformed version of X.
##' @examples ct.softLog(runif(20))
##' @export
ct.softLog <- function(x) {
    stopifnot(is.numeric(x), all(!is.na(x)), !all(x == 0))
    out <- -log10((x + (min(x[x != 0])/2)))
    out[out <= 0] <- 0 #Sometimes correction can deflate below zero 
    return(out)
}


##' @title Extract gCrisprTools objects from a `SummarizedExperiment` 
##' @description Utility function to enable gCrisprTools functions to take `SummarizedExperiment` class objects (or subclasses thereof) as input. 
##' Throws errors when unexpected things happen. 
##' @param what What gCrisprTools-friendly object to compile. options are: 
##'   - `es`: an `ExpressionSet`, compiled from the `exprs` slot of the `assayData`, with `colData` and `rowData` saved as the `fData` and `pData`
##'   - `ann`: a gCrisprTools annotation from the `rowData`
##' @param se The `SummarizedExperiment` object. 
##' @return The specified gCrisprTools-friendly object
##' @importClassesFrom SummarizedExperiment
##' @author Russell Bainer
##' @export
ct.extractSE <- function(what, se){
  
  library('SummarizedExperiment', quietly = TRUE)
  if(!is(se, 'SummarizedExperiment')){
    stop('I tried to extract ', deparse(substitute(what)), ' from ', deparse(substitute(se)), " but it's not a SummarizedExperiment.")
  }
  
  extractable <- c('ann', 'es')
  if(!(what %in% extractable)){
    stop("I don't know how to extract ", deparse(substitute(what)), ' from a SummarizedExperiment.')
  }
  
  return(switch(what, 
                ann = ct.prepareAnnotation(rowData(se)), 
                es = as(se, 'ExpressionSet'))
  )
}
  

