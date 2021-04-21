##' @title Identify Replicated Signals in Pooled Screens Using Conditional Scoring
##' @description This function identifies signals in a provided screening experiment contrast that are also present in a 
##' second contrast sing a conditional strategy. Specifically, this function identifies all significant signals (according to user definitions) 
##' in a provided results DF and returns a `simplifiedResult` dataframe with an appended logical column indicating whether there 
##' is evidence for signal replication in a second provided resultsDF. 
##' 
##' Signals are considered replicated in the second screen contrast if the best signal associated with the same same target is similarly enriched 
##' or depleted in the replicating screen, and the target signal is associated with a P-value at or below a user-indicated threshold. 
##' 
##' Signals are compared across screens on the basis of \code{\link{ct.regularizeContrasts}}, so users must provide an identifier 
##' with which to standardize targets across the specified contrasts (`geneID` by default). 
##' 
##' @param mainresult A (possibly simplified) results data.frame containing the main contrast to be analyzed. See \code{\link{ct.generateResults}}. 
##' @param validationresult Column of the provided \code{summaryDF} to consider. Must be \code{geneID} or \code{geneSymbol}.
##' @param statistics Statistics to use to define congruence; may be a single value, but internally coerced to a vector of length 2 where the first 
##' value corresponds to the stringent cutoff annd the second value is used for the relaxed cutoff. Must be "best.p" or "best.q". 
##' @param cutoffs Numeric value(s) corresponding to the significance cutoff(s) used to define stringent and relaxed values of `statistics`. 
##' Internally coerced to a vector of length 2.
##' @param ... Other arguments to \link{`ct.simpleResult`}, especially `collapse`.
##' @param same.dir Logical indicating whether replicating signals are expected to go in the same direction (e.g., enrich/deplete in both screens)
##' @param return.stats When TRUE, return the significance of overlap instead of the logical vector.
##' @return If `return.stats` is `FALSE`, returns the simplified `mainresults` data.frame, with a `replicated` logical column indicating whether a 
##' signal replicates. Incomparable elements (e.g., targets not assayed in the provided `validationresult`) are set to `NA`.  If `return.stats` is 
##' `TRUE`, returns a dataframe indicating the hypergeometric test statistics summarizing the evidence for significantly enriched signal 
##' replication across the provided contrasts (enrich, deplete, and all together).   
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' summary(ct.compareContrasts(resultsDF, resultsDF[1:5000,])$replicated)
##' ct.compareContrasts(resultsDF, resultsDF[1:5000,], return.stats = TRUE)
##' @export
ct.compareContrasts  <- 
  function(mainresult,
           validationresult,
           statistics = c('best.q', 'best.p'),
           cutoffs = c(0.1, 0.1), 
           same.dir = TRUE,
           return.stats = FALSE,
           ...) {

    #Check the input: 
    #statistics <- match.arg(statistics)
    stopifnot(all(statistics %in% c('best.p', 'best.q')), (length(statistics) <= 2), (length(statistics) > 0), 
              is.numeric(cutoffs), (length(cutoffs) <= 2), (length(cutoffs) > 0), is.logical(same.dir))
    if(length(statistics) == 1){ 
      statistics <- rep(statistics, 2)
      }
    if(length(cutoffs) == 1){ 
      cutoffs <- rep(cutoffs, 2)
    }
    
    mainresult <- ct.simpleResult(mainresult, ...)
    validationresult <- ct.simpleResult(validationresult, ...)
    
    mainresult$replicated <- rep(NA, nrow(mainresult))
 
    shared <- ct.regularizeContrasts(dflist = list(df1 = mainresult, df2 = validationresult), ...)
    if(same.dir){
      valid <- ((shared$df1[,statistics[1]] < cutoffs[1]) & (shared$df2[,statistics[2]] < cutoffs[2]) & (shared$df1$direction == shared$df2$direction))
    } else {
      valid <- ((shared$df1[,statistics[1]] < cutoffs[1]) & (shared$df2[,statistics[2]] < cutoffs[2]) & (shared$df1$direction != shared$df2$direction))
    }
 
    mainresult[(row.names(mainresult) %in% row.names(shared$df1)), 'replicated'] <- FALSE
    mainresult[(row.names(mainresult) %in% row.names(shared$df1)[valid]), 'replicated'] <- TRUE
    
    if(return.stats){
      #calculate a hypergeometric test P-value for enrichment/depletion
      p.up <- gCrisprTools:::.doHyperGInternal(numW = sum(shared$df2[(shared$df2$direction %in% 'enrich'),statistics[2]] <= cutoffs[2]), 
                                numB = nrow(shared$df2), 
                                numDrawn = sum(shared$df1[(shared$df1$direction %in% 'enrich'),statistics[1]] <= cutoffs[1]), 
                                numWdrawn = sum(((mainresult$direction %in% 'enrich') & (mainresult$replicated)), na.rm = TRUE))
      p.dn <- gCrisprTools:::.doHyperGInternal(numW = sum(shared$df2[(shared$df2$direction %in% 'deplete'),statistics[2]] <= cutoffs[2]), 
                                numB = nrow(shared$df2), 
                                numDrawn = sum(shared$df1[(shared$df1$direction %in% 'deplete'),statistics[1]] <= cutoffs[1]), 
                                numWdrawn = sum(((mainresult$direction %in% 'deplete') & (mainresult$replicated)), na.rm = TRUE))
      p.all <- gCrisprTools:::.doHyperGInternal(numW = sum(shared$df2[,statistics[2]] <= cutoffs[2]), 
                                numB = nrow(shared$df2), 
                                numDrawn = sum(shared$df1[,statistics[1]] <= cutoffs[1]), 
                                numWdrawn = sum(((mainresult$replicated)), na.rm = TRUE))
      out <- data.frame('p' = c(p.up$p, p.dn$p, p.all$p), 
                        'odds.ratio' = c(p.up$odds, p.dn$odds, p.all$odds), 
                        'expected' = c(p.up$expected, p.dn$expected, p.all$expected), 
                        'observed' = c(sum(((mainresult$direction %in% 'enrich') & (mainresult$replicated)), na.rm = TRUE), 
                                       sum(((mainresult$direction %in% 'deplete') & (mainresult$replicated)), na.rm = TRUE),
                                       sum(((mainresult$replicated)), na.rm = TRUE)))
      row.names(out) <- c('enrich', 'deplete', 'all') 
      return(out)
    }
    return(mainresult)
  }


##' @title Consolidate shared signals across many contrasts in an UpSet Plot
##' @description This function takes in a named list of `results` dataframes produced by `ct.generateResults()` or similar, 
##' harmonizes them, and identifies overlaps between them using the logic implemented in `ct.compareContrasts()`. It then uses the
##' Overlaps of these sets to compose an UpSet plot summarizing shared overlaps of the provided contrasts. These overlaps can be 
##' specified with some detail via arguments passed to the `ct.compareContrasts()` function; see documentation for more details.
##' 
##' Note that the UpSet plot is constructed to respect signal directionality, and by default constructs overlaps conditionally, 
##' but in a *bidirectional* manner. That is, a signal is considered observed in two (or more) contrasts regardless of the 
##' contrast from which the stringent signal is observed, so a signal replicated in three contrasts is interpreted as a target 
##' for which the evidence crosses the stringent threshold in one or more of the contrasts and passes the lax contrast in the others. 
##' 
##' @param dflist a named list of (possibly simplified) `resultsDf`s. 
##' @param orientation Optionally, a numeric vector encoding the orientation of the contrasts relative to one another. 
##' @param mode Mode of intersection. "intersect" by default; see `ComplexHeatmap::make_comb_mat()` for details.
##' @param ... Other named arguments to `ComplexHeatmap::UpSet()`, `ct.compareContrasts`, or `ct.simpleResult()`. 
##' @return Silently, a named list indicating the set of targets shared within each pair of contrasts.    
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' sets <- ct.contrastUpset(list('first' = resultsDF, 'second' = resultsDF[1:5000,]))
##' @export
ct.contrastUpset <- function(dflist, 
                             orientation = NULL, 
                             statistic = c('best.q', 'best.p'), 
                             cutoff = 0.1,
                             mode = c("intersect", "union", "distinct"), 
                             ...){
  suppressPackageStartupMessages(library(ComplexHeatmap, quietly = TRUE))
  
  dflist <- ct.regularizeContrasts(dflist, ...)
  if(is.null(names(dflist))){
    stop('The names() attribute must be set on the provided dflist for this to make any sense.')
  }
  
  mode <- match.arg(mode)
  statistic <- match.arg(statistic)
  stopifnot(is(cutoff, 'numeric'), cutoff <= 1, cutoff >= 0)
  if(!is.null(orientation)){
    stopifnot(length(orientation) == length(dflist), )
  }
  
  #prep the combinatorial matrix
  m <- ComplexHeatmap::make_comb_mat(setNames(as.list(1:length(dflist)), names(dflist)), mode = mode)
  
  #Compile a list of comparisons
  combos <- combn(names(dflist), 2)
  combos <- cbind(combos, combos[2:1,])
  message(paste0(ncol(combos), ' conditional comparisons defined. Compiling lists.'))
  
  hits <- vapply(1:ncol(combos), 
                 function(x){
                   ct.compareContrasts(dflist[combos[1,x]], dflist[combos[2,x]], return.stats = FALSE, ...)$replicated
                 }, logical(nrow(dflist[[1]])))
  
  
  
  
  
  
}










