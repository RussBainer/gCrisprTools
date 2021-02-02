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
##' `TRUE`, returns a named list indicating the hypergeometric test *P*-values summarizing the evidence for significantly enriched signal 
##' replication across screens (enrich, deplete, and all together).   
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
    statistics <- match.arg(statistics)
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
      valid <- ((shared$df1[,statistics[1]] < cutoffs[1]) & (shared$df1[,statistics[2]] < cutoffs[2]) & (shared$df1$direction == shared$df2$direction))
    } else {
      valid <- ((shared$df1[,statistics[1]] < cutoffs[1]) & (shared$df1[,statistics[2]] < cutoffs[2]) & (shared$df1$direction != shared$df2$direction))
    }
 
    mainresult[row.names(shared$df1), 'replicated'] <- FALSE
    mainresult[row.names(shared$df1)[valid], 'replicated'] <- TRUE
    
    if(return.stats){
      #calculate a hypergeometric test P-value for enrichment/depletion
      p.up <- .doHyperGInternal(numW = sum(shared$df2[(shared$df2$direction %in% 'enrich'),statistics[2]] <= cutoffs[2]), 
                                numB = nrow(shared$df2), 
                                numDrawn = sum(shared$df1[(shared$df1$direction %in% 'enrich'),statistics[1]] <= cutoffs[1]), 
                                numWdrawn = sum(((mainresult$direction %in% 'enrich') & (mainresult$replicated)), na.rm = TRUE))
      p.dn <- .doHyperGInternal(numW = sum(shared$df2[(shared$df2$direction %in% 'deplete'),statistics[2]] <= cutoffs[2]), 
                                numB = nrow(shared$df2), 
                                numDrawn = sum(shared$df1[(shared$df1$direction %in% 'deplete'),statistics[1]] <= cutoffs[1]), 
                                numWdrawn = sum(((mainresult$direction %in% 'deplete') & (mainresult$replicated)), na.rm = TRUE))
      p.all <- .doHyperGInternal(numW = sum(shared$df2[,statistics[2]] <= cutoffs[2]), 
                                numB = nrow(shared$df2), 
                                numDrawn = sum(shared$df1[,statistics[1]] <= cutoffs[1]), 
                                numWdrawn = sum(((mainresult$replicated)), na.rm = TRUE))
      return(list('hypergeom.p.enrich' = p.up, 
                  'hypergeom.p.deplete' = p.up, 
                  'hypergeom.p.all' = p.all))
    }
    

    return(mainresult)
  }













