##' @title Identify Replicated Signals in Pooled Screens Using Conditional Scoring
##' @description This function identifies signals that are present in multiple screening experiment contrasts using a conditional 
##' strategy. Specifically, this function identifies all significant signals (according to user definitions) in a set of provided 
##' results DF and returns a `simplifiedResult` dataframe derived from the first provided contrast with an appended logical column 
##' indicating whether there is evidence for signal replication in the other provided resultsDFs. 
##' 
##' Signals are considered replicated if they cross the specified stringent threshold (default: 10% FDR) in one or more of the provided 
##' contrasts, and are similarly enriched or depleted at the relaxed threshold (default: P = 0.1) in all of the remaining contrasts. 
##' 
##' Signals are compared across screens on the basis of \code{\link{ct.regularizeContrasts}}, so users must provide an identifier 
##' with which to standardize targets (`geneID` by default). 
##' 
##' @param dflist A list of (possibly simplified) results data.frames produced by \code{\link{ct.generateResults}}. 
##' @param statistics Statistics to use to define congruence; may be a single value, but internally coerced to a vector of length 2 where the first 
##' value corresponds to the stringent cutoff annd the second value is used for the relaxed cutoff. Must be "best.p" or "best.q". 
##' @param cutoffs Numeric value(s) corresponding to the significance cutoff(s) used to define stringent and relaxed values of `statistics`. 
##' Internally coerced to a vector of length 2.
##' @param same.dir Logical vector of the same length as `dflist` indicating whether replicating signals are expected to go in the same direction 
##' (e.g., enrich/deplete in their respective screens). For example, a `dflist` of length 3 could be specified as c(TRUE, TRUE, FALSE), indicating 
##' that replicating signals should be enriched in both of the first two contrasts and depleted in the third to be considered replicated (or 
##' vise-versa). Default is `rep(TRUE, length(dflist))`.
##' @param return.stats When TRUE, return the significance of overlap instead of the logical vector (by permutation).
##' @param nperm numeric indicating number of permutations when `return.stats` is true (default 10000).  
##' @param ... Other arguments to \link{`ct.simpleResult`}, especially `collapse`.
##' @return If `return.stats` is `FALSE`, returns the first contrast as a `simplifiedResult` data.frame, with a `replicated` logical column 
##' indicating whether each signal replicates. 
##' 
##' If `return.stats` is `TRUE`, returns a dataframe indicating the hypergeometric test statistics summarizing the evidence for significantly 
##' enriched signal replication across the provided contrasts (enrich, deplete, and all together).   
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' summary(ct.compareContrasts(list(resultsDF, resultsDF[1:5000,]))$replicated)
##' ct.compareContrasts(list(resultsDF, resultsDF[1:5000,]), return.stats = TRUE)
##' @export
ct.compareContrasts  <- 
  function(dflist, 
           statistics = c('best.q', 'best.p'),
           cutoffs = c(0.1, 0.1), 
           same.dir = rep(TRUE, length(dflist)),
           return.stats = FALSE,
           nperm=10000,
           ...) {

    #Check the input: 
    dflist <- ct.regularizeContrasts(dflist, ...)
    stopifnot(all(statistics %in% c('best.p', 'best.q')), (length(statistics) <= 2), (length(statistics) > 0), is(return.stats, 'logical'),
              is(cutoffs, 'numeric'), (length(cutoffs) <= 2), (length(cutoffs) > 0), is.logical(same.dir), length(same.dir) == length(dflist))
    if(length(statistics) == 1){ 
      statistics <- rep(statistics, 2)
      }
    if(length(cutoffs) == 1){ 
      cutoffs <- rep(cutoffs, 2)
    }
 
    if(return.stats){
      stopifnot(is(nperm, 'numeric'))
    }
    
    #Find validated signals. 
    stringent <- vapply(dflist, function(x){x[,statistics[1]] <= cutoffs[1]}, logical(nrow(dflist[[1]])))
    lax <- vapply(dflist, function(x){x[,statistics[2]] <= cutoffs[2]}, logical(nrow(dflist[[1]])))
    dirs <- vapply(dflist, function(x){x$direction %in% 'enrich'}, logical(nrow(dflist[[1]])))
    dirs[,!same.dir] <- !dirs[,!same.dir]
    valid <- (rowSums(dirs) %in% c(ncol(dirs), 0)) & (rowSums(stringent) > 0) & (rowSums(lax) == ncol(lax))
    
    mainresult <- dflist[[1]]
    mainresult$replicated <- valid
    
    if(return.stats){
      #calculate a hypergeometric test P-value for enrichment/depletion
      #p.up <- gCrisprTools:::.doHyperGInternal(numW = sum(shared$df2[(shared$df2$direction %in% 'enrich'),statistics[2]] <= cutoffs[2]), 
      #                          numB = nrow(shared$df2), 
      #                          numDrawn = sum(shared$df1[(shared$df1$direction %in% 'enrich'),statistics[1]] <= cutoffs[1]), 
      #                          numWdrawn = sum(((mainresult$direction %in% 'enrich') & (mainresult$replicated)), na.rm = TRUE))
      #p.dn <- gCrisprTools:::.doHyperGInternal(numW = sum(shared$df2[(shared$df2$direction %in% 'deplete'),statistics[2]] <= cutoffs[2]), 
      #                          numB = nrow(shared$df2), 
      #                          numDrawn = sum(shared$df1[(shared$df1$direction %in% 'deplete'),statistics[1]] <= cutoffs[1]), 
      #                          numWdrawn = sum(((mainresult$direction %in% 'deplete') & (mainresult$replicated)), na.rm = TRUE))
      #p.all <- gCrisprTools:::.doHyperGInternal(numW = sum(shared$df2[,statistics[2]] <= cutoffs[2]), 
      #                          numB = nrow(shared$df2), 
      #                          numDrawn = sum(shared$df1[,statistics[1]] <= cutoffs[1]), 
      #                          numWdrawn = sum(((mainresult$replicated)), na.rm = TRUE))
      
      obs <- c(sum((rowSums(dirs) == ncol(dirs)) & valid), sum((rowSums(dirs) == 0) & valid), sum(valid))
      
      perm <- t(replicate(nperm, 
                        expr={
                          pr <- lapply(1:ncol(stringent), function(x){sample(1:nrow(stringent))})
                          nr <- nrow(stringent)

                          p.stringent <- vapply(1:length(pr), function(x){stringent[pr[[x]], x]}, logical(nr))
                          p.lax <- vapply(1:length(pr), function(x){lax[pr[[x]], x]}, logical(nr))
                          p.dirs <- vapply(1:length(pr), function(x){dirs[pr[[x]], x]}, logical(nr))
                          
                          dn <- sum((rowSums(p.dirs) == 0) & (rowSums(p.stringent) > 0) & (rowSums(p.lax) == ncol(p.lax))) 
                          up <- sum((rowSums(p.dirs) == ncol(p.dirs)) & (rowSums(p.stringent) > 0) & (rowSums(p.lax) == ncol(p.lax)))
                          return(c(up, dn, up+dn))
                        }, 
                        simplify = TRUE))
      out <- data.frame('expected' = colMeans(perm), 
                        'observed' = obs, 
                        'p' = c(sum(perm[,1] > sum((rowSums(dirs) == 0) & valid))/nperm, 
                                sum(perm[,2] > sum((rowSums(dirs) == ncol(dirs)) & valid))/nperm, 
                                sum(perm[,3] > sum(valid))/nperm)) 
      
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
ct.upSet <- function(dflist, 
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










