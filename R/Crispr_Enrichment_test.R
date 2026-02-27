##' @title Test Whether a Specified Target Set is Enriched Within a Pooled Screen
##' @description This function takes in a \code{resultsDF} and a vector of \code{targets} (contained in the \code{geneID} column of
##' \code{resultsDF}) and determines whether the specified targets are enriched within the set of all significantly altered targets.
##' It does this by iteratively testing whether \code{targets} are more likely to be among the set of enriched or depleted targets
##' at various significance thresholds using a hypergeometric test. Note that the returned Hypergeometric P-values are not corrected
##' for multiple testing.
##'
##' Returns a list detailing the \code{targets} used in the tests, and tables indicating the results of the hypergeometric test
##' at various significance thresholds.
##' @param summaryDF A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' Internally coerced via `ct.simpleResult()`.
##' @param targets A character vector containing the names of the targets to be tested; by default these are assumed to be `geneID`s, 
##' but specifying `collapse=geneSymbol` enables setting on `geneSymbol` by passing that value through to `ct.simpleResult`.
##' @param enrich Logical indicating whether to consider guides that are enriched (default) or depleted within the screen.
##' @param ignore Optionally, a character vector containing elements of the \code{summaryDF} that should be ignored in the analysis 
##' (e.g., unassignable or nonfunctional targets, such as nontargeting controls). By default, this function omits targets with 
##' \code{geneSymbol} 'NoTarget'.
##' @param ... Additional parameters to pass to `ct.simpleResult`.
##' @return A named list containing the tested target set and tables detailing the hypergeometric test results using various P-value and
##' Q-value thresholds.
##' @examples data(resultsDF)
##' tar <-  sample(unique(resultsDF$geneSymbol), 20)
##' res <- ct.targetSetEnrichment(resultsDF, tar)
##' @author Russell Bainer
##' @export
ct.targetSetEnrichment <- function(summaryDF, targets, enrich = TRUE, ignore = 'NoTarget', ...){

  #input checks
  stopifnot(is(enrich, 'logical'), is(ignore, 'character'))
  
  if(!is.character(targets)){
    warning("Supplied targets are not a character vector. Coercing.")
    targets <- as.character(targets)
  }
  
  #Infer whether Gsdb is ID or feature centric
  gids <- sum(targets %in% summaryDF$geneID)
  gsids <- sum(targets %in% summaryDF$geneSymbol)
  
  if(all(c(gsids, gids) == 0)){
    stop('None of the targets are present in either the geneID or geneSymbol slots of the first provided result.')
  }
  
  collapse <- ifelse(gids > gsids, 'geneID', 'geneSymbol')
  summaryDF <- ct.simpleResult(summaryDF, collapse)

  valid <- intersect(targets, row.names(summaryDF))

  if(length(setdiff(targets, valid)) != 0){
    warning('Not all of the supplied targets are present in the summary dataframe. Proceeding with ',
            length(valid), ' targets.')
  }

  #Condense the summary frame to gene-level estimates and isolate the ones that we are testing
  summaryDF <- summaryDF[setdiff(row.names(summaryDF), ignore),]  #Remove NoTarget Genes rather than constraining to Entrez

  #Pull out the P-values
  selected <- summaryDF[,c("best.p", "best.q")]
  
  if(enrich){
    selected[(summaryDF$direction %in% 'deplete'), ] <- c(1,1)
  } else {
    selected[(summaryDF$direction %in% 'enrich'), ] <- c(1,1)
  }

  #Set the cutoffs
  cuts <- c(0,1/(10^(5:1)), 0.5, 1)

  out <- list('targets' = valid)
  out$P.values <- cbind(cuts,
                        vapply(cuts, function(x){sum(selected[valid,1] <= x)}, numeric(1)),
                        vapply(cuts, function(x){
                                            .doHyperGInternal(length(valid),
                                                nrow(selected),
                                                sum(selected[,1] <= x),
                                                sum(selected[valid,1] <= x),
                                                TRUE)$p},
                                      numeric(1))

                          )
  out$Q.values <- cbind(cuts,
                        vapply(cuts, function(x){sum(selected[valid,2] <= x)}, numeric(1)),
                        vapply(cuts, function(x){
                          .doHyperGInternal(length(valid),
                                            nrow(selected),
                                            sum(selected[,2] <= x),
                                            sum(selected[valid,2] <= x),
                                            TRUE)$p},
                          numeric(1))
                        )


  colnames(out$P.values) <- colnames(out$Q.values) <- c('Cutoff', 'Sig', 'Hypergeometric_P')
  return(out)
}






## -----------------------------------------------------------------------------
## These were copied from the Bioconductor collection package

##' We envision the test as follows:
##'
##' The urn contains genes from the gene universe.  Genes annotated at a
##' given collection term are white and the rest black.
##'
##' The number drawn is the size of the selected gene list.  The
##' number of white drawn is the size of the intersection of the
##' selected list and the genes annotated at the collection.
##' Here's a diagram based on using GO as the collection:
##'
##'          inGO    notGO
##'          White   Black
##' selected  n11     n12
##' not       n21     n22
##'
##' numW: number of genes in GO category
##' numB: size of universe
##' numDrawn: number of differentially expressed genes
##' numWdrawn: the number of genes differentially expressed in category
##' @keywords internal
##' @noRd

.doHyperGInternal <- function(numW, numB, numDrawn, numWdrawn, over = TRUE) {
  n21 <- numW - numWdrawn
  n12 <- numDrawn - numWdrawn
  n22 <- numB - n12

  odds_ratio <-  (numWdrawn * n22) / (n12 * n21)

  expected <- (numWdrawn + n12) * (numWdrawn + n21)
  expected <- expected / (numWdrawn + n12 + n21 + n22)

  if (over) {
    ## take the -1 because we want evidence for as extreme or more
    pvals <- phyper(numWdrawn - 1L, numW, numB,
                    numDrawn, lower.tail=FALSE)
  } else {
    pvals <- phyper(numWdrawn, numW, numB,
                    numDrawn, lower.tail=TRUE)
  }

  list(p=pvals, odds=odds_ratio, expected=expected)
}
