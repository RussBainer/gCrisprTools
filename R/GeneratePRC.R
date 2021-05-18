##' @title Generate a Precision-Recall Curve from a CRISPR screen  
##' @description Given a set of targets of interest, this function generates a Precision Recall curve from the results of 
##' a CRISPR screen. Specifically, it orders the target elements in the screen in the specified direction, and then plots the recall 
##' rate (proportion of true targets identified) against the precision (proportion of identified targets that are true targets). 
##' 
##' Note that ranking statistics in CRISPR screens are (usually) permutation-based, and so some granularity in the rankings is expected. This 
##' function does a little extra work to ensure that hits are counted as soon as the requisite value of the ranking statistic is reached 
##' regardless of where the gene is located within the block of equally-significant genes. Functionally, this means that the drawn curve is
##' somewhat anticonservative in cases where the gene ranks are not well differentiated.  
##'
##' @param summaryDF A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param target.list A character vector containing the names of the targets to be tested; by default these are assumed to be `geneID`s, 
##' but specifying `collapse=geneSymbol` enables setting on `geneSymbol` by passing that value through to `ct.simpleResult`.
##' @param direction Direction by which to order target signals (`enrich` or `deplete`).  
##' @param plot.it Logical value indicating whether to plot the curves. 
##' @return A list containing the the x and y coordinates of the curve.
##' @author Russell Bainer
##' @examples data('resultsDF')
##' data('essential.genes') #Note that this is an artificial example.
##' pr <- ct.PRC(resultsDF, essential.genes, 'enrich')
##' str(pr)
##' @export
ct.PRC <-
  function(summaryDF,
           target.list,
           direction = c("enrich", "deplete"), 
           plot.it = TRUE) {

    direction <- match.arg(direction)
    stopifnot(is(plot.it, 'logical'))
    
    collapse <- ifelse(sum(target.list %in% summaryDF$geneID) > sum(target.list %in% summaryDF$geneSymbol), 'geneID', 'geneSymbol')
    simpleDF <- ct.simpleResult(summaryDF, collapse)
    
    if(!is.character(target.list)){
      warning("Supplied target.list is not a character vector. Coercing.")
      target.list <- as.character(target.list)
    }
    present <- intersect(target.list, row.names(simpleDF))
    
    if(length(present) != length(target.list)){
      if(length(present) < 1){
        stop(paste0("None of the genes in the input list are present in the ", collapse, " column of the input data.frame."))
        }
      warning(paste(length(present), "of", length(target.list), "genes are present in the supplied results data.frame. Ignoring the remainder of the target.list."))
    }
    
    #Subset signals
    values <- simpleDF[simpleDF$direction == direction,]
    values <- c(values$best.p[order(values$best.p, decreasing =FALSE)],  rep(1, times = sum(!(simpleDF$direction %in% direction))))

    targvals <- vapply(target.list, function(x){ifelse(simpleDF[x,'direction'] %in% direction, simpleDF[x, 'best.p'], 1)}, numeric(1))
    
    out <- list()
    out$precision <- c(1, unlist(lapply(unique(values), function(x){sum(targvals <= x, na.rm = TRUE)/sum(values <= x, na.rm= TRUE)})), 0)
    out$recall <- c(0, unlist(lapply(unique(values), function(x){sum(targvals <= x, na.rm = TRUE)/length(targvals)})), 1)
    
    enrich <- switch(direction, 
                     enrich = ct.targetSetEnrichment(simpleDF, target.list, enrich = TRUE, collapse = collapse),
                     deplete = ct.targetSetEnrichment(simpleDF, target.list, enrich = FALSE, collapse = collapse)
    )
    out <- c(out, enrich)
    
    #Plot it?
    if(plot.it){
      plot(out$recall, out$precision, xlim = c(0, 1), ylim = c(0,1), 
           type = "l", ylab = "Precision", xlab = "Recall", 
           main = paste("Precision and Recall of", deparse(substitute(target.list))), col = "blue", lwd = 3)
      }
    return(invisible(out))
    }


  
  
  
  
  
  
  
  



