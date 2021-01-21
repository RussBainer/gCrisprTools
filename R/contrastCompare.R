##' @title Compare Two CRISPR Screen Contrasts via a Scatter Plot
##' @description This is a function for comparing the results of two screening experiments. Given two \code{summaryDF}, 
##' the function places them in register with one another, generates a simplified scatter plot where enrichment or depletion 
##' in each contrast is represented by the associated -log10 *P*-value, and returns an 
##' invisible list of targets present in the four significance quadrants. 
##' 
##' This is a target-level analysis, and some minor simplifications are introduced to screen signals for the sake of clarity. 
##' Principal among these is the decision to collapse gene signals to a single directional enrichment statistic; as target-level
##' signals aretypically aggregates of many guide-level signals, it is formally possible for targets to be both significantly 
##' enriched and significantly depleted within a single screen contrast as a result of substntially divergent reagent activity. 
##' This behavior is uncommon, however, and so targets are represented by selecting the direction of enrichment or depletion 
##' associated with the most significant *P*-value. 
##'
##' @param df1 A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param df2 A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param targets Column of the provided \code{summaryDF} to consider. Must be \code{geneID} or \code{geneSymbol}.
##' @param statistic Statistic to plot on each axis (after -log10 transformation). Must be "p", "q", or "rho".
##' @param cutoff significance cutoff used to define the significance quadrants (cannot be exactly zero).
##' @param plot.it Logical indicating whether to compose the plot on the default device. 
##' @return Invisibly, a list of length 4 containing the genes passing significance for the respective quadrants.
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' cat <- ct.scatter(resultsDF[,100:2100], resultsDF[1:2000,])
##' head(cat)
##' @export
ct.scatter <- 
  function(df1, df2,
           targets = c('geneSymbol', 'geneID'),
           statistic = c('best.p', 'best.q'), 
           cutoff = 0.05,
           plot.it = TRUE) {

    #Check the input: 
    targets <- match.arg(targets)
    statistic <- match.arg(statistic)
    stopifnot(is.numeric(cutoff), is.logical(plot.it), cutoff > 0)
    cutoff <- -log10(cutoff)
    
    dfs <- ct.regularizeContrasts(df1, df2, collapse = targets)
    dfs <- sapply(dfs, 
                  function(x){
                    x[,c("Rho_enrich", "Rho_deplete", "best.p", "best.q")] <- apply(x[,c("Rho_enrich", "Rho_deplete", "best.p", "best.q")], 2, ct.softLog)
                    return(x)
                  }, simplify = FALSE)

    #Form output & divide into quadrants: 
    out <- cbind(dfs$df1[,1:6], dfs$df2[,3:6])
    colnames(out) <- c('geneID', 'geneSymbol', 
                       'Rho_enrich.c1', 'Rho.deplete.c1', 'log10p.c1', 'log10q.c1',
                       'Rho_enrich.c2', 'Rho.deplete.c2', 'log10p.c2', 'log10q.c2')
    c1.sig <- (dfs$df1[,statistic] >= cutoff)
    c2.sig <- (dfs$df2[,statistic] >= cutoff)
    
    out$quadrant <- 5
    out$quadrant[(c1.sig & c2.sig & (dfs$df1$direction == 'deplete') & (dfs$df2$direction == 'enrich'))] <- 1
    out$quadrant[(!c1.sig & c2.sig & (dfs$df2$direction == 'enrich'))] <- 2
    out$quadrant[(c1.sig & c2.sig & (dfs$df1$direction == 'enrich') & (dfs$df2$direction == 'enrich'))] <- 3
    out$quadrant[(c1.sig & !c2.sig & (dfs$df1$direction == 'deplete'))] <- 4
    out$quadrant[(c1.sig & !c2.sig & (dfs$df1$direction == 'enrich'))] <- 6
    out$quadrant[(c1.sig & c2.sig & (dfs$df1$direction == 'deplete') & (dfs$df2$direction == 'deplete'))] <- 7
    out$quadrant[(!c1.sig & c2.sig & (dfs$df2$direction == 'deplete'))] <- 8
    out$quadrant[(c1.sig & c2.sig & (dfs$df1$direction == 'enrich') & (dfs$df2$direction == 'deplete'))] <- 9
    
    #Plot it
    if(plot.it){
      
      plot(dfs$df1[,statistic] * vapply(dfs$df1$direction, function(x){ifelse(x == 'enrich', 1, -1)}, numeric(1)),  
           dfs$df2[,statistic] * vapply(dfs$df2$direction, function(x){ifelse(x == 'enrich', 1, -1)}, numeric(1)),  
           main = paste0(names(dfs)[1], ' vs ', names(dfs)[2], ' (', statistic, ')'), 
           xlab = paste0(names(dfs)[1], ' -log10(', statistic, ')'), 
           ylab = paste0(names(dfs)[2], ' -log10(', statistic, ')'), 
           pch = 19, col = rgb(0,0,0, 0.4))
      abline(v = c(-cutoff, cutoff), h = c(-cutoff, cutoff), lty = 2)
      
      points((dfs$df1[,statistic] * vapply(dfs$df1$direction, function(x){ifelse(x == 'enrich', 1, -1)}, numeric(1)))[out$quadrant %in% c(1,3,7,9)],  
             (dfs$df2[,statistic] * vapply(dfs$df2$direction, function(x){ifelse(x == 'enrich', 1, -1)}, numeric(1)))[out$quadrant %in% c(1,3,7,9)],  
             pch = 19, col = rgb(0.5,0,0, 0.4))
      points((dfs$df1[,statistic] * vapply(dfs$df1$direction, function(x){ifelse(x == 'enrich', 1, -1)}, numeric(1)))[out$quadrant %in% c(2,4,6,8)],  
             (dfs$df2[,statistic] * vapply(dfs$df2$direction, function(x){ifelse(x == 'enrich', 1, -1)}, numeric(1)))[out$quadrant %in% c(2,4,6,8)],  
             pch = 19, col = rgb(0,0,0.6, 0.4))
    }
    return(invisible(out))
  }













