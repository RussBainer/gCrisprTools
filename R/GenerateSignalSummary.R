##' @title Generate a Figure Summarizing Overall Signal for One or More Targets
##' @description Given one or more targets of interest, this function generates a summary image contextualizing the 
##' corresponding signals within the contest of the provided contrast. This takes the form of an annotated ranking 
##' curve of target-level signals, supplemented with horizontal Q-value cutoffs and an inset volcano plot of gRNA 
##' behavior. 
##' 
##' Limited annotation is provided for the specified targets using the following logic: 
##' 
##' - If a character vector is provided, up to five targets are annotated; longer lists are highlighted without specifying individual elements.  
##' - If a list is provided, the `names` element is used as the annotation. This is similarly constrained to a total of 5 annotated elements. 
##' 
##' @param summaryDF A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param targets A list or character vector containing the names of the targets to be displayed. Only targets contained in the \code{geneSymbol} 
##' column of the provided \code{summaryDF} are considered. Plotting priority (e.g., the points to plot last in the case 
##' of overlapping signals) is given to earlier elements in the list. 
##' @param direction Should enrichment or depletion be considered? Must be one of \code{"enrich"} or \code{"deplete"}.
##' @param callout Logical indicating whether lines should be plotted indicating individual gene sets to augment the point highlighting.
##' @return A summary plot on the current device. 
##' @author Russell Bainer
##' @examples data('resultsDF')
##' ct.signalSummary(resultsDF, list('CandidateA' = 'Target229', 'Pathway3' = resultsDF$geneSymbol[c(42,116,1138,5508)]), 'enrich')
##' @export
ct.signalSummary <-
  function(summaryDF,
           targets,
           direction = c("enrich", "deplete"), 
           callout = FALSE) {

    #Check the input: 
    direction <- match.arg(direction)
    if(!ct.resultCheck(summaryDF)){
      stop("Execution halted.")
    }

    if(!is.list(targets)){
      targets <- as.list(targets)
      names(targets) <- targets
    }
    
    if(length(targets) > 5){
      stop('Too many targets specified; suppressing annotation.')
      targets <- list(unique(unlist(targets)))
    }
    targets <- rev(targets)
    
    bad <- setdiff(unlist(targets), summaryDF$geneSymbol)
    if(length(bad) > 0){
      stop(paste0('Cannot find the following supplied targets in the geneSymbol column of the supplied DF: ', paste(bad, collapse = ', ')))
    }
    
    #Prep data
    summaryDF <- switch(direction, 
                        'enrich' = summaryDF[order(summaryDF$`Target-level Enrichment P`, decreasing = FALSE),], 
                        'deplete' = summaryDF[order(summaryDF$`Target-level Depletion P`, decreasing = FALSE),])
    
    p <- switch(direction, 
                   'enrich' = -log10(summaryDF$`Target-level Enrichment P`), 
                   'deplete' = -log10(summaryDF$`Target-level Depletion P`))
    p[is.infinite(p)] <- max(p[is.finite(p)]) + 0.5

    genewise <- !duplicated(summaryDF$geneSymbol)
    gwp <- p[genewise]
    names(gwp) <- summaryDF$geneSymbol[genewise]
    exes <- (1:length(gwp))/length(gwp)
    
    qcut <- min(p[switch(direction, 
                              'enrich' = summaryDF$`Target-level Enrichment Q`[genewise] < 0.1, 
                              'deplete' = summaryDF$`Target-level Depletion Q`[genewise] < 0.1)])
    
    #Compose Plot
    plot(exes, gwp, 
         ylab = 'Target -log10 P', xaxt = 'n', xlab = 'Signal Rank', 
         pch = 19, cex = 0.5, col = rgb(14/255,41/255,56/255))
    #inset
    maxval <- max(gwp)
    glfc <- summaryDF$`gRNA Log2 Fold Change`
    gp <- switch(direction, 
                 'enrich' = -log10(summaryDF$`gRNA Enrichment P`), 
                 'deplete' = -log10(summaryDF$`gRNA Depletion P`))
    inset.x <- ((0.29/(max(glfc) - min(glfc))) * (glfc - min(glfc))) + 0.7
    inset.y <- ((((maxval - 0.1) - ((maxval)/2))/(max(gp) - min(gp))) * (gp - min(gp))) + ((maxval + 0.05)/2)
    inset.zero <- inset.x[which(abs(glfc) == min(abs(glfc), na.rm = TRUE))][1]
    
    polygon(x = c(0.7,1,1,0.7), y = (rep(maxval, 4) - rep(c(maxval/2, 0), each = 2)), col = 'white') 
    lines(rep(inset.zero, 2), c(maxval/2, maxval), lty = 2, col = 'lightgrey')
    points(inset.x, inset.y, pch = 19, cex = 0.2, col= rgb(14/255,41/255,56/255))
    graphics::text(0.85, maxval, 'gRNA', adj = c(0.5,1.5), cex = 1)
    graphics::text(0.85, (maxval/2), 'Log2 Fold Change', adj = c(0.5, 1.5), cex = 0.7)
    graphics::text(0.7, 3*(maxval/4), '-log10P', srt = 90, adj = c(0.7, -0.5), cex = 0.7)
    
    #add annotation
    t.col <- c(rgb(218/255, 111/255, 90/255), 
               'white', 
               rgb(111/255, 130/255, 138/255), 
               rgb(100/255, 190/255, 203/255),
               rgb(127/255, 47/255, 105/255))
   
    #Optionally add annotations
    if(length(targets) <= 5){ 
      ylocs <- rev(seq((maxval * 0.95), (maxval/2), length.out = 5)[1:length(targets)]) 
      ybuff <- (maxval/20)/2
      
      invisible(
        lapply(1:length(targets),
               function(x){
                 selected <- (summaryDF$geneSymbol %in% targets[[x]])
                 all.targ <- (selected & genewise) 
                 gw.ranks <- vapply(summaryDF$geneSymbol[all.targ], 
                                    #function(x){grep(x, summaryDF$geneSymbol[genewise], fixed = TRUE)},
                                    function(x){which(summaryDF$geneSymbol[genewise] == x)},
                                    integer(1))
                 
                 if(callout){
                   segments(exes[gw.ranks], gwp[gw.ranks], 0.3, ylocs[x], col = rgb(78/255, 78/255, 76/255, 0.4))
                 }
                 graphics::text(0.3, ylocs[x], names(targets)[x], pos = 4)
                 points(0.3, ylocs[x], pch = 22, cex = 2, bg = t.col[x])
               })
      )
    }
    
    #Highlight Points
    invisible(lapply(1:length(targets), 
                     function(x){
                       selected <- (summaryDF$geneSymbol %in% targets[[x]])
                       all.targ <- (selected & genewise) 
                       gw.ranks <- vapply(summaryDF$geneSymbol[all.targ], 
                                          #function(x){grep(x, summaryDF$geneSymbol[genewise], fixed = TRUE)},
                                          function(x){which(summaryDF$geneSymbol[genewise] == x)},
                                          integer(1))
                       points(inset.x[selected], inset.y[selected], col = rgb(14/255,41/255,56/255), bg = t.col[x], pch = 21, cex = 0.7, lwd = 1.2)
                       points(exes[gw.ranks], gwp[gw.ranks], pch = 21, bg = t.col[x], cex = 1.2, lwd = 1.2)
                     }))
  }


##' @title Visualize Signal Across A List of Contrasts 
##' @description Given a list of provided results `data.frame`s summarizing a series of contrasts from one or more pooled screens, 
##' this function visualizes the signal associated with their shared targets as a series of stacked barcharts. Enriched signals 
##' are represented in the positive direction, and depleted signals are represented in the negative direction. 
##' 
##' @param dflist A named list of `data.frame`s summarizing the results of one or more screen contrasts, returned by the function 
##' \code{\link{ct.generateResults}}. 
##' @param background Logical indicating whether to represent the nonsignificant hits in the barchart.  
##' @param ... Additional parameters for internal functions, such as `ct.simpleResult()`
##' @param statistic Should cutoffs be calculated based on FDR (`best.q`) or P-value (`best.p`)?
##' @param ... Other parameters to lower functions, especially `ct.simpleResult()`
##' @return A summary plot on the current device. Invisibly, the data.frame tallying signals at various thresholds. 
##' @author Russell Bainer
##' @examples data('resultsDF')
##' # Not so interesting b/c of limited weak signal in example
##' ct.contrastBarchart(list('First Result' = resultsDF, 'Second Result' = resultsDF))
##' @export
ct.contrastBarchart <- function(dflist, background = TRUE, statistic = c('best.q', 'best.p'), ...){

  #Check input
  dflist <- ct.regularizeContrasts(dflist, ...)
  stopifnot(is(background, 'logical'))
  statistic <- match.arg(statistic)
  
  colors <- c('grey', colorRampPalette(c('white', 'orange', 'red', 'darkred'))(5)[2:5])
  names(colors) <- c('N/S', '< 0.1', '< 0.01', '< 0.001', '< 0.00001')
  
  #Collect values
  vals <- vapply(dflist, 
                 function(x){
                   c(sum((x[,statistic] > 0.1) & (x$direction == 'enrich')), 
                     sum((x[,statistic] > 0.01) & (x[,statistic] <= 0.1) & (x$direction == 'enrich')),
                     sum((x[,statistic] > 0.001) & (x[,statistic] <= 0.01) & (x$direction == 'enrich')), 
                     sum((x[,statistic] > 0.00001) & (x[,statistic] <= 0.001) & (x$direction == 'enrich')), 
                     sum((x[,statistic] < 0.00001) & (x$direction == 'enrich')), 
                     -sum((x[,statistic] < 0.00001) & (x$direction == 'deplete')), 
                     -sum((x[,statistic] > 0.00001) & (x[,statistic] <= 0.001) & (x$direction == 'deplete')), 
                     -sum((x[,statistic] > 0.001) & (x[,statistic] <= 0.01) & (x$direction == 'deplete')),
                     -sum((x[,statistic] > 0.01) & (x[,statistic] <= 0.1) & (x$direction == 'deplete')),
                     -sum((x[,statistic] > 0.1) & (x$direction == 'deplete'))) 
                 }, numeric(10))
  
  row.names(vals) <- paste0(rep(c('enrich', 'deplete'), each = 5), 
                            '_p', 
                            c('Other', '1', '01', '001', '00001', '00001', '001', '01', '1', 'Other'))
  
  if(!background){
    vals <- vals[2:9,]
    colors <- colors[2:5]
  }
  
  plot(NA, xlim = range(vals), ylim = c(0,(length(dflist))), 
       yaxt = "n",  
       ylab = '', xlab = 'Significant Target Signals', 
       main = switch(statistic, 'best.p' = 'Target P Values', 'best.q' = 'Target Q Values'))

  for(j in 1:length(dflist)){
    for(k in 1:5){
     polygon(x = c(vals[k,j], vals[k,j], vals[(nrow(vals) - k + 1),j], vals[(nrow(vals) - k + 1),j]), 
             y = c((j-0.8), (j-0.2), (j-0.2), (j - 0.8)), col = colors[k])
    }
    text(range(vals)[1], j, names(dflist)[j], adj = c(0,1))
  }
  abline(v = 0, lty = 2, col = 'black')
  legend('bottom', names(colors), fill = colors, horiz = TRUE, cex = 0.7)
  text(range(vals)[1], 0, 'Depleted', col = 'grey', adj = c(0,0))
  text(range(vals)[2], 0, 'Enriched', col = 'grey', adj = c(1,0))
  
  return(invisible(vals))
}
  





