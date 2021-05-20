##' @title Generate a Figure Summarizing Overall Signal for One or More Targets
##' @description Given one or more targets of interest, this function generates a summary image contextualizing the 
##' corresponding signals within the provided contrast. This takes the form of an annotated ranking 
##' curve of target-level signals, supplemented with horizontal Q-value cutoffs and an inset volcano plot of gRNA 
##' behavior. 
##' 
##' Limited annotation is provided for the specified targets using the following logic: 
##' 
##' - If a character vector is provided, up to five targets are annotated; longer lists are highlighted without specifying individual elements.  
##' - If a list is provided, the `names` element is used as the annotation. This is similarly constrained to a total of 5 annotated elements. 
##' 
##' @param summaryDF A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param targets A list or character vector containing the names of the targets to be displayed. Only targets contained in the column 
##' specified by the `collapse` parameter to `ct.simpleResult()` will be displayed; default is `geneSymbol`. Plotting priority (e.g., 
##' the points to plot last in the case of overlapping signals) is given to earlier elements in the list. 
##' @param callout Logical indicating whether lines should be plotted indicating individual gene sets to augment the point highlighting.
##' @param ... Additional optional arguments to `ct.simpleResult()`
##' @return A summary plot on the current device. 
##' @author Russell Bainer
##' @examples data('resultsDF')
##' ct.signalSummary(resultsDF, list('CandidateA' = 'Target229', 'Pathway3' = resultsDF$geneSymbol[c(42,116,1138,5508)]))
##' @export
ct.signalSummary <-
  function(summaryDF,
           targets,
           callout = FALSE, 
           ...) {

    #Check the input: 
    ct.resultCheck(summaryDF)
    simpleDF <- ct.simpleResult(summaryDF, ...)

    if(!is.list(targets)){
      targets <- as.list(targets)
      names(targets) <- targets
    }
    
    if(length(targets) > 5){
      stop('Too many targets specified; suppressing annotation.')
      targets <- list(unique(unlist(targets)))
    }
    targets <- rev(targets)
    
    bad <- setdiff(unlist(targets), row.names(simpleDF))
    if(length(bad) > 0){
      stop(paste0('Cannot find the following supplied targets in the geneSymbol column of the supplied DF: ', paste(bad, collapse = ', ')))
    }
    id <- ifelse(unlist(targets)[1] %in% simpleDF$geneID, 'geneID', 'geneSymbol')
    
    #Prep data
    top <- simpleDF[simpleDF$direction == 'enrich',]
    top <- top[order(top$best.p, top$Rho_enrich, decreasing = FALSE),]
    bot <- simpleDF[simpleDF$direction == 'deplete',]
    bot <- bot[order(bot$best.p, bot$Rho_deplete, decreasing = TRUE),]
    simpleDF <- rbind(top, bot)

    gwp <- ct.softLog(simpleDF$best.p)
    gwp <- gwp*vapply(simpleDF$direction, function(x){ifelse(x == 'enrich', 1, -1)}, numeric(1))
    names(gwp) <- row.names(simpleDF)

    exes <- (1:length(gwp))/length(gwp)
    
    #qcut <- c(min(gwp[(simpleDF$direction == 'enrich') & (simpleDF$best.q <= 0.1)]), min(gwp[(simpleDF$direction == 'deplete') & (simpleDF$best.q <= 0.1)]))

    
    #Compose Plot
    plot(exes, gwp, 
         ylab = 'Target -log10 P', xaxt = 'n', xlab = 'Signal Rank, Most Enriched to Most Depleted', 
         pch = 19, cex = 0.5, col = rgb(14/255,41/255,56/255))
    #inset
    maxval <- max(gwp)
    glfc <- summaryDF$`gRNA Log2 Fold Change`
    gp <- rowMax(cbind(-log10(summaryDF$`gRNA Enrichment P`), -log10(summaryDF$`gRNA Depletion P`)))
    
    
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
                 selected <- intersect(row.names(simpleDF),targets[[x]])
                 gw.ranks <- vapply(selected, 
                                    #function(x){grep(x, summaryDF$geneSymbol[genewise], fixed = TRUE)},
                                    function(y){which(row.names(simpleDF) %in% y)},
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
                       picked <- intersect(row.names(simpleDF),targets[[x]])
                       gw.ranks <- vapply(picked, 
                                          #function(x){grep(x, summaryDF$geneSymbol[genewise], fixed = TRUE)},
                                          function(y){which(row.names(simpleDF) %in% y)},
                                          integer(1))
                       reagents <- which(summaryDF[,id] %in% picked)
                       points(inset.x[reagents], inset.y[reagents], col = rgb(14/255,41/255,56/255), bg = t.col[x], pch = 21, cex = 0.7, lwd = 1.2)
                       points(exes[gw.ranks], gwp[gw.ranks], pch = 21, bg = t.col[x], cex = 1.2, lwd = 1.2)
                     }))
  }


##' @title Visualize Signal Across A List of Contrasts 
##' @description Given a list of provided results `data.frame`s summarizing a series of contrasts from one or more pooled screens, 
##' this function visualizes their respective signals as a series of stacked barcharts. Enriched signals 
##' are represented in the positive direction, and depleted signals are represented in the negative direction.  Note that the 
##' provided contrast results are not regularized by this function.
##' 
##' This function may be used to compare signals across different screen contrasts, or to compare signals within interesting 
##' subsets of targets ascertained within a single experiment.  
##' 
##' @param dflist A named list of `data.frame`s summarizing the results of one or more screen contrasts, returned by the function 
##' \code{\link{ct.generateResults}}. 
##' @param background Logical indicating whether to represent the nonsignificant hits in the barchart.  
##' @param statistic Should cutoffs be calculated based on FDR (`best.q`) or P-value (`best.p`)?
##' @param ... Other parameters to lower functions, especially `ct.simpleResult()`
##' @return A summary plot on the current device. Invisibly, the data.frame tallying signals at various thresholds. 
##' @author Russell Bainer
##' @examples data('resultsDF')
##' ct.contrastBarchart(list('FirstResult' = resultsDF, 'SecondResult' = resultsDF))
##' ct.contrastBarchart(list('FirstResult' = resultsDF, 'SecondResult' = resultsDF), background = FALSE)
##' ct.contrastBarchart(list('FirstResult' = resultsDF[1:1000,], 'SecondResult' = resultsDF))
##' @export
ct.contrastBarchart <- function(dflist, background = TRUE, statistic = c('best.q', 'best.p'), ...){

  #Check input
  stopifnot(is(background, 'logical'))
  statistic <- match.arg(statistic)

  #listify results as needed
  if((!is(dflist, 'list'))){
    if(ct.resultCheck(dflist)){
      dflist <- list('result' = dflist)
    }
  }
  dflist <- sapply(dflist, 
                   function(x){ct.simpleResult(x, ...)}, 
                   simplify =  FALSE)
  
  colors <- c('grey', colorRampPalette(c('white', 'orange', 'red', 'darkred'))(5)[2:5])
  names(colors) <- c('N/S', '< 0.1', '< 0.01', '< 0.001', '< 0.00001')
  
  #Collect values
  vals <- vapply(dflist, 
                 function(x){
                   c(sum((x$direction == 'enrich')), 
                     sum((x[,statistic] <= 0.1) & (x$direction == 'enrich')),
                     sum((x[,statistic] <= 0.01) & (x$direction == 'enrich')), 
                     sum((x[,statistic] <= 0.001) & (x$direction == 'enrich')), 
                     sum((x[,statistic] <= 0.00001) & (x$direction == 'enrich')), 
                     -sum((x[,statistic] <= 0.00001) & (x$direction == 'deplete')), 
                     -sum((x[,statistic] <= 0.001) & (x$direction == 'deplete')), 
                     -sum((x[,statistic] <= 0.01) & (x$direction == 'deplete')),
                     -sum((x[,statistic] <= 0.1) & (x$direction == 'deplete')),
                     -sum((x$direction == 'deplete'))) 
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
  





