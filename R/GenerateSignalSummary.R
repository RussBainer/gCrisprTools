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
##' @param targets A list or character vector containing the names of the targets to be displayed. 
##' 
##' Only targets contained in the \code{geneSymbol} 
##' column of the provided \code{summaryDF} are considered.
##' @param direction Should enrichment or depletion be considered? Must be one of \code{"enrich"} or \code{"deplete"}.
##' @return A summary plot on the current device. 
##' @author Russell Bainer
##' @examples data('resultsDF')
##' data('essential.genes') #Note that this is an artificial example.
##' pr <- ct.PRC(resultsDF, essential.genes, 'enrich.p')
##' str(pr)
##' @export

ct.signalSummary <-
  function(summaryDF,
           targets,
           direction = c("enrich", "deplete")) {

    #Check the input: 
    direction <- 
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
    
    bad <- setdiff(unlist(targets), summaryDF$geneSymbol)
    if(length(bad) > 0){
      stop(paste0('Cannot find the following supplied targets in the geneSymbol column of the supplied DF: ', paste(bad, collapse = ', ')))
    }
    
    #Prep data
    summaryDF <- switch(match.arg(direction), 
                        'enrich' = summaryDF[order(summaryDF$`Target-level Enrichment P`, decreasing = FALSE),], 
                        'deplete' = summaryDF[order(summaryDF$`Target-level Depletion P`, decreasing = FALSE),])
    
    p <- switch(match.arg(direction), 
                   'enrich' = -log10(summaryDF$`Target-level Enrichment P`), 
                   'deplete' = -log10(summaryDF$`Target-level Depletion P`))
    p[is.infinite(p)] <- max(p[is.finite(p)]) + 0.5

    genewise <- !duplicated(summaryDF$geneSymbol)
    gwp <- p[genewise]
    exes <- order(gwp, decreasing = TRUE)/length(gwp)
    
    qcut <- min(p[switch(match.arg(direction), 
                              'enrich' = summaryDF$`Target-level Enrichment Q`[genewise] < 0.1, 
                              'deplete' = sumamryDF$`Target-level Depletion Q`[genewise] < 0.1)])
    
    #Compose Plot
    plot(exes, gwp, 
         ylab = '-log10 P', xaxt = 'n', xlab = 'Signal Rank', 
         pch = 19, cex = 0.5, col = rgb(14/255,41/255,56/255))
    #inset
    maxval <- max(gwp)
    glfc <- summaryDF$`gRNA Log2 Fold Change`
    gp <- switch(match.arg(direction), 
                 'enrich' = -log10(summaryDF$`gRNA Enrichment P`), 
                 'deplete' = -log10(summaryDF$`gRNA Enrichment P`))
    g.big <-  switch(match.arg(direction), 
                     'enrich' = glfc > 0.1, 
                     'deplete' = glfc < 0.1)
    inset.x <- ((0.29/(max(glfc))) * (glfc)) + 0.7
    inset.y <- ((((maxval - 0.05) - (maxval/2))/(max(gp) - min(gp))) * (gp - min(gp))) + (maxval/2)
    
    polygon(x = c(0.7,1,1,0.7), y = (rep(maxval, 4) - rep(c(maxval/2, 0), each = 2)), col = 'white') 
    points(inset.x[g.big], inset.y[g.big], pch = 19, cex = 0.2, col= rgb(14/255,41/255,56/255))
    text(0.7, maxval, 'gRNA', adj = c(-0.2,1.5), cex = 0.5)
    text(0.85, (maxval/2), 'Log2 Fold Change', pos = 1, cex = 0.7)
    text(0.7, 3*(maxval/4), '-log10P', srt = 90, adj = c(0.7, -0.5), cex = 0.7)
    
    #add annotation
    

    
        
    
    #Convert to gene-level stats
    summaryDF <- summaryDF[!duplicated(summaryDF$geneID),]
    row.names(summaryDF) <- summaryDF$geneID
    
    if(!is.character(target.list)){
      warning("Supplied target.list is not a character vector. Coercing.")
      target.list <- as.character(target.list)
    }
    present <- intersect(target.list, summaryDF$geneID)
    if(length(present) != length(target.list)){
      if(length(present) < 1){
        stop("None of the genes in the input list are present in the geneSymbol column of the input data.frame.")
        }
      warning(paste(length(present), "of", length(target.list), "genes are present in the supplied results data.frame. Ignoring the remainder of the target.list."))
    }
    
    #Gather the values for the targets: 
    stat <- match.arg(stat)
    targvals <- switch(stat, 
         enrich.p = (summaryDF[(summaryDF$geneID %in% present),"Target-level Enrichment P"]), 
         deplete.p = (summaryDF[(summaryDF$geneID %in% present),"Target-level Depletion P"]), 
         enrich.fc = (-summaryDF[(summaryDF$geneID %in% present),"Median log2 Fold Change"]), 
         deplete.fc = (summaryDF[(summaryDF$geneID %in% present),"Median log2 Fold Change"]),
         enrich.rho = (summaryDF[(summaryDF$geneID %in% present),"Rho_enrich"]),
         deplete.rho = (summaryDF[(summaryDF$geneID %in% present),"Rho_deplete"])
    )   
    #Extract the appropriate stat. 
    values <- switch(stat, 
        enrich.p = sort(summaryDF[,"Target-level Enrichment P"]), 
        deplete.p = sort(summaryDF[,"Target-level Depletion P"]), 
        enrich.fc = sort(-summaryDF[,"Median log2 Fold Change"]), 
        deplete.fc = sort(summaryDF[,"Median log2 Fold Change"]),
        enrich.rho = sort(summaryDF[,"Rho_enrich"]), 
        deplete.rho = sort(summaryDF[,"Rho_deplete"])
    )

    out <- list()
    out$precision <- c(1, unlist(lapply(unique(values), function(x){sum(targvals <= x, na.rm = TRUE)/sum(values <= x, na.rm= TRUE)})), 0)
    out$recall <- c(0, unlist(lapply(unique(values), function(x){sum(targvals <= x, na.rm = TRUE)/length(targvals)})), 1)
    
    enrich <- switch(stat, 
                     enrich.p = ct.targetSetEnrichment(summaryDF, target.list, enrich = TRUE),
                     deplete.p =  ct.targetSetEnrichment(summaryDF, target.list, enrich = FALSE),
                     enrich.fc =  ct.targetSetEnrichment(summaryDF, target.list, enrich = TRUE),
                     deplete.fc =  ct.targetSetEnrichment(summaryDF, target.list, enrich = FALSE),
                     enrich.rho = ct.targetSetEnrichment(summaryDF, target.list, enrich = TRUE),
                     deplete.rho = ct.targetSetEnrichment(summaryDF, target.list, enrich = FALSE)
    )
    out <- c(out, enrich)
    
    #Plot it?
    if(plot.it){
      plot(out$recall, out$precision, xlim = c(0, 1), ylim = c(0,1), 
           type = "l", ylab = "Precision", xlab = "Recall", 
           main = paste("Precision and Recall of", deparse(substitute(target.list))), col = "blue", lwd = 3)
      }
    return(out)
    }


  
  
  
  
  
  
  
  



