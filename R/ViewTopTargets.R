##' @title Display the log2 fold change estimates and associated standard deviations of the guides targeting the top 
##' candidates in a crispr screen
##' @description This is a function for displaying candidates from a crispr screen, using the information summarized 
##' in the corresponding \code{fit} and the output from \code{ct.generateResults()}. The fold change and standard deviation 
##' estimates for each gRNA associated with each target (extracted from the \code{coefficients} and \code{stdev.unscaled} slot 
##' of \code{fit}) are plotted on the y axis. Targets are selected on the basis of their gene-level enrichment or depletion 
##' P-values; in the case of ties, they are ranked on the basis of their corresponding Rho statistics. 
##' @param fit An object of class \code{MArrayLM} containing, at minimum, a \code{coefficents} slot with coefficients from the comparison, 
##' and a \code{stdev.unscaled} slot with the corresponding standard deviation of the coefficent estimates. The \code{row.names} attribute 
##' should ideally match that which is found in \code{annotation}.
##' @param summaryDF A data.frame summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param annotation An annotation file for the experiment, usually extracted with \code{ep.load.annot()} in ExpressionPlot. gRNAs are annotated by 
##' row, and must minimally contain a column \code{geneSymbol}. 
##' @param targets Either the number of top targets to display, or a list of \code{geneSymbol}s contained in the \code{geneSymbol} 
##' slot of the \code{annotation} object. 
##' @param enrich Logical indicating whether to display guides that are enriched (default) or depleted within the screen. If a vector of 
##' \code{geneSymbol}s is specified, this controls the left-t0-right ordering of the corresponding gRNAs. 
##' @param contrast.term If a fit object with multiple coefficients is passed in, a string indiating the coefficient of interest.   
##' @return An image on the default device indicating each gRNA's log2 fold change and the unscaled standard deviation of the effect estimate, 
##' derived from the \code{MArrayLM} object.
##' @author Russell Bainer
##' @examples 
##' data('fit')
##' data('resultsDF')
##' data('ann')
##' 
##' ct.topTargets(fit, resultsDF, ann) 
##' @export

ct.topTargets <- function(fit, summaryDF, annotation, targets = 10, enrich = TRUE, contrast.term = NULL){
  current.graphic.params <- par(no.readonly = TRUE)
  on.exit(suppressWarnings(par(current.graphic.params)))

  if(ncol(fit$coefficients) > 1){
    if(is.null(contrast.term)){
      stop("The fit object contains multiple coefficients. Please specify a contrast.term.")
    }
    fit <- ct.preprocessFit(fit, contrast.term)
  }
  
  
    #Test input: 
    #testing
    if(class(fit) != "MArrayLM"){stop(paste(deparse(substitute(eset)), "is not an MArrayLM."))}
    
    if(!setequal(row.names(annotation), row.names(fit))){
      warning("row.names of the fit object and the annotation file are not identical. Using the intersection only.")
        grnas <- intersect(row.names(fit), row.names(annotation))
        fit <- fit[grnas,]
        annotation <- annotation[grnas,]
    }

  if(!(enrich %in% c(TRUE, FALSE))){
    stop('enrich must be either TRUE or FALSE.')
  }
  
    if(!ct.resultCheck(summaryDF)){
      stop("Execution halted.")
    }
  

    #Identify the top targets from the summary DF; order and group gRNAs within a target, then rank targets
    summaryDF$geneSymbol <- as.character(summaryDF$geneSymbol)
    summaryDF <- summaryDF[with(summaryDF, 
                               order(summaryDF[,"Target-level Enrichment P"], 
                                      summaryDF[,"RhoRank_enrich"], 
                                      summaryDF[,"geneSymbol"], 
                                      -summaryDF[,"gRNA Log2 Fold Change"])),]   
    plottitle <- "Enriched Targets"
    
    if(enrich == FALSE){
      summaryDF <- summaryDF[with(summaryDF, 
                                  order(summaryDF[,"Target-level Depletion P"], 
                                        summaryDF[,"RhoRank_deplete"], 
                                        summaryDF[,"geneSymbol"], 
                                        summaryDF[,"gRNA Log2 Fold Change"])),]   
      plottitle <- "Depleted Targets"
    }

    if(is.character(targets)){
      toptargets <- intersect(targets, annotation$geneSymbol)
      ntargets <- length(toptargets)
      plottitle <- ''
    } else {
      if((length(targets) != 1) | !is.numeric(targets)){stop('"targets" must be specified as a single number or a vector of elements in the geneSymbol column of the annotation object.')}
      ntargets <- targets
      plottitle <- paste('Top', ntargets, plottitle)
      toptargets <- unique(summaryDF$geneSymbol)[1:ntargets]
    } 
    if(ntargets <= 0){stop("No valid targets were specified.")}
    
    targetrows <- row.names(summaryDF)[(summaryDF$geneSymbol %in% toptargets)]
    nguides <- unlist(lapply(toptargets, function(x){sum(summaryDF$geneSymbol %in% x, na.rm = TRUE)}))

    lfc <- fit$coefficients[targetrows,1]
    sdu <- fit$stdev.unscaled[targetrows,1] + fit$coefficients[targetrows,1]
    sdl <- fit$coefficients[targetrows,1] - fit$stdev.unscaled[targetrows,1]
    
    ylimit <- c(min(0, min(sdl)), max(0, max(sdu)))
    
    #Compose a vector of the x locations
    xloc <- as.vector(unlist(mapply(function(start, numguides){if(numguides > 1){seq(start, start+1, length.out = numguides)} else{start + 0.5}}, 
                   start = seq(1, by = 2, length.out = ntargets), 
                   numguides = nguides, 
                   SIMPLIFY = TRUE)))
  
    sddf <- cbind(xloc, sdu, sdl)  
    
    #make the plot
    plot(xloc, lfc, xaxt='n', main = plottitle, ylim = ylimit, xlab = "", ylab = "Log2 gRNA Abundance Change", pch = 18, col = "darkred")
    suppress <- apply(sddf, 1, function(x){lines(c(x[1], x[1]), c(x[2], x[3]), col = rgb(0,0,1,0.3), lwd = 4)})
    points(xloc, lfc, pch = 18, col = "darkred")
    suppress <- capture.output(lapply(seq(2.5, by = 2, length.out = (ntargets - 1)), function(x){abline(v = x, lty = "dotted", col= "gray")}))
    abline(h = 0)
    #Add the labels
    axis(1, at = seq(1.5, by = 2, length.out = ntargets), labels = toptargets, las = 3)
    if (!is.null(contrast.term)) {
        mtext(contrast.term, side=3, line=0)
    }
  }

  
  
  
  
  
  
  
  
  
  



