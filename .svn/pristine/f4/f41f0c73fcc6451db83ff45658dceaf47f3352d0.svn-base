##' @title Assign Colors Based on the Position of a Value in a Distribution
##' @description This is a function to generate colors for plot elements on the basis of the position of a value within a distribution. 
##' Called internally by ct.viewGuides. 
##' @param exprs The value whose color is to be returned. 
##' @param rankedexprs A vector of values the length of cols that corresponds to the values in the distribution.
##' @param colors The vector of colors to be used as a reference. 
##' @return A value contained in cols. 
##' @keywords internal
##' @author Russell Bainer
ct.exprsColor <- function(exprs, rankedexprs, colors){
    colors[which(abs(rankedexprs-exprs) == min(abs(rankedexprs-exprs)))]
  }
  
##' @title Draw a horizontal line of a specified color.
##' @description This is a function called internally by ct.viewGuides to generate the color legend. End users should not use it. 
##' @param x,y Minimal coordinates to specify the line, which is drawn from the Y axis. 
##' @param color Guess!
##' @param width Line width. 
##' @return A line on an open device. 
##' @keywords internal
##' @author Russell Bainer
ct.drawFlat <- function(x, y, color, width = 1){
    lines(c(0, x), c(y, y), col = color, lwd = width);
  }

##' @title Draw a density color legend. 
##' @description This is a function called internally by ct.viewGuides to generate the color legend. End users should not use it. 
##' @param dens A density object. 
##' @param colorscale A vector of colors to draw behind the density. 
##' @return A color legend on the current graphics device. 
##' @keywords internal
##' @author Russell Bainer
ct.drawColorLegend <- function(dens, colorscale){
    yrange <- (rep(max(abs(range(dens$x))), 2) * c(-1, 1))
  
    plot(NA, 
         xlim = c(0,max(dens$y)), 
         ylim = range(dens$x), 
         xaxt = "n", xlab = "", 
         ylab = "Log2 Normalized Reads", 
         bty="n")
    invisible(mapply(ct.drawFlat, 
                     x = rep(max(dens$y), length(colorscale)), 
                     y = seq(yrange[2], yrange[1], length.out = length(colorscale)), 
                     width = 2, 
                     color = colorscale))
    polygon(x=dens$y, y= dens$x, col = "black")
  }
  
##' @title Generate a Plot of individual gRNA Pair Data in a Crispr Screen
##' @description This function generates a visualization of the effect estimates from a MArrayLM model result for all of the 
##' individual guides targeting a particular element, specified somewhere in the library annotation file. The estimated effect 
##' size and variance is plotted relative to zero for the specified contrast, with the color of the dot indicating the relative scale of the
##' of the guide intercept within the model framework, with warmer colors indicating lowly expressed guides. For comparison, the density of gRNA 
##' fold change estimates is privided in a pane on the right, with white lines indicating the exact levels of the individual guides. 
##' @param gene the name of the target element of interest, contained within the "type" column of the annotation file.
##' @param fit An object of class MArrayLM containing, at minimum, an "Amean" slot containing the guide level abundances, 
##' a "coefficients" slot containing the effect estimates for each guide, and an "stdev.unscaled" slot giving the coefficient standard Deviations. 
##' @param ann A data.frame object containing the gRNA annotations, usually returned from the ExpressionPlot function 'ep.load.annot()'. 
##' At mimimum, it should have a column with the name specified in the type argument, containing the element targeted by each guide. 
##' @param type A character string indicating the column in ann containing the target of interest. 
##' @param contrast.term If a fit object with multiple coefficients is passed in, a string indiating the coefficient of interest.   
##' @return An image summarizing gRNA behavior within the specifed gene on the default device. 
##' @author Russell Bainer
##' @examples
##' data('fit')
##' data('ann')
##' ct.viewGuides('Target1633', fit, ann)
##' @export

ct.viewGuides <- function(gene, fit, ann, type = "geneSymbol", contrast.term = NULL){
  current.graphic.params <- par(no.readonly = TRUE)
  on.exit(suppressWarnings(par(current.graphic.params)))
  
  par(xpd = FALSE)

  if(ncol(fit$coefficients) > 1){
    if(is.null(contrast.term)){
      stop("The fit object contains multiple coefficients. Please specify a contrast.term.")
    }
    fit <- ct.preprocessFit(fit, contrast.term)
  }
  
  
  
  #testing
  if(class(fit) != "MArrayLM"){stop(paste(deparse(substitute(eset)), "is not an MArrayLM."))}
 
  #Find the gRNAs targeting the gene from the annotation, and order them
  options(warn=-1)
  ann <- ct.prepareAnnotation(ann, fit)
  options(warn=0)
  
  if(!(sum(ann[,type] %in% gene))){stop(paste(gene, "is not present in the annotation file."))}
  
  grna.inx <- row.names(ann)[(ann[,type] %in% gene)]
  grna.inx <- sort(grna.inx)

  #Set up the color palette and the legend first: 
  day3Density <- density(fit$coefficients)
  ylimit <- max(fit$coefficients)  
  ylims <- range(c(min(0,(fit$coefficients[grna.inx,1]-fit$stdev.unscaled[grna.inx,1])), 
                   max(0,(fit$coefficients[grna.inx,1]+fit$stdev.unscaled[grna.inx,1]))))

  #Set up the guide colors based on their overall expression
  cols <- colorRampPalette(c("red", "orange", "white", "blue", "purple"))(1000)
  rankexprscolors <- seq(from = max(fit$Amean), to = min(fit$Amean), length.out = 1000)

  #Set up a couple of panes: 
  layout(matrix(c(1,2),ncol=2,byrow=TRUE),widths=c(5,0.5),heights=3,respect=TRUE)
  par(mar = c(5,4,4,2))
  
  plot(1:length(grna.inx), 
       fit$coefficients[grna.inx,1],
       ylab = "Log2(Fold Change)", 
       main = gene, xlab = "", ylim = ylims,
       col = vapply(fit$coefficients[grna.inx,1], ct.exprsColor, character(1), rankedexprs = rankexprscolors, colors = cols), 
       xaxt = "n", pch = 18)
  segments(seq_len(length(grna.inx)), 
           fit$coefficients[grna.inx,1], 
           seq_len(length(grna.inx)), 
           (fit$coefficients[grna.inx,1]-fit$stdev.unscaled[grna.inx,1]))
  segments(seq_len(length(grna.inx)), 
           fit$coefficients[grna.inx,1],
           seq_len(length(grna.inx)), 
           (fit$coefficients[grna.inx,1]+fit$stdev.unscaled[grna.inx,1]))
  
  
  axis(1, at=1:length(grna.inx), labels=grna.inx, las = 2)
  abline(h = 0, lty = 4, col = "red")
  abline(v = (seq_len(length(grna.inx)) + 0.5), col = "grey")
  #title(main = )
  if (!is.null(contrast.term)) {
      mtext(contrast.term, side=3, line=0)
  }
  #Plop in the legend:
  par(mar = c(5,1,4,1))
  
  ct.drawColorLegend(dens = day3Density, colorscale = cols)
  mapply(ct.drawFlat, x = rep(max(day3Density$y), length(grna.inx)),
         y = fit$coefficients[grna.inx,1], width = 0.5, color = "black")
  lines(c(0, max(day3Density$y)), c(0, 0), col = 'red', lty = 'dashed')
  return(TRUE)
}
















