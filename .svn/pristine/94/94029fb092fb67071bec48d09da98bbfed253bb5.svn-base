##' @title View a Barchart Summarizing Alignment Statistics for a Crispr Screen
##' @description This function displays the alignemnt statistics for a pooled Crispr screen, reported directly from an alignment statistic matrix. 
##' @param aln A numeric matrix of alignment statistics for a Crispr experiment, typically generated from ExpressionPlot with \code{ep.alignment.class.counts()}. 
##' Alternatively, a 4xN matrix of read counts, with columns indicating samples and rows indicating the number of "targets", "nomatch", "rejections", 
##' and "double_match" reads. Details about these classes may be found in the best practices vignette or as part of the report generated with 
##' \code{ct.makeReport()}. 
##' @param sampleKey An optional ordered factor linking the samples to experimental variables. The \code{names} attribute should exactly match those present in \code{aln}. 
##' @return A grouped barplot displaying the alignment statistics for each sample included in the alignment matrix, which usually corresponds to all of 
##' the samples in the experiment.  
##' @author Russell Bainer
##' @examples 
##' data('aln')
##' ct.alignmentChart(aln)
##' @export

ct.alignmentChart <- function(aln, sampleKey = NULL){

  #input checks
  if(!is.matrix(aln) | !setequal(row.names(aln), c("targets", "nomatch", "rejections", "double_match"))){
    stop("I don't think that the provided alignment matrix is actually an alignment matrix.")
  }  
  
  plotnames <- NULL
  
  if(!is.null(sampleKey)){    
    if(ct.inputCheck(sampleKey, aln)){
      aln <- aln[,names(sampleKey)[order(sampleKey)]]     
      plotnames <- rep('', ncol(aln))
    } 
  }
  
  current.graphic.params <- par(no.readonly = TRUE)
  on.exit(suppressWarnings(par(current.graphic.params)))
  
  par(mar = c(7, 5, 4, 2))
  barplot(aln, main="Read Alignment Statistics by Sample",
          xlab="", ylab ="Reads", col=c("darkblue","darkred", "grey", "purple"), 
          names.arg = plotnames, legend= row.names(aln),
          beside=TRUE, las = 3, cex.names = 0.8, ylim = c(0, 1.3*max(aln)))

  if(!is.null(sampleKey)){    
    axiscolors <- colorRampPalette(c("brown", "green", "orange", "black"))(length(levels(sampleKey)))
    axisloc <- seq(from = 3, by = 5, length.out = length(sampleKey))
    names(axisloc) <- colnames(aln) 
    
    invisible(lapply(levels(sampleKey), function(x){
    axis(1, axisloc[names(sampleKey)[sampleKey == x]], 
         labels = names(sampleKey)[sampleKey == x], tick = FALSE, cex = 0.8, 
         col.axis = axiscolors[match(x, levels(sampleKey))], las = 2)}))
    legend("topleft", levels(sampleKey), fill = axiscolors)
    
    } 
  }
  




















