##' @title View nontargeting guides within an experiment
##' @description This function tries to identify, and then plot the abundance of, the full set of non-targeting controls from an ExpressionSet
##' object. Ideally, the user will supply a geneSymbol present in the appropriate annotation file that uniquely identifies the nontargeting gRNAs.
##' Absent this, the the function will search for common identifier used by nontargeting controls (geneID "no_gid", or geneSymbol NA).
##' @param eset An ExpressionSet  object containing, at minimum, a matrix of gRNA abundances extractable with the \code{exprs} function.
##' @param annotation An annotation data.frame for the experiment. gRNAs are annotated by
##' row, and must minimally contain columns \code{geneSymbol} and \code{geneID}.
##' @param sampleKey A sample key, supplied as an ordered factor linking the samples to experimental
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control condition is assumed to be
##' the first \code{level}.
##' @param geneSymb The \code{geneSymbol} identifier in \code{annotation} that corresponds to nontargeting gRNAs. If absent, \code{ct.ViewControls} will
##' attempt to infer nontargeting guides by searching for \code{"no_gid"} or \code{NA} in the appropriate columns.
##' @param normalize Logical indicating whether to attempt to normalize the data in the \code{eset} by DESeq size factors present in the metadata. If \code{TRUE},
##' then the metadata must contain a column containing these factors, named \code{sizeFactor.crispr-gRNA}.
##' @param lib.size An optional vector of voom-appropriate library size adjustment factors, usually calculated with \code{\link[edgeR]{calcNormFactors}} 
##' and transformed to reflect the appropriate library size. These adjustment factors are interpreted as the total library sizes for each sample, 
##' and if absent will be extrapolated from the columnwise count sums of the \code{exprs} slot of the \code{eset}.
##' @return An image of nontargeting control gRNA abundances on the default device.
##' @author Russell Bainer
##' @examples 
##' data('es')
##' data('ann')
##' 
##' #Build the sample key
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' 
##' ct.viewControls(es, ann, sk, geneSymb = NULL, normalize = FALSE)
##' ct.viewControls(es, ann, sk, geneSymb = NULL, normalize = TRUE)
##' @export

ct.viewControls <- function(eset, annotation, sampleKey, geneSymb = NULL, normalize = TRUE, lib.size = NULL){
  #current.graphic.params <- par(no.readonly = TRUE)
  #on.exit(suppressWarnings(par(current.graphic.params)))

  annotation <- ct.prepareAnnotation(annotation, object = eset, throw.error = FALSE)
  invisible(ct.inputCheck(sampleKey, eset))

  if(is.null(lib.size)){
    lib.size <- colSums(exprs(eset))
  } else if(!is.numeric(lib.size) | (length(lib.size) != ncol(eset))){
    stop('If specified, lib.size must be a numeric vector of length equal to the number of samples.')
  }
  
  #Find the row.names that correspond to the guides.
  if(!is.null(geneSymb)){
    if(!(geneSymb %in% annotation$geneSymbol)){geneSymb <- NULL}
    }

  if(is.null(geneSymb)){
    message("No control geneSymbol supplied, so I'll use the default of 'NoTarget'.")
    ntc <- row.names(annotation)[annotation$geneSymbol == "NoTarget"]
    } else{
      ntc <- row.names(annotation)[annotation$geneSymbol == geneSymb]
      }

  if(length(ntc) < 1){
    stop('No controls detected.')
  }
  
  colorSpace <- colorRampPalette(c("lightblue", "darkred"))(length(ntc))
  plottitle <- paste0(geneSymb, " Guide Abundance (Raw Reads)")

  if(normalize){
    message('Normalizing gRNAs by median scaling.')
    if(is.null(annotation)){
        eset <- ct.normalizeGuides(eset, 'scale', sampleKey = sampleKey, lib.size = lib.size, plot.it = FALSE)
    } else {
      eset <- ct.normalizeGuides(eset, 'scale', annotation = annotation, sampleKey = sampleKey, lib.size = lib.size, plot.it = FALSE)
      }
    plottitle <- paste0(geneSymb, " Guide Abundance (Median Scaled)")
  }

  counts <- exprs(eset)[ntc,names(sampleKey)[order(sampleKey)]]
  counts <- (log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  counts <- counts[,order(colMeans(counts))]

    #Set up and draw the plot
    ylimit <- range(counts) + c(-1,1)
    xlimit <- nrow(counts)

par(mar = c(10,4,4,6), xpd = TRUE)
plot(1:xlimit, counts[,1],
         main = plottitle,
         ylab = "Log2 Normalized Counts",
         xlab = "",
         xaxt = 'n',
         xlim = c(1,xlimit),
         ylim = ylimit,
         type = "l", lwd = 2,
         col = colorSpace[1])
for(q in 2:ncol(counts)){lines(1:xlimit, counts[,q], lwd = 2, col = colorSpace[q])}
axis(side = 1, labels = sampleKey[order(sampleKey)], at = 1:xlimit, las = 3)
legend(xlimit + 1, ylimit[2], legend = colnames(counts), fill = colorSpace, cex = 0.5)
}




















