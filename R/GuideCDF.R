##' @title Generate a cumulative tally of reads by guide rank
##' @description This function returns a numeric vector of the same length as the input, where each element \code{n} contains the proportion
##' the sum of the full vector that is captured by its first \code{1-n} elements (arranged in descending order).  
##' @param vector An input numeric vector to be aggregated. 
##' @return A CDF plot displaying the appropriate CDF curves on the default device. 
##' @author Russell Bainer
##' @keywords internal
##' @examples v <- sort(sample(1:100, 30, replace = TRUE, 100:1))
##' ct.ecdf(v)
##' @export

ct.ecdf <- function(vector){
  sorted <- as.numeric(sort(vector, decreasing = TRUE))
  out <- unlist(lapply(seq_along(sorted), function(x){sum(sorted[1:x])}))/sum(sorted)
  return(out)
  }

##' @title View CDFs of the ranked gRNAs or Targets present in a crispr screen
##' @description This function generates a plot relating the cumulative proportion of reads in each sample of a crispr screen to the abundance rank of the 
##' underlying guides (or Targets). The purpose of this algorithm is to detect potential distortions in the library composition 
##' that might not be properly controlled by sample normalization (see also: \code{ct.stackedGuides()}).
##' @param eset An ExpressionSet  object containing, at minimum, a matrix of gRNA abundances extractable with the exprs() function.
##' @param sampleKey An optional sample key, supplied as an ordered factor linking the samples to experimental 
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control set is assumed to be 
##' the first \code{level}.
##' @param plotType A string indicating whether the individual guides should be displayed ("\code{gRNA}"), or if they should be aggregated into target-level
##' estimates ("\code{Target}") according to the \code{geneSymbol} column in the \code{annotation} object.  
##' @param annotation An optional data.frame containing an annotation object to be used to aggregate the guides into targets. gRNAs are annotated by row, 
##' and must minimally contain a column \code{geneSymbol} indicating the target elements.
##' @return A CDF plot displaying the appropriate CDF curves on the default device. 
##' @author Russell Bainer
##' @examples data('es')
##' ct.guideCDF(es)
##' @export

ct.guideCDF <- function(eset, sampleKey = NULL, plotType = "gRNA", annotation = NULL){

  current.graphic.params <- par(no.readonly = TRUE)
  on.exit(suppressWarnings(par(current.graphic.params)))

  if(class(eset) != "ExpressionSet"){stop('eset must be an expressionset object.')}  
  if(!(plotType %in% c("gRNA", "Target"))){stop('Please specify "gRNA" or "Target" to be displayed.')}
  
  #Extract and preprocess data  
  d <- exprs(eset)  
  legendnames <- colnames(d)
  
  if(!is.null(sampleKey)){
    if(ct.inputCheck(sampleKey, eset)){
      d <- d[,names(sampleKey)[order(sampleKey)]]
      legendnames <- paste(sampleKey[colnames(d)], colnames(d), sep = "_")
    } else {
      warning('The supplied sampleKey is incompatible with the eset. Ignoring.')  
    }
    
  }
  if(plotType == "Target"){
    if(is.null(annotation) | !("geneSymbol" %in% names(annotation))){
      stop('An annotation object containing a "geneSymbol" column must be supplied 
           to display target-level representation.')
      }    
    
    annotation <- ct.prepareAnnotation(annotation, eset, throw.error = FALSE)
    
    message('Summarizing gRNA counts into targets.')
    
    genects <- lapply(unique(annotation$geneSymbol), function(x){
      if(sum(annotation$geneSymbol %in% x) > 1){
        colSums(d[row.names(annotation)[annotation$geneSymbol %in% x],])
      } else {
        d[row.names(annotation)[annotation$geneSymbol %in% x],]
      }})
    d <- data.frame(t(simplify2array(genects)))
    row.names(d) <- unique(annotation$geneSymbol)    
  }
  
  d.rank <- apply(d, 2, ct.ecdf)
  
  #plot it
  colorSpace <- colorRampPalette(c("blue", "lightblue", "darkred"))(ncol(d.rank)) 
  
  plot(seq_len(nrow(d)), d.rank[,1], 
       pch = 20, col = colorSpace[1], 
       xlim = c(0, nrow(d.rank)), 
       xlab = paste(plotType, "Abundance Rank"), 
       ylab = "Cumulative Proportion of Reads", 
       main = paste("Cumulative Proportion of Reads by", plotType, "Rank"))
  
  for(z in 2:ncol(d.rank)){
    lines(seq_len(nrow(d)), d.rank[,z], 
          pch = 20, 
          col = colorSpace[z], 
          new = FALSE)
    }
  legend("bottomright", legend = legendnames, fill = colorSpace, ncol = 2, cex = 0.5)  
  }
  


