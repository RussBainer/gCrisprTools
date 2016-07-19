##' @title Visualization of Raw gRNA Count Densities
##' @description This function plots the per-sample densities of raw gRNA read counts on the log10 scale. The curve colors are assigned based on a user-
##' specified sampleKey. This function is primarily useful to determine whether libraries are undersequenced (low mean raw gRNA counts), 
##' contaminated (many low-abundance gRNAs present), or if PCR artifacts may be present (subset of extremely abundant guides, multiple gRNA distribution modes). In most 
##' well-executed experiments the majority of gRNAs will form a tight distribution around some reasonably high average read count (hundreds of reads), 
##' at least among the control samples. Excessively low raw count values can compromise normalization steps and subsequent estimation of gRNA levels, especially 
##' in screens in which most gRNAs have minimal effects on cell viability. 
##' @param eset An ExpressionSet object containing, at minimum, count data accessible by exprs() and some phenoData. 
##' @param sampleKey A sample key, supplied as a (possibly ordered) factor linking the samples to experimental 
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control set 
##' is assumed to be the first \code{level}.
##' @return A density plot as specified on the default device. 
##' @author Russell Bainer
##' @examples 
##' data('es')
##' 
##' #Build the sample key
##' library(Biobase)
##' sk <- relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference")
##' names(sk) <- row.names(pData(es))
##' 
##' ct.rawCountDensities(es, sk)
##' @export
ct.rawCountDensities <- function(eset, sampleKey){

  if(class(eset) != "ExpressionSet"){stop(paste(deparse(substitute(eset)), "is not an ExpressionSet."))}
  
  if (is.null(sampleKey)) {
    sampleKey <- as.factor(colnames(eset))
    names(sampleKey) <- sampleKey
  } else {
    ct.inputCheck(sampleKey, eset)
    sampleKey <- sampleKey[order(sampleKey)]
  }
    
  e.dat <- log10(exprs(eset) + 1)
  densities <- apply(e.dat, 2, density)
  
  y <- c(0, max(unlist(lapply(densities, function(dens){max(dens$y)}))))
  x <- c(0,  max(ceiling(unlist(lapply(densities, function(dens){max(dens$x)})))));  
  plot(x[1], y[1], 
       xlim = x, ylim = y, 
       xlab = "gRNA Read Counts (Raw)", 
       ylab = "Density", 
       pch = NA, 
       main = "Raw gRNA Counts (Log Scale)", 
       xaxt = "n")
  axis(1, at= 0:x[2], labels = 10^(0:x[2]))
 
  #Set up colors for the factor
  colors <- colorRampPalette(c("blue", "red"), alpha = TRUE)(length(levels(sampleKey)))
  colors <- gsub("FF$", "99", colors, perl = TRUE)
  
  invisible(lapply(seq_len(length(densities)), 
                   function(x){lines(densities[[x]], 
                                     col = colors[as.numeric(sampleKey[colnames(e.dat)[x]])])}
                   )
            )
  legend("topleft", legend = levels(sampleKey), fill = colors)

}


