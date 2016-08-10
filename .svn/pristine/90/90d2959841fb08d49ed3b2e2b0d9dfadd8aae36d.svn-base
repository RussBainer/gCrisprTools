##' @title Visualization of Ranked gRNA Abundances by Replicate
##' @description This function median scales and log2 transforms the raw gRNA count data contained in an ExpressionSet, 
##' and then plots the ordered expression values within each replicate. The curve colors are assigned based on a user-
##' specified column of the pData contained in the ExpressionSet. Optionally, this function can plot the location of Nontargeting control 
##' guides (or any guides, really) within the distribution. 
##' @param eset An ExpressionSet object containing, at minimum, count data accessible by exprs() and some phenoData. 
##' @param sampleKey A sample key, supplied as a (possibly ordered) factor linking the samples to experimental 
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control set 
##' is assumed to be the first \code{level}.
##' @param annotation An annotation dataframe indicating the nontargeting controls in the geneID column. 
##' @param geneSymb The \code{geneSymbol} identifier(s) in \code{annotation} that corresponds to gRNAs to be plotted on the curves. 
##' If the provided value is not present in the \code{geneSymbol}, nontargeting controls will be plotted instead.
##' @param lib.size An optional vector of voom-appropriate library size adjustment factor, usually calculated with \code{edgeR::calcNormFactors}.
##' @return A waterfall plot as specified, on the default device.
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
##' ct.gRNARankByReplicate(es, sk, ann, 'Ripk3')
##' @export
ct.gRNARankByReplicate <- function(eset, sampleKey, annotation = NULL, geneSymb = NULL, lib.size = NULL){
  #current.graphic.params <- par(no.readonly = TRUE)
  #on.exit(suppressWarnings(par(current.graphic.params)))

  if(class(eset) != "ExpressionSet") {
    stop(paste(deparse(substitute(eset)), "must be an ExpressionSet."))
  }
  
  if (is.null(sampleKey)) {
    sampleKey <- as.factor(colnames(eset))
    names(sampleKey) <- sampleKey
  } else {
    ct.inputCheck(sampleKey, eset)
    sampleKey <- sampleKey[order(sampleKey)]
  }
  
  counts <- exprs(eset)
  
  if (is.null(lib.size)){
    lib.size <- colSums(counts)
  }else if(!is.numeric(lib.size) | length(lib.size) != ncol(counts)){
    stop('If specified, lib.size must be a numeric vector of the same 
         length as the number of samples in the eset.')
  } 
  
  #Convert to CPM
  e.dat <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))

  y <- range(e.dat)
  x <- c(1, nrow(e.dat)) 
  if(is.null(annotation) | is.null(geneSymb)){
    plottitle <- '' 
    } else if(geneSymb %in% annotation$geneSymbol){
      plottitle <- geneSymb
      } else {
        plottitle <- "Nontargeting Controls"
        }
    
  plot(x[1], y[1], 
       xlim = x, ylim = y, 
       xlab = "gRNA Abundance Rank", ylab = "Log2 Counts", 
       pch = NA, 
       main = plottitle)
  
  colors <- colorRampPalette(c("blue", "red"), alpha = TRUE)(length(levels(sampleKey)))
  colors <- gsub("FF$", "99", colors, perl = TRUE)
  
  invisible(lapply(seq_len(ncol(e.dat)), 
                   function(x){lines(seq_len(nrow(e.dat)), 
                                     sort(e.dat[,names(sampleKey)[x]], 
                                          decreasing = TRUE), 
                                     col = colors[as.numeric(sampleKey)[x]])}
                   )
            )

  #Add the NTC locations if requested. 
    if(!is.null(geneSymb)){

      if(is.null(annotation) | !("geneSymbol" %in% names(annotation))){
        stop("An annotation dataframe must be supplied if geneSymb is not NULL.")
        }

      annotation <- ct.prepareAnnotation(annotation, eset, throw.error = FALSE)

      if(any(geneSymb %in% annotation$geneSymbol)){
        ntc <- row.names(annotation)[annotation$geneSymbol %in% geneSymb]      
        } else {
          ntc <- row.names(annotation)[(annotation$geneSymbol %in% "NoTarget")]
          } 
      
      if(length(ntc) == 0){
            stop('No suitable elements are present in the supplied annotation file. 
                 Please specify a geneSymbol if you want to display individual guides.')
          }

    #make a table of the ntc locations and ranks
    ntc.ranks <- (apply(e.dat, 2, rank))
    ntc.ranks <- ntc.ranks[ntc,]
    ntc.ranks <- nrow(e.dat) - ntc.ranks
    ntc.dat <- e.dat[ntc,]
    
    rimcolors <- colorRampPalette(c("black", "blue", "green", "red", "yellow", "white"))(nrow(ntc.dat))
    
    invisible(lapply(seq_len(ncol(ntc.dat)), 
                     function(x){points(ntc.ranks[,names(sampleKey)[x]], 
                                        ntc.dat[,names(sampleKey)[x]], 
                                        bg = rimcolors, 
                                        col = colors[as.numeric(sampleKey[x])], 
                                        pch = 23, lwd = 4)
                                  }
                     )
              )
    }
  legend("topright", legend = levels(sampleKey), fill = colors)
}



