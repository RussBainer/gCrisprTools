##' @title View a stacked representation of the most variable targets or individual guides within an experiment, 
##' as a percentage of the total aligned reads
##' @description This function identifies the gRNAs or targets that change the most from sample to sample within an experiment as a percentage of 
##' the entire library. It then plots the abundance of the top \code{nguides} as a stacked barplot for all samples in the experiment. The purpose of this 
##' algorithm is to detect potential distortions in the library composition that might not be properly controlled by sample normalization, and so 
##' the most variable entites are defined by calculating the percent of aligned reads that they contribute to each sample, and then ranking each entity
##' by the range of these percentages across all samples. Consequently, gRNAs or Targets that are highly abundant in at least one condition will be 
##' are more likely to be identified. 
##' @param eset An ExpressionSet  object containing, at minimum, a matrix of gRNA abundances extractable with the exprs() function, and a metadata 
##' object containing a column named \code{SAMPLE_LABEL} containing unique identifers for each sample. The \code{colnames} should be syntactically 
##' @param sampleKey An optional sample key, supplied as an ordered factor linking the samples to experimental 
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control set is assumed to be 
##' the first \code{level}.
##' @param nguides The number of guides (or targets) to display.    
##' @param plotType A string indicating whether the individual guides should be displayed ("\code{gRNA}"), or if they should be aggregated into target-level
##' estimates ("\code{Target}") according to the \code{geneSymbol} column in the \code{annotation} object.  
##' @param annotation An optional data.frame containing an annotation object to be used to aggregate the guides into targets. gRNAs are annotated by row, 
##' and must minimally contain a column \code{geneSymbol} indicating the target elements.
##' @param ylimit An optional numeric vector of length 2 specifying the y limits for the plot, useful in comparin across studies. 
##' @param subset An optional character vector containing the sample labels to be used in the analysis; all elements must be contained in the \code{colnames} of the specified \code{eset}. 
##' @return A stacked barplot displaying the appropriate entities on the default device. 
##' @author Russell Bainer
##' @import ggplot2
##' @examples 
##' data('es')
##' data('ann')
##' ct.stackGuides(es, nguides = 20, plotType = "Target", annotation = ann, ylimit = NULL, subset = NULL)
##' @export

ct.stackGuides <- function(eset, sampleKey = NULL, nguides = 20, plotType = "gRNA", annotation = NULL, ylimit = NULL, subset = NULL){
  current.graphic.params <- par(no.readonly = TRUE)
  on.exit(suppressWarnings(par(current.graphic.params)))

  if (!requireNamespace("ggplot2")) {
    stop("The ggplot2 package is required")
  }
  if(!is.numeric(nguides)){stop('Please specify a numeric number of guides to display.')}  
  if(!methods::is(eset, "ExpressionSet")){stop('eset must be an expressionset object.')}  
  if(!(plotType %in% c("gRNA", "Target"))){stop('Please specify "gRNA" or "Target" to be displayed.')}

  #Check eset colnames
  d <- exprs(eset)  
  
  if(is.null(sampleKey)){
    sampleKey <- ordered(rep('', ncol(d)))
    names(sampleKey) <- colnames(d)
  }
  
  if(any(make.names(colnames(d)) != colnames(d))){
    warning("Some of the sample names are not syntactically valid. Coercing.")  
    sampleKey <- sampleKey[colnames(d)]
    names(sampleKey) <- make.names(colnames(d))
    colnames(d) <- make.names(colnames(d))
  }  
  
  #If a fit is included, check to see if it's valid and then subset the expression and pheno data appropriately. 
  if(!is.null(subset)){
    if(!is.character(subset)){
      stop("subset should be a character vector containing the labels of the samples that you wish to analyze.")
    }
    subset <- make.names(subset)
    
    if(length(setdiff(subset, colnames(d))) != 0){
      stop("Not all of the samples in subset are present in the specified eset.")
    }

    d <- d[,subset]
    sampleKey <- sampleKey[subset]
  }  
  
  #idiosyncracies of ggplot forces rearrangement of factor labels for proper plotting. 
  sampleKey <- sampleKey[order(sampleKey)]

  invisible(ct.inputCheck(sampleKey, d))

  #Everything ok, moving on. 
  plottitle <- paste0("Top ", nguides, " Most Variable ", plotType, "s Across Experimental Condition")
  
  if(plotType == "Target"){
    if(is.null(annotation)){
      stop('An annotation object containing a "geneSymbol" column must be supplied to 
           display target-level representation.')
      }    

    annotation <- ct.prepareAnnotation(annotation, throw.error = FALSE)
    
    if(sum(is.na(annotation$geneSymbol)) > 0){
      message('Converting missing values in the annotation file to "NoTarget".')
      annotation$geneSymbol[is.na(annotation$geneSymbol)] <- "NoTarget"
      }

    #Convert all gRNA abundance to % representation
    message('Summarizing gRNA counts into targets.')
    d <- t(vapply(levels(annotation$geneSymbol), 
                        function(x){if(sum(annotation$geneSymbol %in% x) > 1){
                                     colSums(d[row.names(annotation)[annotation$geneSymbol %in% x],])
                                     } else {
                                       d[row.names(annotation)[annotation$geneSymbol %in% x],]
                                       }
                                    }, 
                        numeric(ncol(d))
                        )
                 )
    
    }
  
  d <- apply(d, 2, function(x){x/sum(x, na.rm = TRUE)})
  d <- d[order(apply(d, 1, function(x){range(x, na.rm = TRUE)[2] - range(x, na.rm = TRUE)[1]}), 
               decreasing = TRUE),
         names(sampleKey)[order(sampleKey)]]
  d <- d[seq_len(nguides),]
    
  plotframe <- data.frame("gRNA" = rep(row.names(d), ncol(d)), 
                          "Condition" = rep(paste0(colnames(d), '_', sampleKey[colnames(d)]), each = nrow(d)), 
                          "ReadProportion" = as.numeric(d))
  
  colorScale <- colorRampPalette(c("white", "red",  "blue", "black"))(nguides)
  
  #as requested by Sarah, highlight the NTCs if they are present.
  if(is.element("NoTarget", levels(plotframe$gRNA))){
    colorScale[is.element(levels(plotframe$gRNA), "NoTarget")] <- "forestgreen"    
    }
  
  #also scale the legend as appropriate for the number of guides: 
  legend.scale.factor <- 15/nguides
  
  if(!is.null(ylimit)){
    if(!is.numeric(ylimit) | (length(ylimit) != 2)){
      stop("The ylimit variable must be NULL, or a numeric vector of length 2.")
    }
    ggplot(plotframe, aes_string(x = 'Condition', y = 'ReadProportion', fill = 'gRNA')) + 
      geom_bar(stat = "identity", position = 'stack') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      scale_fill_manual(values=colorScale) + ggtitle(plottitle) + 
      coord_cartesian(ylim = c(ylimit[1], ylimit[2])) + 
      ylab("Proportion of Total Reads") + 
      theme(legend.key.size = unit(legend.scale.factor, "cm"), legend.title=element_blank()) 
    
  } else{
    ggplot(plotframe, aes_string(x = 'Condition', y = 'ReadProportion', fill = 'gRNA')) + 
      geom_bar(stat = "identity", position = 'stack') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      scale_fill_manual(values=colorScale) + 
      ggtitle(plottitle) + 
      ylab("Proportion of Total Reads") + 
      theme(legend.key.size = unit(legend.scale.factor, "cm"), legend.title=element_blank())
  }
}





  
  
  
  
  
  
  
  



