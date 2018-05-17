##' @title Calculate results of a crispr screen from a contrast
##' @description This is a wrapper function that enables direct generation of target-level p-values from a crispr screen.  
##' @param fit An object of class \code{MArrayLM} containing, at minimum, a \code{t} slot with t-statistics from the comparison, 
##' a \code{df.residual} slot with the corresponding residuals fo the model fits, and an \code{Amean} slot with the respective mean abundances. 
##' @param annotation An annotation file for the experiment. gRNAs are annotated by 
##' row, and must minimally contain columns \code{geneSymbol} and \code{geneID}.
##' @param RRAalphaCutoff A cutoff to use when defining gRNAs with significantly altered abundance during the RRAa aggregation step. If \code{scoring} is set 
##' to \code{pvalue} or \code{combined}, this parameter is interpreted as the maximum nominal p-value required to consider a gRNA's abundance meaningfully 
##' altered during the aggregation step. If \code{scoring} is \code{fc}, this parameter is interpreted as the proportion of the list to be considered 
##' meaningfully altered in the experiment (e.g., if \code{RRAalphaCutoff} is set to 0.05, only consider the rankings of the 5% most upregulated 
##' (or downregulated) gRNAs for the purposes of RRAa calculations).
##' @param permutations The number of permutations to use during the RRAa aggregation step.
##' @param multi.core Logical indicating whether to attempt to parallelize the analysis on multiple cores. 
##' @param contrast.term If a fit object with multiple coefficients is passed in, a string indiating the coefficient of interest.   
##' @param scoring The gRNA ranking method to use in RRAa aggregation. May take one of three values: \code{pvalue}, \code{fc},
##' or '\code{combined}'. \code{pvalue} indicates that the gRNA ranking statistic should be created from the (one-sided) p-values in the 
##' fit object. \code{fc} indicates that the ranks of the gRNA coefficients should be used instead, and \code{combined} indicates that 
##' that the coefficents should be used as the ranking statistic but gRNAs are discarded in the aggregation step based on the corresponding nominal 
##' p-value in the fit object. 
##' @param permutation.seed numeric seed for permutation reproducibility. 
##'   Default: \code{NULL} means to not set any seed. This argument is passed
##'   through to \code{\link{ct.RRAaPvals}}.
##' @return A dataframe containing gRNA-level and target-level statistics. In addition to the information present in the supplied annotation object, 
##' the returned object indicates P-values and Q-values for the depletion and enrichment of each gRNA and associated target, the median log2 fold 
##' change estimate among all gRNAs associated with the target, and Rho statistics that are calculated internally by the RRAa algorithm that may be 
##' useful in ranking targets that are considered significant at a given alpha or false discovery threshold.
##' @author Russell Bainer
##' @examples data('fit')
##' data('ann')
##' output <- ct.generateResults(fit, ann, permutations = 10)
##' head(output)
##' @export

ct.generateResults <- function(fit,
                               annotation,
                               RRAalphaCutoff = 0.1,
                               permutations = 1000,
                               multi.core = TRUE,
                               contrast.term = NULL,
                               scoring = c("combined", "pvalue", "fc"),
                               permutation.seed = NULL) {
  
  #figure out the scoring method
  methods <- c("combined", "pvalue", "fc")
  scoring <- match.arg(scoring, methods, several.ok = FALSE)
  
  #make sure that the Fit has P-values. 
  if(!('p.value' %in% names(fit))){
    stop('Evidence for differential gRNA abundance (p-values) has not been estimated for the provided fit object. Please do so using eBayes(), treat(), or similar method.')
  }
  
  #Check/omit the extra coefficient columns as necessary. 
  if(ncol(fit$coefficients) > 1){
    if(is.null(contrast.term)){
      stop("The fit object contains multiple coefficients. Please specify a contrast.term.")
    }
    fit <- ct.preprocessFit(fit, contrast.term)
  }
 

  #The code below is now handled by the prepareAnnotation function but is retained here for legacy purposes. 
  #if(!setequal(row.names(fit), row.names(annotation))){
  #  if(length(setdiff(row.names(fit), row.names(annotation))) > 0){
  #    stop("fit contains elements not present in the annotation file.")
  #    }
  #  message("The annotation file contains elements not present in the fit object. They will be discarded for downstream analyses.")
  #  annotation <- annotation[row.names(fit),]
  #}
  key <- ct.prepareAnnotation(annotation, fit, throw.error = FALSE)

  #Prepare the ranking values and calculate p-values
  pvals <- ct.DirectionalTests(fit)
  foldchange <- cbind(fit$coefficients[,1], -fit$coefficients[,1]) 

  if(scoring %in% "fc"){
    #Normalize values to rank scores
    scores.deplete <- as.matrix(rank(foldchange[,1])/nrow(foldchange))
    scores.enrich <- as.matrix(rank(foldchange[,2])/nrow(foldchange))

    #determine the fold-change significance cutoffs
    cut.deplete <- sort(scores.deplete[,1])[round(nrow(scores.deplete) * RRAalphaCutoff)]
    cut.enrich <- sort(scores.enrich[,1])[round(nrow(scores.enrich) * RRAalphaCutoff)]    
    
    }else if(scoring %in% "pvalue"){
          
      #Normalize values to rank scores
      scores.deplete <- as.matrix(rank(pvals[,1])/nrow(pvals))
      scores.enrich <- as.matrix(rank(pvals[,2])/nrow(pvals))
      
      #determine the significance cutoffs for the rank statistics on the basis of p-values
      cut.deplete <- sort(scores.deplete[,1])[sum(pvals[,1] <= 0.05, na.rm = TRUE)]
      cut.enrich <- sort(scores.enrich[,1])[sum(pvals[,2] <= 0.05, na.rm = TRUE)]
          
    } else {
  
      cut.deplete <- (pvals[,1] <= RRAalphaCutoff)
      is.na(cut.deplete) <- FALSE
      scores.deplete <- as.matrix(rank(foldchange[,1])/nrow(foldchange))
      #scores.deplete[!cut.deplete,] <- 1     #set the nonsignificant scores to 1. 
      
      cut.enrich <- (pvals[,2] <= RRAalphaCutoff)
      is.na(cut.enrich) <- FALSE
      scores.enrich <- as.matrix(rank(foldchange[,2])/nrow(foldchange))
      #scores.enrich[!cut.enrich,] <- 1     #set the nonsignificant scores to 1. 
  
      }

  geneP.depletion <-
    ct.RRAaPvals(
      scores.deplete,
      g.key = key,
      alpha = cut.deplete,
      multicore = multi.core,
      permute = permutations,
      core.perm = 100,
      permutation.seed = permutation.seed
    )
  geneQ.depletion <- p.adjust(geneP.depletion, method = "fdr")

  geneP.enrichment <-
    ct.RRAaPvals(
      scores.enrich,
      g.key = key,
      alpha = cut.enrich,
      multicore = multi.core,
      permute = permutations,
      core.perm = 100,
      permutation.seed = permutation.seed
    )
  geneQ.enrichment <- p.adjust(geneP.enrichment, method = "fdr")

  #generate the Rho Ranks: 
  rhoEnrich <- ct.RRAalpha(scores.enrich, 
                           g.key = key, 
                           alpha = cut.enrich, 
                           shuffle = FALSE, 
                           return.obj = TRUE)
  #rhoEnrich <- rank(rho)
  #names(rhoEnrich) <- names(rho)
  
  rhoDeplete <- ct.RRAalpha(scores.deplete, 
                           g.key = key, 
                           alpha = cut.deplete, 
                           shuffle = FALSE, 
                           return.obj = TRUE)
  #rhoDeplete <- rank(rho)
  #names(rhoDeplete) <- names(rho)
  
  annotFields <- c("ID", "target", "geneID", "geneSymbol")  
  if(!all(annotFields %in% names(key))){
    message(paste("Some expected columns are not present in the supplied annotation file.", call. = FALSE))
    annotFields <- intersect(annotFields, names(key))
    message(paste("Only the following information will be included in the output:", paste(annotFields, collapse = ',')))  
  } 

  #make the DF
  summaryDF <- key[,annotFields]
  summaryDF$geneSymbol <- as.character(summaryDF$geneSymbol)
  summaryDF["gRNA Log2 Fold Change"] <- fit$coefficients[row.names(summaryDF),1]
  summaryDF["gRNA Depletion P"] <- signif(pvals[row.names(summaryDF),1], 5)
  summaryDF["gRNA Depletion Q"] <- signif(p.adjust(pvals[,1], "fdr")[row.names(summaryDF)], 5)
  summaryDF["gRNA Enrichment P"] <- signif(pvals[row.names(summaryDF),2], 5)
  summaryDF["gRNA Enrichment Q"] <- signif(p.adjust(pvals[,2], "fdr")[row.names(summaryDF)], 5)
  summaryDF["Target-level Enrichment P"] <- geneP.enrichment[summaryDF$geneSymbol]
  summaryDF["Target-level Enrichment Q"] <- geneQ.enrichment[summaryDF$geneSymbol]
  summaryDF["Target-level Depletion P"] <- geneP.depletion[summaryDF$geneSymbol]
  summaryDF["Target-level Depletion Q"] <- geneQ.depletion[summaryDF$geneSymbol]
  
  #Add a column for the median FC for each target: 
  targets <- unique(summaryDF$geneSymbol)
  medianfc <- vapply(targets, 
                     function(x){median(summaryDF[(summaryDF$geneSymbol == x),"gRNA Log2 Fold Change"], 
                                        na.rm = TRUE)}, 
                     numeric(1), 
                     USE.NAMES = TRUE)
  summaryDF["Median log2 Fold Change"] <- vapply(summaryDF$geneSymbol, 
                                                 function(x){medianfc[x][1]}, 
                                                 numeric(1))
  summaryDF["Rho_enrich"] <-rhoEnrich[summaryDF$geneSymbol]
  summaryDF["Rho_deplete"] <- rhoDeplete[summaryDF$geneSymbol]
  
  #order them
  summaryDF <- summaryDF[order(summaryDF[,"Rho_enrich"], decreasing = FALSE),]
  summaryDF <- summaryDF[order(summaryDF[,"Target-level Enrichment P"], decreasing = FALSE),]
  return(summaryDF)
}
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



