##' @title Calculate results of a crispr screen from a contrast
##' @description This is a wrapper function that enables direct generation of target-level p-values from a crispr screen.  
##' @param fit An object of class \code{MArrayLM} containing, at minimum, a \code{t} slot with t-statistics from the comparison, 
##' a \code{df.residual} slot with the corresponding residuals fo the model fits, and an \code{Amean} slot with the respective mean abundances. 
##' @param annotation An annotation file for the experiment. gRNAs are annotated by 
##' row, and must minimally contain columns \code{geneSymbol} and \code{geneID}.
##' @param RRAalphaCutoff A cutoff to use when defining gRNAs with significantly altered abundance during the RRAa aggregation step, which may be specified
##' as a single numeric value on the unit interval or as a logical vector. When supplied as a logical vector (of length equal to \code{nrows(fit)}), 
##' this parameter directly indicates the gRNAs to include during RRAa aggregation. Otherwise, if \code{scoring} is set 
##' to \code{pvalue} or \code{combined}, this parameter is interpreted as the maximum nominal p-value required to consider a gRNA's abundance meaningfully 
##' altered during the aggregation step. If \code{scoring} is \code{fc}, this parameter is interpreted as the proportion of the list to be considered 
##' meaningfully altered in the experiment (e.g., if \code{RRAalphaCutoff} is set to 0.05, only consider the rankings of the 5% most upregulated 
##' (or downregulated) gRNAs for the purposes of RRAa calculations).
##' 
##' Note that this function uses directional tests to identify enriched or depleted targets, and when RRAalphaCutoff is 
##' provided as a logical vector the interpretation of the various aggregation statistics is going to be dependent on the specific criteria 
##' used to select reagents for inclusion. 
##' 
##' @param permutations The number of permutations to use during the RRAa aggregation step.
##' @param contrast.term If a fit object with multiple coefficients is passed in, a string indiating the coefficient of interest.   
##' @param scoring The gRNA ranking method to use in RRAa aggregation. May take one of three values: \code{pvalue}, \code{fc},
##' or '\code{combined}'. \code{pvalue} indicates that the gRNA ranking statistic should be created from the (one-sided) p-values in the 
##' fit object. \code{fc} indicates that the ranks of the gRNA coefficients should be used instead, and \code{combined} indicates that 
##' that the coefficents should be used as the ranking statistic but gRNAs are discarded in the aggregation step based on the corresponding nominal 
##' p-value in the fit object. 
##' @param alt.annotation Libraries targeting ambiguous biological elements (e.g., alternative promoters to a gene where the boundaries between
##' elelments is contested) may contain reagents that are plausibly annotated to a finite set of possible targets. To accomodate this, users may 
##' supply an alternative reagent annotation in the form of a named list of vectors, where each list element corresponds something coercible to a to
##' a character vector of associated targets that will ultimately be assembled into the `geneSymbol` column of the `resultsDF` object. Each of these 
##' character vectors should be named identically to a row of the supplied fit object (e.g., the `row.names`). It is assumed that the `geneID` values 
##' are assigned unambiguously to the reagents, and are passed through directly. 
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
##' @return A `resultsDF` formatted dataframe containing gene-level statistics.
##' @examples
##'   p = seq(0, 1, length.out=20)
##'   fc = seq(-3, 3, length.out=20)
##'   fc[2] = NA
##'   fc[3] = -20
##'   stats = data.frame(
##'     Depletion.P=p,
##'     Enrichment.P=rev(p),
##'     fc=fc
##'   )
##'   ct.applyAlpha(stats,scoring="combined")
##' @export 
ct.generateResults <- function(fit,
                               annotation,
                               RRAalphaCutoff = 0.1,
                               permutations = 1000,
                               contrast.term = NULL,
                               scoring = c("combined", "pvalue", "fc"),
                               alt.annotation = NULL,
                               permutation.seed = NULL) {

    ## figure out the scoring method
    scoring <- match.arg(scoring)

    #make sure that the Fit has P-values. 
    if(!('p.value' %in% names(fit))){
      stop('Evidence for differential gRNA abundance (p-values) has not been estimated for the provided fit object. Please do so using eBayes(), treat(), or similar method.')
    }
 
    #Select relevant contrast
    if(ncol(fit$coefficients) > 1){
        if(is.null(contrast.term)){
            stop("The fit object contains multiple coefficients. Please specify a contrast.term.")
        }
        fit <- ct.preprocessFit(fit, contrast.term)
    }
    
    #Apply directional test to the model estimates
    pvals <- ct.DirectionalTests(fit)
    
    #prep the annotation
    key <- ct.prepareAnnotation(annotation, fit, throw.error = FALSE)
    if(!('ID' %in% names(key))){key$ID <- row.names(key)}

    #If provided, check the alt.annotation and expand the fit/key as needed. 
    if(!is.null(alt.annotation)){
      stopifnot(is.list(alt.annotation), all(names(alt.annotation) %in% row.names(fit)))

      #Expand the annotation object
      key <- ct.expandAnnotation(key, alt.annotation)
    }

    ## Prepare the ranking values 
    rra.input <- ct.applyAlpha(cbind(pvals[key$ID,], FC=fit$coefficients[key$ID,1]), RRAalphaCutoff, scoring )
    
    
    ## Add pvalues and qvalues
    geneP.depletion <-
        ct.RRAaPvals(
            rra.input[,"scores.deplete", drop=FALSE],
            g.key = key,
            permute = permutations,
            permutation.seed = permutation.seed
        )
    
    geneP.enrichment <-
        ct.RRAaPvals(
            rra.input[,"scores.enrich", drop=FALSE],
            g.key = key,
            permute = permutations,
            permutation.seed = permutation.seed
        )

    ## generate the Rho values and ranks:
    ## These are redundant with work done in ct.RRAPvals ...
    rhoEnrich <- ct.RRAalpha(rra.input[,"scores.enrich", drop=FALSE],
                             g.key = key, 
                             shuffle = FALSE)
    
    rhoDeplete <- ct.RRAalpha(rra.input[,"scores.deplete", drop=FALSE], 
                              g.key = key, 
                              shuffle = FALSE)
    
    annotFields <- c("ID", "target", "geneID", "geneSymbol")  
    if(!all(annotFields %in% names(key))){
        message(paste("Some expected columns are not present in the supplied annotation file.", call. = FALSE))
        annotFields <- intersect(annotFields, names(key))
        message(paste("Only the following information will be included in the output:", paste(annotFields, collapse = ',')))  
    } 

    ## make the DF
    summaryDF <- key[,annotFields]
    summaryDF$geneSymbol <- as.character(summaryDF$geneSymbol)
    summaryDF["gRNA Log2 Fold Change"] <- fit$coefficients[summaryDF$ID,1]
    summaryDF["gRNA Depletion P"] <- signif(pvals[summaryDF$ID,1], 5)
    summaryDF["gRNA Depletion Q"] <- signif(p.adjust(pvals[,'Depletion.P'], "fdr")[summaryDF$ID], 5)
    summaryDF["gRNA Enrichment P"] <- signif(pvals[summaryDF$ID,2], 5)
    summaryDF["gRNA Enrichment Q"] <- signif(p.adjust(pvals[,'Enrichment.P'], "fdr")[summaryDF$ID], 5)
    summaryDF["Target-level Enrichment P"] <- geneP.enrichment[summaryDF$geneSymbol]
    summaryDF["Target-level Enrichment Q"] <- p.adjust(geneP.enrichment,"fdr")[summaryDF$geneSymbol]
    summaryDF["Target-level Depletion P"] <- geneP.depletion[summaryDF$geneSymbol]
    summaryDF["Target-level Depletion Q"] <- p.adjust(geneP.depletion,"fdr")[summaryDF$geneSymbol]

    ## Add a column for the median FC for each target:
    medianfc <- tapply(summaryDF[,"gRNA Log2 Fold Change"], summaryDF[,"geneSymbol"], median, na.rm=TRUE)
    summaryDF["Median log2 Fold Change"] <- as.numeric(medianfc[ summaryDF$geneSymbol ]) # One per gene rep'd out to one per guide
    summaryDF["Rho_enrich"] <-rhoEnrich[summaryDF$geneSymbol]
    summaryDF["Rho_deplete"] <- rhoDeplete[summaryDF$geneSymbol]
    
    ## order them
    #summaryDF <- summaryDF[order(summaryDF[,"Rho_enrich"], decreasing = FALSE),]
    #summaryDF <- summaryDF[order(summaryDF[,"Target-level Enrichment P"], decreasing = FALSE),]
    return(summaryDF)
}

##' Apply RRA 'alpha' cutoff to RRAalpha input
##'
##' The 'alpha' part of RRAalpha is used to consider only the top guide-level scores for gene-level
##'     statistics. Practically, all guides failing the cutoff get a pvalue of 1.  There are three ways of
##'     determining which guides fail. See 'scoring' below.
##' @param stats three-column numeric matrix with pvalues for down and up one-sided test with guide-level fold
##'     changes (coefficients from the relevant contrast).
##' @param RRAalphaCutoff A cutoff to use when defining gRNAs with significantly altered abundance during the RRAa aggregation step, which may be specified
##' as a single numeric value on the unit interval or as a logical vector. When supplied as a logical vector (of length equal to \code{nrows(fit)}), 
##' this parameter directly indicates the gRNAs to include during RRAa aggregation. Otherwise, if \code{scoring} is set 
##' to \code{pvalue} or \code{combined}, this parameter is interpreted as the maximum nominal p-value required to consider a gRNA's abundance meaningfully 
##' altered during the aggregation step. If \code{scoring} is \code{fc}, this parameter is interpreted as the proportion of the list to be considered 
##' meaningfully altered in the experiment (e.g., if \code{RRAalphaCutoff} is set to 0.05, only consider the rankings of the 5% most upregulated 
##' (or downregulated) gRNAs for the purposes of RRAa calculations).
##' @param scoring The gRNA ranking method to use in RRAa aggregation. May take one of three values: \code{pvalue}, \code{fc},
##' or '\code{combined}'. \code{pvalue} indicates that the gRNA ranking statistic should be created from the (one-sided) p-values in the 
##' fit object. \code{fc} indicates that the ranks of the gRNA coefficients should be used instead, and \code{combined} indicates that 
##' that the coefficents should be used as the ranking statistic but gRNAs are discarded in the aggregation step based on the corresponding nominal 
##' p-value in the fit object. 
##' @return data.frame with guide-level pvals, fold change, and scores.deplete and scores.enrich which are the input the RRAalpha
##' @examples 
##' fakestats <- matrix(runif(300), ncol = 3)
##' colnames(fakestats) = c('Depletion.P', 'Enrichment.P', 'lfc')
##' ct.applyAlpha(fakestats)
##' @author Russell Bainer
##' @export 
ct.applyAlpha <- function(stats, RRAalphaCutoff=0.1, scoring = c("combined", "pvalue", "fc")) {
    scoring <- match.arg(scoring)
    pvals = stats[,1:2]
    foldchange <- cbind(stats[,3], -stats[,3])

    ##RRAalphaCutoff format 
    if (length(RRAalphaCutoff) == 1) {
        rra.logic <- FALSE    
        if(!is.numeric(RRAalphaCutoff) | (RRAalphaCutoff < 0) | (RRAalphaCutoff > 1)){
            stop('When provided as a single value, RRAalphaCutoff must be a numeric value equal to 0, 1, or something in between.')
        }
    } else {
        rra.logic <- TRUE    
        if((sum(RRAalphaCutoff %in% c(TRUE, FALSE)) != nrow(pvals)) | (length(RRAalphaCutoff) != nrow(pvals))){
            stop('When provided as a vector, RRAalphaCutoff must be the exact length of the p-value vectors in the fit object and must only contain TRUE or FALSE values.')
        }
    }
    
    if (scoring %in% "fc") {
        ##Normalize values to rank scores
        scores.deplete <- as.matrix(rank(foldchange[,1])/nrow(foldchange))
        scores.enrich <- as.matrix(rank(foldchange[,2])/nrow(foldchange))

        ##determine the fold-change significance cutoffs
        cut.deplete <- sort(scores.deplete[,1])[round(nrow(scores.deplete) * RRAalphaCutoff)]
        cut.enrich <- sort(scores.enrich[,1])[round(nrow(scores.enrich) * RRAalphaCutoff)]    
        
        if(rra.logic){
            cut.deplete <- RRAalphaCutoff
            cut.enrich <- RRAalphaCutoff
        }
        
    } else if (scoring %in% "pvalue") {
        
        ## Normalize values to rank scores
        scores.deplete <- as.matrix(rank(pvals[,'Depletion.P'])/nrow(pvals))
        scores.enrich <- as.matrix(rank(pvals[,'Enrichment.P'])/nrow(pvals))
        
        ## determine the significance cutoffs for the rank statistics on the basis of p-values
        cut.deplete <- (pvals[,'Depletion.P'] <= RRAalphaCutoff)
        cut.enrich <- (pvals[,'Enrichment.P'] <= RRAalphaCutoff)
        
        if (rra.logic) {
            cut.deplete <- RRAalphaCutoff
            cut.enrich <- RRAalphaCutoff
        }
        
    } else {
        
        cut.deplete <- (pvals[,'Depletion.P'] <= RRAalphaCutoff)
        cut.deplete[ is.na(cut.deplete) ] <- FALSE
        scores.deplete <- as.matrix(rank(foldchange[,1])/nrow(foldchange))

        cut.enrich <- (pvals[,'Enrichment.P'] <= RRAalphaCutoff)
        cut.enrich[ is.na(cut.enrich) ] <- FALSE
        scores.enrich <- as.matrix(rank(foldchange[,2])/nrow(foldchange))

        if (rra.logic) {
            cut.deplete <- RRAalphaCutoff
            cut.enrich <- RRAalphaCutoff
        }
    }
    ## nonsignificant gRNAs are set to 1. Implicitly, this means that the significance 
    ## relative to alpha is treated as an inherent property of the gRNA.

    if(length(cut.deplete) == 1)
        pass <- scores.deplete[,1] <= cut.deplete
    else
        pass <- cut.deplete
    scores.deplete[!pass,1] <- 1
    if(length(cut.enrich) == 1)
        pass <- scores.enrich[,1] <= cut.enrich
    else
        pass <- cut.enrich
    scores.enrich[!pass,1] <- 1
    out = cbind(stats, scores.deplete=scores.deplete[,1], scores.enrich=scores.enrich[,1])
    return(out)
}
