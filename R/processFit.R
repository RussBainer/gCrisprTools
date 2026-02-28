##' @title Preprocess a 'MArrayLM' model fit object to include only one contrast. 
##' @description This function preprocesses a fit object returned from eBayes to include only the values relevant to the
##' \code{modelTerm} specified.  
##' @param fit An object of class MArrayLM to be processed. 
##' @param modelTerm The model coefficient to be isolated for downstream analyses. 
##' @param fill.missing Logical indicating whether missing elements of `fit` (usually due to incomplete or unreplicated data)
##' should be created with all elements set to zero. This is not recommended, but sometimes is necessary to use specific 
##' `gCrisprTools` QC or visualization tooling.  
##' @return A \code{MArrayLM} object for downstream processing. 
##' @author Russell Bainer
##' @import limma
##' @importFrom methods is
##' @keywords internal
##' @examples 
##' 
##' #Load and preprocess data
##' data('es')
##' library(Biobase)
##' library(limma)
##' 
##' #Make a multi-level contrast
##' design <- model.matrix(~ 0 + TREATMENT_NAME, pData(es))
##' colnames(design) <- gsub('TREATMENT_NAME', '', colnames(design))
##' contrasts <- makeContrasts((ControlExpansion - ControlReference), (DeathExpansion - ControlExpansion), levels = design)
##' 
##' #Make a multi-level fit object
##' vm <- voom(exprs(es), design)
##' fit <- lmFit(vm, design)
##' fit <- contrasts.fit(fit, contrasts)
##' fit <- eBayes(fit)  
##' 
##' #And trim it
##' fit2  <- ct.preprocessFit(fit, modelTerm = '(DeathExpansion - ControlExpansion)')
##' 
##' ncol(fit)
##' ncol(fit2)
##' @export 

ct.preprocessFit <- function(fit, modelTerm, fill.missing = FALSE) {
    stopifnot(methods::is(fill.missing, 'logical'), length(fill.missing) == 1)
  
    if (!methods::is(fit, "MArrayLM")) {
        stop(deparse(substitute(fit)), " is not an MArrayLM object.")
    }
    if (!(modelTerm %in% colnames(fit$coefficients))) {
        stop("Specified coefficient is not present in the fit object.")
    }
  
    missingCols <- setdiff(c('p.value', 't', 'stdev.unscaled', 'lods'), names(fit))
    if (!("p.value" %in% missingCols)) {
      warning(deparse(substitute(fit)), " does not contain p-values quantifying the evidence for differential gRNA abundance. Eventually, you will need to process it with eBayes(), treat(), or a similar function.")
    }
    
    if(length(missingCols) > 0){
      if(!fill.missing){
        stop("Columns missing from ", deparse(substitute(fit)), "!\nIf you want me to make something up set fill.missing=TRUE.")
      } else {
        newset <- fit$coefficients
        newset[,] <- 0
        for(missingCol in missingCols){
          fit[[missingCol]] <- newset
          }
        }
      warning(deparse(substitute(fit)), " is missing the following expected columns: ", paste0(missingCols, collapse = ', '))
      }

    fit$coefficients <- as.matrix(fit$coefficients[, modelTerm])
    colnames(fit$coefficients) <- modelTerm
    fit$stdev.unscaled <- as.matrix(fit$stdev.unscaled[, modelTerm])
    colnames(fit$stdev.unscaled) <- modelTerm
    fit$t <- as.matrix(fit$t[, modelTerm])
    colnames(fit$t) <- modelTerm
    fit$p.value <- as.matrix(fit$p.value[, modelTerm])
    colnames(fit$p.value) <- modelTerm
    if ("lods" %in% names(fit)) {
        fit$lods <- as.matrix(fit$lods[, modelTerm])
        colnames(fit$lods) <- modelTerm
    }

    return(fit)
}















