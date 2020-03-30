##' @title Preprocess a "MArrayLM" model fit object to include only one contrast. 
##' @description This function preprocesses a fit object returned from eBayes to include only the values relevant to the
##' \code{modelTerm} specified.  
##' @param fit An object of class MArrayLM to be processed. 
##' @param modelTerm The model coefficient to be isolated for downstream analyses. 
##' @return A \code{MArrayLM} object for downstream processing. 
##' @author Russell Bainer
##' @import limma
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

ct.preprocessFit <- function(fit, modelTerm){
  if(!methods::is(fit, "MArrayLM")){stop(paste(deparse(substitute(fit)), "is not an MArrayLM object."))}
  if(!(modelTerm %in% colnames(fit$coefficients))){stop("Specified coefficient is not present in the fit object.")}
  if(!('p.value' %in% names(fit))){warning(paste(deparse(substitute(fit)), ' does not contain p-values quantifying the evidence for differential gRNA abundance. Eventually, you will need to process it with eBayes(), treat(), or a similar function.'))}
  
  fit$coefficients <- as.matrix(fit$coefficients[,modelTerm])
  colnames(fit$coefficients) <- modelTerm
  fit$stdev.unscaled <- as.matrix(fit$stdev.unscaled[,modelTerm])
  colnames(fit$stdev.unscaled) <- modelTerm
  fit$t <- as.matrix(fit$t[,modelTerm])
  colnames(fit$t) <- modelTerm
  fit$p.value <- as.matrix(fit$p.value[,modelTerm])
  colnames(fit$p.value) <- modelTerm
  if ("lods" %in% names(fit)) {
      fit$lods <- as.matrix(fit$lods[,modelTerm])
      colnames(fit$lods) <- modelTerm
  }
  
  return(fit)
}















