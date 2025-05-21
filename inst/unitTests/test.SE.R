library(gCrisprTools)
library(RUnit)
library(SummarizedExperiment)
suppressPackageStartupMessages(library('limma'))
data('es')
data('ann')
data('fit')
data('resultsDF')
data('se')

design <- model.matrix(~ 0 + REPLICATE_POOL + TREATMENT_NAME, pData(es))
colnames(design) <- gsub('TREATMENT_NAME', '', colnames(design))
vm <- vm <- voom(exprs(es), design)

sk <- relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference")
names(sk) <- row.names(pData(es))
  
  
test.SEbuild <- function() {
  
  newse <- ct.buildSE(es = es, 
                      sampleKey = sk, 
                      ann = ann, 
                      vm = vm, 
                      fit = fit, 
                      summaryList = list('resA' = resultsDF, 
                                         'resB' = resultsDF))
  
    checkIdentical( newse, se )
    checkIdentical( es, ct.extractSE('es', newse) )
    checkIdentical( ann, ct.extractSE('ann', newse) )

    return(TRUE)
}
