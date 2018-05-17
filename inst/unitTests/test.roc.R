library(gCrisprTools)
library(RUnit)

test.RocPRC <- function() {
  data('resultsDF')
  data('essential.genes')

  rhp <- c(1.000000000, 1.000000000, 1.000000000, 1.000000000, 1.000000000, 0.977446507, 0.959707507, 0.009076336)
  
  prc <- ct.PRC(resultsDF, essential.genes, 'enrich.p')
  roc <- ct.ROC(resultsDF, essential.genes, 'enrich.p')  
    
  checkIdentical(round(roc$P.values[,3], 9), unlist(rhp))
  checkIdentical(round(prc$P.values[,3], 9), unlist(rhp))

  checkIdentical(prc$recall, roc$sensitivity)
}






