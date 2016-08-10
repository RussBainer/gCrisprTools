## ---- eval = FALSE-------------------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("gCrisprTools")

## ---- message=FALSE, warning=FALSE---------------------------------------
library(Biobase)
library(limma)
library(gCrisprTools)

## ------------------------------------------------------------------------
data("es", package = "gCrisprTools")
es
head(exprs(es))

## ------------------------------------------------------------------------
data("ann", package = "gCrisprTools")
head(ann)

## ------------------------------------------------------------------------
sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
names(sk) <- row.names(pData(es))
sk

## ---- fig.width = 7, fig.height = 5--------------------------------------
data("aln", package = "gCrisprTools")
head(aln)
ct.alignmentChart(aln, sk)

## ---- fig.width=6, fig.height = 8----------------------------------------
es.floor <- ct.filterReads(es, read.floor = 30, sampleKey = sk)
es <- ct.filterReads(es, trim = 1000, log2.ratio = 4, sampleKey = sk)

##Convenience function for conforming the annotation object to exclude the trimmed gRNAs
ann <- ct.prepareAnnotation(ann, es, controls = "NoTarget")

## ---- fig.width=6, fig.height = 8----------------------------------------
es <- ct.normalizeGuides(es, 'scale', annotation = ann, sampleKey = sk, plot.it = TRUE)
es.norm <- ct.normalizeGuides(es, 'slope', annotation = ann, sampleKey = sk, plot.it = TRUE)
es.norm <- ct.normalizeGuides(es, 'controlScale', annotation = ann, sampleKey = sk, plot.it = TRUE, geneSymb = 'NoTarget')
es.norm <- ct.normalizeGuides(es, 'controlSpline', annotation = ann, sampleKey = sk, plot.it = TRUE, geneSymb = 'NoTarget')

## ---- eval= FALSE--------------------------------------------------------
#  #Not run:
#  path2QC <- ct.makeQCReport(es,
#                   trim = 1000,
#                   log2.ratio = 0.05,
#                   sampleKey = sk,
#                   annotation = ann,
#                   aln = aln,
#                   identifier = 'Crispr_QC_Report',
#                   lib.size = NULL
#                   )

## ---- fig.width=6, fig.height=6------------------------------------------
ct.rawCountDensities(es, sk)

## ---- fig.width=6, fig.height = 6----------------------------------------
ct.gRNARankByReplicate(es, sk)  #Visualization of gRNA abundance distribution

## ---- fig.width=6, fig.height = 6----------------------------------------
ct.gRNARankByReplicate(es, sk, annotation = ann, geneSymb = "Target1633")

## ---- fig.width=6, fig.height = 6----------------------------------------
ct.viewControls(es, ann, sk, normalize = FALSE, geneSymb = 'NoTarget')

## ---- fig.width=6, fig.height = 4----------------------------------------
ct.guideCDF(es, sk, plotType = "gRNA")

## ------------------------------------------------------------------------
design <- model.matrix(~ 0 + sk, sk)
colnames(design) <- gsub('sk', '', colnames(design))
contrasts <- makeContrasts(DeathExpansion - ControlExpansion, levels = design)

vm <- voom(exprs(es), design)

fit <- lmFit(vm, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)

## ---- message=FALSE, warning=FALSE---------------------------------------
resultsDF <-
  ct.generateResults(
    fit,
    annotation = ann,
    RRAalphaCutoff = 0.1,
    permutations = 1000,
    scoring = "combined"
  )

## ---- fig.width=6, fig.height = 6----------------------------------------
ct.topTargets(fit,
              resultsDF,
              ann,
              targets = 10,
              enrich = FALSE)

## ---- fig.width=6, fig.height = 8----------------------------------------
ct.stackGuides(
  es,
  sk,
  plotType = "Target",
  annotation = ann,
  subset = names(sk)[grep('Expansion', sk)]
)

## ---- fig.width=6, fig.height = 4----------------------------------------
ct.viewGuides("Target1633", fit, ann)

## ---- eval= FALSE--------------------------------------------------------
#  #Not run:
#  path2Contrast <-
#    ct.makeContrastReport(eset = es,
#                          fit = fit,
#                          sampleKey = sk,
#                          results = resultsDF,
#                          annotation = ann,
#                          comparison.id = NULL,
#                          identifier = 'Crispr_Contrast_Report')

## ---- eval=FALSE---------------------------------------------------------
#  #Not run:
#  path2report <-
#    ct.makeReport(fit = fit,
#                  eset = es,
#                  sampleKey = sk,
#                  annotation = ann,
#                  results = resultsDF,
#                  aln = aln,
#                  outdir = ".")

## ---- eval= FALSE--------------------------------------------------------
#  #Not run:
#  enrichmentResults <-
#    ct.PantherPathwayEnrichment(
#      resultsDF,
#      pvalue.cutoff = 0.01,
#      enrich = TRUE,
#      organism = 'mouse'
#    )
#  
#  > head(enrichmentResults)   #Note: Pathway names have been edited for display purposes.
#                           PATHWAY nGenes sigGenes expected     odds            p        FDR
#  1 EGF receptor signaling pathway    200       14 5.240550 3.332647 0.0004498023 0.03958260
#  2          FGF signaling pathway    230       14 5.949958 2.869779 0.0016284304 0.07165094
#  3     Insulin/MAP kinase cascade    138        9 3.714916 2.785331 0.0101632822 0.20272465
#  4          CCKR signaling map ST    331       15 8.211061 2.148459 0.0126368744 0.20272465
#  5               p38 MAPK pathway    145        9 3.891333 2.641968 0.0135928688 0.20272465
#  6              B cell activation    204       11 5.336192 2.368705 0.0146442541 0.20272465

## ---- fig.width=6, fig.height = 6, warning=FALSE-------------------------
data("essential.genes", package = "gCrisprTools")  #Artificial list created for demonstration
data("resultsDF", package = "gCrisprTools")
ROC <- ct.ROC(resultsDF, essential.genes, stat = "deplete.p")
str(ROC)

## ---- fig.width=6, fig.height = 6, warning=FALSE-------------------------
PRC <- ct.PRC(resultsDF, essential.genes, stat = "deplete.p")
str(PRC)

## ---- fig.width=6, fig.height = 6, warning=FALSE-------------------------
targetsTest <- ct.targetSetEnrichment(resultsDF, essential.genes, enrich = FALSE)
str(targetsTest)

## ------------------------------------------------------------------------
sessionInfo()

