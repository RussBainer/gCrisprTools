##' gCrisprTools
##'
##' Pipeline for using CRISPR data at Genentech
##'
##' @docType package
##' @name gCrisprTools-package
NULL

##' @import Biobase
##' @import limma
##' @import RobustRankAggreg
##' @import ggplot2
##' @import parallel
##' @import BiocParallel
##' @import PANTHER.db
##' @importFrom grDevices colorRampPalette rgb
##' @importFrom graphics abline axis barplot layout legend lines mtext par plot points polygon segments hist
##' @importFrom stats density lm median na.omit p.adjust pbeta phyper predict pt smooth.spline spline
##' @importFrom utils capture.output getFromNamespace 
##' @importFrom rmarkdown render draft 
NULL



##' @name es
##' @aliases es
##' @docType data
##' @title ExpressionSet of count data from a Crispr screen with strong selection
##' @description
##' Expressionset of raw counts from a screen in mouse cells performed at Genentech, Inc. 
##' All sample, gRNA, and Gene information has been anonymized and randomized. 
##' @examples
##' data("es")
##' print(es)
##' @seealso Please see \file{vignettes/Crispr_example_workflow.R} for details.
##' @source Genentech, Inc.
NULL

##' @name ann
##' @aliases ann
##' @docType data
##' @title Annotation file for a mouse Crispr library
##' @description
##' Example annotation file for the screen data provided in \code{es}. 
##' All sample, gRNA, and Gene information has been anonymized and randomized. 
##' @examples
##' data("ann")
##' head(ann)
##' @seealso Please see \file{vignettes/Crispr_example_workflow.R} for details.
##' @source Genentech, Inc.
NULL

##' @name fit
##' @aliases fit
##' @docType data
##' @title Precalculated contrast fit from a Crispr screen 
##' @description A precalculated fit object (class \code{MArrayLM}) comparing the death 
##' and control expansion arms of a crispr screen performed at Genentech, Inc. 
##' All sample, gRNA, and Gene information has been anonymized and randomized. 
##' @examples
##' data("fit")
##' show(fit)
##' @seealso Please see \file{vignettes/Crispr_example_workflow.R} for model details.
##' @source Genentech, Inc.
NULL

##' @name resultsDF
##' @aliases resultsDF
##' @docType data
##' @title Precalculated gene-level summary of a crispr screen 
##' @description A precalculated summary Dataframe comparing the death and control expansion arms of 
##' the provided example Crispr screen (using 8 cores, seed = 2). 
##' All sample, gRNA, and Gene information has been anonymized and randomized. 
##' @examples
##' data("resultsDF")
##' head(resultsDF)
##' @seealso Please see \file{vignettes/Crispr_example_workflow.R} for model details.
##' @source Genentech, Inc.
NULL

##' @name aln
##' @aliases aln
##' @docType data
##' @title Precalculated alignment statistics of a crispr screen 
##' @description Example alignment matrix file for the provided example Crispr screen.
##' All sample, gRNA, and Gene information has been anonymized and randomized. 
##' @examples
##' data("aln")
##' head(aln)
##' @seealso Please see \file{vignettes/Crispr_example_workflow.R} for details.
##' @source Genentech, Inc.
NULL

##' @name essential.genes
##' @aliases essential.genes
##' @docType data
##' @title Artificial list of 'essential' genes in the example Crispr screen 
##' included for plotting purposes 
##' @description Example gene list, designed to demonstrate ROC and PRC functions
##' All sample, gRNA, and Gene information has been anonymized and randomized. 
##' @examples
##' data("essential.genes")
##' essential.genes
##' @seealso Please see \file{vignettes/Crispr_example_workflow.R} for details.
##' @source Russell Bainer
NULL


