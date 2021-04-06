##' Update a gene-centric msdb object for GREAT-style enrichment analysis using a specified CRISPR annotation.
##'
##' @param annotation an annotation object returned by \code{ct.prepareAnnotation()}. 
##' @param gsdb A gene-centric \code{GeneSetDb} object to conform to the relevant peakwise dataset.
##' @param minsize Minimum number of targets required to consider a geneset valid for analysis.
##' @return A new \code{GeneSetDb} object with the features annotated genewise to pathways.
ct.GREATdb <- function(annotation, 
                       gsdb = getMSigGeneSetDb(c('h', 'c2'), 'human'), 
                       minsize = 10){
  
  #check the anbnotation object
  annotation <- ct.prepareAnnotation(annotation)
  
  #Check to make sure that the geneID column values are present in the provided gsdb
  genes <-na.omit(unique(annotation$geneID))
  message(paste0(length(intersect(multiGSEA::featureIds(gsdb), genes)), ' valid genes found in the supplied gene database (of ', length(genes), ')'))

  #Tabulate targets per gene
  targetsPerGene <- lapply(genes,
                         function(x){
                           return(as.character(unique(annotation$geneSymbol[annotation$geneID %in% x])))
                         })
  names(targetsPerGene) <- genes
  
  #Pull the lists from the gsdb
  gsdb_list <- as.list(gsdb)
  
  new_gs <- lapply(gsdb_list,
                   function(x){
                     return(as.character(unique(unlist(targetsPerGene[x]))))
                   })
  new_gs <- new_gs[unlist(lapply(new_gs, length)) >= minsize]
  message(paste0(length(new_gs), ' adjusted genesets created.'))
  if(length(new_gs) == 0){stop('Exiting.')}
  
  collections <- vapply(names(new_gs), function(x){strsplit(x, split = ';;')[[1]][1]}, character(1))
  sets <- vapply(names(new_gs), function(x){strsplit(x, split = ';;')[[1]][2]}, character(1))
  names(new_gs) <- sets
  
  out <- lapply(unique(collections),
                function(z){
                  new_gs[collections %in% z]
                })
  names(out) <- unique(collections)
  return(multiGSEA:::GeneSetDb.list(out))
}

##' Prepare a resultsDF object for analysis via MultiGSEA
##'
##' @param dflist A list of gCrisprTools results `data.frames` to be formatted.
##' @param cutoff Numeric maximum value of `statistic` to define significance.
##' @param statistic Should cutoffs be calculated based on FDR (`best.q`) or P-value (`best.p`)?
##' @param gdb Optionally, a `GeneSetDb` object to enable proper registration of the output. If provided, the 
##' collapsing features in the provided `simpleDF`s must be present in the `gsd@db$feature_id` slot. 
##' @param ... Other parameters to lower functions, especially `ct.simpleResult()`.
##' @return A list of `data.frames` formatted for evaluation with `multiGSEA`. 
##' @export
ct.mgseaPrep <- function(dflist, 
                         cutoff = 0.1, 
                         statistic = c('best.q', 'best.p'), 
                         gdb = NULL,
                         ...){
  #Input check
  dflist <- ct.regularizeContrasts(dflist, ...)
  stopifnot(is(cutoff, 'numeric'), cutoff <= 1, cutoff >= 0)
  statistic <- match.arg(statistic)
  
  if(!is.null(gsd)){
    dflist <- sapply(dflist, 
                     function(x){
                       x[(row.names(x) %in% gsd@db$feature_id),]
                     }, simplify = FALSE)
  }
  
  out <- sapply(dflist, 
                function(x){
                  data.frame('feature_id' = x$geneID, 
                             'logFC' = 1,
                             'selected' = (x[,statistic] <= cutoff),
                             'direction' = x$direction,
                             stringsAsFactors = FALSE)
                },simplify = FALSE)

  names(out) <- names(dflist)
  return(out)
  
}
                         
##' @title Geneset Enrichment within a CRISPR screen using multiGSEA
##' 
##' This function identifies differentially enriched/depleted ontological categories within the hits of a CRISPR screen 
##' given a provided `GenseSetDb()` and a results `data.frame` created by `ct.generateResults()`. Testing is performed using 
##' a Hypergeometric test, and results are returned as a `MultiGSEAResult` object defined in the `multiGSEA` package. Note that the 
##' `@logFC` slot in the returned object will contain the median gRNA lfc across all associated guides, which in some cases may 
##' have dubious interpretive value. 
##' 
##' This method used overrepresentation analysis, derived from `limma::kegga()`, and incorporates the number of gRNAs associated with 
##' each Target (inferred from the `geneSymbol` column of the `resultsDF`)  as the bias vector (because standard aggregation methods 
##' should be underpowered for targets with few guides). Setting `unbiased` = `TRUE` suppresses this behavior, which is identical to 
##' a hypergeometric test.  
##' 
##' @param resultsDF `data.frame` returned by `ct.generateResults()`.
##' @param gsdb `GenseSetDb` object containing annotations.
##' @param stat Statistic to be used in calling enrichment/depletion in the screen. Must be one of 'q', 'p', or 'rho'. 
##' @param cutoff Q, P, or Rho statistic cutoff defining significant enrichment/depletion in the screen. Default is 0.1. 
##' @param unbiased Logical indicating whether to estimate bias on the basis of the number of gRNAs associated with each target.
##' @examples 
##' data('ann')
##' data('resultsDF')
##' #gsd <- multiGSEA::getMSigGeneSetDb(c('h', 'c2'), 'mouse', id.type = 'entrez')
##' 
##' @author Steve Lianoglou for multiGSEA; Russell Bainer for wrapping functions.
ct.multiGSEA <- function(resultsDF, 
                         gsdb,
                         cutoff = 0.1, 
                         stat = c('q', 'p', 'rho'),
                         unbiased = FALSE){
  #input check
  if (!is(gsdb, "GeneSetDb")) {
    if (is(gsdb, "GeneSetCollection") || is(gsdb, "GeneSet")) {
      stop("A GeneSetDb is required. GeneSetCollections can be converted to a ", 
           "GeneSetDb via the GeneSetDb constructor, or a call to ", 
           "`as(gene_set_collection, 'GeneSetDb')`. See ?GeneSetDb for more ", 
           "details.")
    }
    stop("GeneSetDb required. Please see `?GeneSetDb` for ways to turn your ", 
         "gene sets into a GeneSetDb")
  }
  stopifnot(ct.resultCheck(resultsDF))
  
  guides <- table(resultsDF$geneSymbol)
  resultsDF <- resultsDF[!duplicated(resultsDF$geneSymbol),]
  stopifnot(is.numeric(cutoff), cutoff <= 1, cutoff >= 0)

  stat <- match.arg(stat)
  id <- match.arg(id)
  
  #Compose the df for feeding multigsea
  
  oraIpt <- data.frame('feature_id' = as.character(resultsDF$geneID), 
                       'significant' = switch(stat,
                                      'q' = (resultsDF$`Target-level Enrichment Q` < cutoff) | (resultsDF$`Target-level Depletion Q` < cutoff), 
                                      'p' = (resultsDF$`Target-level Enrichment P` < cutoff) | (resultsDF$`Target-level Depletion P` < cutoff), 
                                      'rho' = (resultsDF$Rho_enrich < cutoff) | (resultsDF$Rho_deplete < cutoff)
                                      ),
                       'direction' = vapply(resultsDF$`Median log2 Fold Change`, function(x){ifelse(x < 0, 'Deplete', 'Enrich')}, character(1)),
                       'guides' = as.numeric(guides[as.character(resultsDF$geneSymbol)]),
                       stringsAsFactors = FALSE)
  row.names(oraIpt) <- oraIpt$feature_id
  
  #Check if the Db should be GREAT
  if(!all(table(oraIpt$feature_id) == 1)){
    warning("Multiple targets are associated with some features. Consider converting to a GREATdb if you haven't done so already; see gCrisprTools::ct.GREATdb().")
  }
                       
fb <- ifelse(unbiased, NULL, 'guides')

  
msdb <- multiGSEA::multiGSEA(gsd = gsdb, x = oraIpt, methods = 'ora', groups = 'direction', rank_by = 'guides', feature.bias = 'guides')                     
                       
# When we run data through the "multiGSEA pipeline" as opposed to calling
# methods like "ora" directlly, the *eventual* data.frame that materializes
# that ora() will execute over expects to have certain things in place.
#
# In this case, it is expected that the features (rows) that are "enriched"
# is a TRUE/FALSE column named 'significant'. So we just have to make it so
# here.
in.df <- transform(in.df, significant = sig)
out <- multiGSEA::multiGSEA(gsdb, oraIpt, methods = c("ora", "cameraPR"),
                 groups = "direction", rank_by = "guides",
                 feature.bias = "guides")
  

  #identify significant targets
  enr <- 
  dep <- switch(stat,
                'q' = resultsDF$geneID[resultsDF$`Target-level Depletion Q` < cutoff], 
                'p' = resultsDF$geneID[resultsDF$`Target-level Depletion P` < cutoff], 
                'rho' = resultsDF$geneID[resultsDF$Rho_deplete < cutoff]
                )
  
  raw.hgt <-  hyperGeometricTest(gsd = genesetDB, selected = up, universe = universe, do.conform = TRUE) 
  
  
}
  
  

