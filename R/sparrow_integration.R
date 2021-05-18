##' Update a gene-centric msdb object for GREAT-style enrichment analysis using a specified CRISPR annotation.
##'
##' Update a gene-centric `GeneSetDb` object for GREAT-style enrichment analysis using a specified annotation.
##' 
##' Often, pooled screening libraries are constructed such that the gene targets of interest are associated 
##' with variable numbers of semi-independent screen signals (associated with, e.g., sets of alternative promoters 
##' or cis regulatory units). Such an arrangement is often unavoidable but produces to complications when 
##' performing gene set enrichent analyses. This function conforms a standard `GeneSetDb` object to appropriately 
##' consider this form of ultiple testing during ontological enrichment analyses according to the GREAT strategy
##' outlined by [McLean et al. (2009)](https://doi.org/10.1038/nbt.1630).
##' 
##' Operationally, this means that genewise sets in the provided object will be translated to the corresponding 
##' `geneSymbol` sets provided in the annotation file. 
##' 
##' @param annotation an annotation object returned by \code{ct.prepareAnnotation()}. 
##' @param gsdb A gene-centric \code{GeneSetDb} object to conform to the relevant peakwise dataset.
##' @param minsize Minimum number of targets required to consider a geneset valid for analysis.
##' @param ... Additional arguments to be passed to `ct.prepareAnnotation()`.
##' @return A new \code{GeneSetDb} object with the features annotated genewise to pathways.
##' @examples 
##' data(resultsDF)
##' data(ann)
##' gsdb <- ct.GREATdb(ann, gsdb = getMSigGeneSetDb(collection = 'h', species = 'human', id.type = 'entrez'))
##' show(featureIds(gsdb))
##' @export
ct.GREATdb <- function(annotation, 
                       gsdb = getMSigGeneSetDb(collection = c('h', 'c2'), 
                                               species = 'human', 
                                               id.type = 'ensembl'), 
                       minsize = 10, 
                       ...){
  
  #input check
  if (!is(gsdb, "GeneSetDb")) {
    if (is(gsdb, "GeneSetCollection") || is(gsdb, "GeneSet")) {
      stop("A GeneSetDb is required. GeneSetCollections can be converted to a ", 
           "GeneSetDb via the GeneSetDb constructor, or a call to ", 
           "`as(gene_set_collection, 'GeneSetDb')`. See ?GeneSetDb for more ", 
           "details.")
    }
    stop("GeneSetDb required. Please see `?GeneSetDb` for ways to turn your ", 
         "gene sets into a GeneSetDb, or ")
  }
  stopifnot(is(gsdb, 'GeneSetDb'), is(minsize, 'numeric'), (length(minsize) == 1))
  
  #check the annotation object
  annotation <- ct.prepareAnnotation(annotation, ...)
  
  #Check to make sure that the geneID column values are present in the provided gsdb
  genes <-na.omit(unique(annotation$geneID))
  message(paste0(length(intersect(sparrow::featureIds(gsdb), genes)), ' valid genes found in the supplied gene database (of ', length(genes), ')'))

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
  return(sparrow::GeneSetDb(out, ...))
}

##' Prepare one or more resultsDF objects for analysis via Sparrow. 
##' 
##' Take in a list of results objects and return an equivalently-named list of input `data.frames` appropriate for `sparrow::seas()`. 
##' By construction, the relevant target unit is extracted from the `geneSymbol` column of the provided results objects, which may. Note that the 
##' genewise `@logFC` slot in the returned object will contain the appropriately-signed Z transformation of the P-value 
##' assigned to the target. In most applications this is arguably more interpretable than e.g., the median gRNA log2 fold change.  
##'  
##' @param dflist A list of gCrisprTools results `data.frames` to be formatted.
##' @param collapse.on Should targets be annotated as `geneSymbol`s or `geneID`s (default)? 
##' @param cutoff Numeric maximum value of `statistic` to define significance.
##' @param statistic Should cutoffs be calculated based on FDR (`best.q`) or P-value (`best.p`)?
##' @param gdb Optionally, a `GeneSetDb` object to enable proper registration of the output. If provided, the 
##' collapsing features in the provided `simpleDF`s must be present in the `gsd@db$feature_id` slot. Note that a GREAT-style `GeneSetDb` that 
##' has been conformed via `ct.GREATdb()` will use `geneID`s as the `feature_id`.
##' @return A list of `data.frames` formatted for evaluation with `sparrow::seas()`. 
##' @examples 
##' data(resultsDF)
##' ct.seasPrep(list('longer' = resultsDF, 'shorter' = resultsDF[1:10000,]), collapse.on = 'geneSymbol')
##' @export
ct.seasPrep <- function(dflist, 
                        collapse.on = c('geneID', 'geneSymbol'),
                        cutoff = 0.1, 
                        statistic = c('best.q', 'best.p'), 
                        regularize = FALSE,
                        gdb = NULL){
  #Input check
  collapse.on <- match.arg(collapse.on)
  stopifnot(is(cutoff, 'numeric'), cutoff <= 1, cutoff >= 0, is(regularize, 'logical'))
  
  if(regularize){
    dflist <- ct.regularizeContrasts(dflist, collapse = collapse.on)
  } else {
    dflist <- sapply(dflist, ct.simpleResult, collapse = collapse.on, simplify = FALSE)
  }

  statistic <- match.arg(statistic)
  
  if(!is.null(gdb)){
    dflist <- sapply(dflist, 
                     function(x){
                       x[(row.names(x) %in% gdb@db$feature_id),]
                     }, simplify = FALSE)
    message('Removed genes absent from the provided GeneSetDb.')
  }

  if(any(unlist(lapply(dflist, nrow)) == 0)){
    stop('No valid `feature_id`s observed for one or more of the provided contrasts! check `collapse.on`?')
  }

  out <- sapply(dflist, 
                function(x){
                  z <- 10^-(gCrisprTools:::ct.softLog(x$best.p))
                  z <- qnorm(z/2) * ifelse(x$direction == 'enrich', -1, 1)
                  
                  df <- data.frame('feature_id' = row.names(x),
                                   'logFC' =  z,
                                   'significant' = (x[,statistic] <= cutoff),
                                   'direction' = x$direction,
                                   'rank_by' = z,
                                   stringsAsFactors = FALSE)
                  return(df)
                },simplify = FALSE)

  names(out) <- names(dflist)
  return(out)
  
}
                         
##' @title Geneset Enrichment within a CRISPR screen using `sparrow`
##' 
##' This function is a wrapper for the `sparrow::seas()` function, which identifies differentially enriched/depleted ontological 
##' categories within the hits identified by a pooled screening experiment, given a provided `GenseSetDb()` object and a list of 
##' results objects created by `ct.generateResults()`. By default testing is performed using `fgsea` and a hypergeometric test 
##' (`sparrow::ora()`), and results are returned as a `SparrowResult` object. 
##' 
##' This function will attempt to coerce them into inputs appropriate for the above analyses via `ct.seasPrep()`, after checking 
##' the relevant parameters within the provided `GeneSetDb`. This is generally easier than going through the individual steps yourself, 
##' especially when 
##' 
##' Note that many pooled libraries specifically target biased sets of genes, often focusing on genes involved
##' in a particular pathway or encoding proteins with a shared biological property. Consequently, the enrichment results
##' returned by this function represent the disproportionate enrichment or depletion of targets annotated to pathways *within the context
##' of the screen*, and may or may not be informative of the underlying biology in question. This means that
##' pathways not targeted by a library will obviously never be enriched a positive target set regardless of
##' their biological relevance, and pathways enriched within a focused library screen are similarly expected to partially
##' reflect the composition of the library and other confounding issues (e.g., number of targets within a pathway).
##' Analysts should therefore use this function with care. For example, it might be unsurprising to detect pathways related
##' to histone modification within a screen employing a crispr library primarily targeting epigenetic regulators.
##'  
##' @param dflist A result object created by `ct.generateResults()`, or a named list containing many of them; will be passed as a list to 
##' `ct.seasPrep()` with the associated `...` arguments.
##' @param gdb A `GenseSetDb` object containing annotations for the targets specified in `result`.
##' @param as.dfs Logical indicating whether to return the various contrast statistics as in-register lists of data.frames to facilitate comparisons.
##' @param ... Additional arguments to pass to `ct.seasPrep()` or `sparrow::seas()`. 
##' @return A named list of `SparrowResults` objects.
##' @examples 
##' data('resultsDF')
##' ct.seas(list('longer' = resultsDF, 'shorter' = resultsDF[1:10000,]), gdb = getMSigGeneSetDb(collection = 'h', species = 'human', id.type = 'entrez'))
##' @author Steve Lianoglou for seas; Russell Bainer for GeneSetDb processing and wrapping functions.
##' @export
ct.seas <- function(dflist,
                    gdb, 
                    as.dfs = FALSE,
                    ...){
  library('sparrow', quietly = TRUE)
  #Check GSDB and determine feature set
  stopifnot(is(gdb, 'GeneSetDb'), is(as.dfs, 'logical'))
  
  #listify results as needed
  if(!is(dflist, 'list')){
    if(ct.resultCheck(dflist)){
      dflist <- list('result' = dflist)
    }
  }

  
  #Infer whether Gsdb is ID or feature centric
  gids <- sum(gdb@featureIdMap$feature_id %in% dflist[[1]]$geneID)
  gsids <- sum(gdb@featureIdMap$feature_id %in% dflist[[1]]$geneSymbol)
  
  if(all(c(gsids, gids) == 0)){
    stop('None of the features in the GeneSetDb are present in either the geneID or geneSymbol slots of the first provided result.')
  }
  
  identifier <- ifelse(gids > gsids, 'geneID', 'geneSymbol')
  message(paste0('GeneSetDb feature_ids coded as ', identifier, 's.'))
  if(identifier %in% 'geneID'){
    message('Depending on the composition of your library, you might consider switching to a target-level analysis; see ?ct.GREATdb() for details.')
    }
  
  #Check that all provided objects are keyed to the proper values
  ipts <- ct.seasPrep(dflist, 
                      collapse.on = identifier, 
                      gdb = gdb)
  
  outs <- sapply(ipts, 
                function(x){
                  seas(gdb, x, methods = c('ora', 'fgsea'), rank_by = 'rank_by', selected = 'selected', groups = 'direction', ...)
                }, 
                simplify = FALSE)
  
  if(as.dfs){
    outs <- ct.compileSparrow(outs)
  }
  return(outs)
}

##' Compile the values from a set of SparrowResult Objects
##' 
##' This function takes in a named list of `SparrowResult` objects and breaks them into in0register data.frames for easy comprisons. 
##' Specifically, the function assembles the values in each of the `@result` slots for each of the provided contrasts into standalone 
##' `data.frame`s with the rows named for the pathways, and returns these objects as a list-of-lists for each result type.
##' 
##' Note that results are compiled as-is, so users need to orient the screen results themselves in whatever manner they deem appropriate. 
##' 
##' @param resultList A named list of `SparrowResult` objects
##' @value A (test-level) list of (output-level) lists of (statistic-level) dataframes, such that category statistics can be easily 
##' compared to one another.   
##' @author Russell Bainer
##' @export
ct.compileSparrow <- function(resultList){
  stopifnot(is(resultList, 'list'), all(sapply(resultList, is, "SparrowResult")), !is.null(names(resultList)))
  
  sapply(names(outs[[1]]@results),                
         function(tests){
           sapply(names(outs[[1]]@results[[tests]])[2:length(names(outs[[1]]@results[[tests]]))],       
                  function(testcols){
                    ret <- sapply(names(outs),
                                  function(sparrowres){
                                    outs[[sparrowres]]@results[[tests]][[testcols]]
                                  }, simplify = TRUE)
                    row.names(ret) <- outs[[sparrowres]]@results[[tests]][[1]]
                    return(ret)
                  }, simplify = FALSE)
         }, simplify = FALSE)

}



