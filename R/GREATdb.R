##' Update a gene-centric msdb object for GREAT-style enrichment analysis using a specified CRISPR annotation.
##'
##' @param annotation an annotation object returned by \code{ct.prepareAnnotation()}. 
##' @param gsdb A gene-centric \code{GeneSetDb} object to conform to the relevant peakwise dataset.
##' @param minsize Minimum number of targets required to consider a geneset valid for analysis.
##' @return A new \code{GeneSetDb} object with the features annotated genewise to pathways.
##' @importFrom multiGSEA getMSigGeneSetDb
##' @export
ct.GREATdb <- function(annotation, gsdb = getMSigGeneSetDb(c('h', 'c2'), 'human'), minsize = 10){
  
  #check the anbnotation object
  annotation <- ct.prepareAnnotation(annotation)
  
  #Check to make sure that the geneID column values are present in the provided gsdb
  genes <-na.omit(unique(annotation$geneID))
  message(paste0(length(intersect(featureIds(gsdb), genes)), ' valid genes found in the supplied gene database (of ', length(genes), ')'))

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
