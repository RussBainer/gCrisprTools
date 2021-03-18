##' @title Check and optionally subset an annotation file for use in a Crispr Screen
##' @description This function processes a supplied annotation object for use in a pooled screening experiment. 
##' Originally this was processed into something special, but now it essentially returns
##' the original annotation object in which the geneSymbol column has been factorized. This is primarily used 
##' internally during a call to the \code{ct.generateResults()} function. Also performs some minor functionality checking.
##' 
##' Valid annotations contain both `geneID` and `geneSymbol` columns. This is because there is often a distinction between
##' the official gene that is being targeted and a coherent set of gRNAs that make up a testing cohort. For example, 
##' multiple sets of guides may target distinct promoters, exons, or other entities that are expected to produce distinct 
##' biological phenomena related to the gene that should be interpreted separately. For this reason, the `geneID` column 
##' encodes the official gene designation (typically an ensembl or entrez gene identifier) while the `geneSymbol` column
##' contains a human-readable descriptor of the gRNA target (such as a gene symbol or promoter name). 
##' 
##' @param ann A \code{data.frame} containing an annotation object with gRNA-level information encoded as rows. The 
##' \code{row.names} attribute should correspond to the individual gRNAs, and it should at minimum contain columns 
##' named "geneID" and "geneSymbol" indicating the corresponding gRNA target gene ID and symbol, respectively. 
##' @param object If supplied, an object with \code{row.names} to be used to subset the supplied annotation frame
##' for downstream analysis.
##' @param controls The name of a value in the \code{geneSymbol} column of \code{ann} that corresponds to nontargeting 
##' control gRNAs. May also be supplied as a logical value, in which case the function will try to identify and format 
##' nontargeting guides. 
##' @param throw.error Logical indicating whether to throw an error when \code{controls} is \code{TRUE} but no nontargeting
##' gRNAs are detected. 
##' @return A new annotation data frame, usually with nontargeting controls and NA values reformatted to \code{NoTarget} 
##' (and \code{geneID} set to \code{'no_gid'}),  and the "geneSymbol" column of \code{ann} factorized. If supplied with 
##' an \code{object}, the gRNAs not present in the \code{object} will be omitted. 
##' @author Russell Bainer
##' @examples 
##' data('ann')
##' data('es')
##' es <- ct.filterReads(es)
##' newann <- ct.prepareAnnotation(ann, es)
##' @export

ct.prepareAnnotation <- function(ann, object = NULL, controls = TRUE, throw.error = TRUE){
  if(!(length(setdiff(c("geneSymbol", "geneID"), colnames(ann))) == 0)){
    stop("The annotation object must contain a geneSymbol column.")
  }
  
  exception <- ifelse(throw.error, stop, warning)

  #Convert the geneID column to a charvec if it isn't already.
  if(!is.character(ann$geneID)){
    ann$geneID <- as.character(ann$geneID)
  }
  if(any(is.na(ann$geneID))){
    warning('NA detected within the geneID column; setting to no_gid...')
    ann$geneID[is.na(ann$geneID)] <- "no_gid"
  }
  
  
  #If supplied, trim out the gRNAs not present in the object. 
  if(!is.null(object)){
    
    if(length(row.names(object)) > 0){

      if(any(!is.element(row.names(object), row.names(ann)))){
        stop("The supplied object contains rows not present in the annotation file. Please omit these elements prior to downstream analyses.")
        }
    
      omit <- length(setdiff(row.names(ann), row.names(object)))
      if(omit > 0){
        message(paste(omit, "elements defined in the annotation file are not present in row names of the specified object. Omitting."))
        }

      ann <- ann[row.names(object),]

      } else {
          warning('The supplied object has no row.names! Ignorning.')
        }
    }
  
  if(!is.character(controls)){
    geneSymb <- "NoTarget" 
    } else {
      geneSymb <- as.character(controls) 
      controls <- TRUE
      }
  
    if(controls){  
      ntc <- NULL
      if(any(geneSymb %in% ann$geneSymbol)){
        ntc <- row.names(ann)[ann$geneSymbol %in% geneSymb]      
        } else {
          message(paste(geneSymb, "is not present in the geneSymbol column of the specified annotation file; trying to find something nontargeting..."))
          
          if((sum(is.na(ann$geneSymbol)) > 0) && (sum(ann$geneID %in% "no_gid") > 0)){
              message('NA and "no_gid" elements are both present in the supplied annotation, so I am using the "no_gid" elements. If you wish to select another set of gRNAs, please change geneSymb.')  
            } 
        
          if(sum(ann$geneID %in% "no_gid") > 0){
            ntc <- row.names(ann)[ann$geneID %in% "no_gid"]
            message(paste('I found some gRNAs targeting "no_gid".', "Let's use that." ))
          } else if (sum(is.na(ann$geneSymbol)) > 0){
              ntc <- row.names(ann)[is.na(ann$geneSymbol)]
              message('No "no_gid" geneIDs in the annotation file, Using gRNAs targeting geneSymbol NA.')
            } else {
              exception('I cannot find any obvious control guides in the supplied annotation, so I am assuming that they are absent. Please specify a geneSymbol if you know something I do not.')
              NULL
              }
          }
      #Compel the NTCs to the proper values
      if (!is.null(ntc)) {
        ann[ntc, "geneSymbol"] <- "NoTarget"
        ann[ntc, "geneID"] <- "no_gid"
      }
      
      }
      
    #Plug any remaining holes.
    if(sum(is.na(ann[,"geneSymbol"])) > 0){
      warning('NAs present in the annotation geneSymbol column. These will be collected into an element called "Unknown".', call. = TRUE)
      ann[is.na(ann[,"geneSymbol"]),"geneSymbol"] <- "Unknown"
    }  
  ann$geneSymbol <- as.factor(ann$geneSymbol)
  ann <- droplevels(ann)
  return(ann)
}


##' @title Expand an annotation object to allow annotations of reagents to multiple targets
##' @description This function takes a gCrisprTools annotation object and expands it to allow 1:many mappings of reagents. 
##' This mostly is used for internal processing, and users should interact with the wrapper functions that call it (e.g., `ct.generateResults`). 
##' 
##' Libraries targeting ambiguous biological elements (e.g., alternative promoters to a gene where the boundaries between
##' elements is contested) may contain reagents that are plausibly annotated to a finite set of possible targets. To accomodate this, users may 
##' supply an alternative reagent annotation in the form of a named list of vectors, where the names correspond to reagent `ID`s in the annotation 
##' object and each list element corresponds something coercible to a to a character vector of associated targets that will ultimately be assembled
##' into the `geneSymbol` column of the annotation object. It is assumed that the `geneID` values are assigned unambiguously to the reagents, 
##' and are passed through directly. 
##'
##' @param ann A \code{data.frame} containing an annotation object with gRNA-level information encoded as rows, typically produced by 
##' `ct.prepareAnnotation`. The `ID` column should correspond to the individual reagent identifiers. 
##' @param alt.annotation A named list of character vectors, which should be named identically to a value in the `ID` column of the supplied 
##' annotation object. The values in the character vectors will eventually form the `geneSymbol` column of the annotation file.
##' @return A new annotation data frame, expanded as described above.
##' @author Russell Bainer
##' @examples 
##' data('ann')
##' alt.annotation <- list("Target2089_b" = c('foo', 'bar'), "Target2089_c" = 'foo', "Target2089_a" = 'bar')
##' ct.expandAnnotation(ann, alt.annotation)
##' @export
ct.expandAnnotation <- function(ann, alt.annotation){
  #input check
  ann <- ct.prepareAnnotation(ann)
  if(!('ID' %in% names(ann))){
    stop('Annotation must contain a reagent-level "ID" column to enable expansion.')
  }
  stopifnot(is.list(alt.annotation), all(names(alt.annotation) %in% ann$ID))
  alt.annotation <- sapply(alt.annotation, as.character, simplify = FALSE)
  
  assigned <- vapply(alt.annotation, length, numeric(1)) > 0
  alt.annotation <- alt.annotation[assigned]
  
  if(length(setdiff(ann$ID, names(alt.annotation))) > 0){
    message(paste0('Removing ', length(setdiff(ann$ID, names(alt.annotation))), ' reagents that were not assigned to a target.'))
  }
  
  subs <- lapply(1:length(alt.annotation), 
                 function(x){
                   out <- ann[ann$ID %in% names(alt.annotation)[x], , drop = FALSE]
                   out <- out[rep(1, times = length(alt.annotation[[x]])),, drop = FALSE]
                   out$geneSymbol <- alt.annotation[[x]]
                   return(out)
                 })
  
  return(do.call(rbind, subs))
}


##' @title Parse `geneSymbol` values to construct an alternative annotation list
##' @description This is an accessory function to \link{\code{ct.expandAnnotation()}} function, which enables users to expand annotation objects 
##' to accomodate reagent libraries where reagents are expected to impact a set of known targets. See documentation for that function for 
##' additional details. 
##' 
##' Often, libraries that contain multiply-targeting reagents are annotated using a structured format that can be decomposed by regex matching. 
##' This function takes in an annotation object containing an `ID` column indicating the reagent ID and a `geneSymbol` column containing the 
##' target mappings, and parses the target mappings according to a known annotation format. Currently supported formats are "cellecta" 
##' (e.g., "TARGET_P1P2P3" indicating multiple promoters associated with a known target), and "underscore", where different targets are 
##' concatenated using the underscore separator (e.g., "TARGET1_TARGET2_TARGET3").
##' 
##' Returns an `alt.annotation`-type list of character vectors encoding the target mappings for each reagent. 
##'
##' @param ann A \code{data.frame} containing reagent-level information encoded as rows. The `ID` column should correspond to the individual 
##' reagent identifiers, and the "geneSymbol" column should contain target annotation strings to be parsed (both are coerced to strings). 
##' Does not, strictly speaking, need to be a proper annotation object, but one of those will work. 
##' @param format Format of the geneSymbol column strings.
##' @return A named `alt.annotation`-type list of character vectors encoding the target mappings for each reagent
##' @author Russell Bainer
##' @examples 
##' fakeann <- data.frame('ID' = LETTERS[1:4], 'geneSymbol' = c('T1_P1', 'T1_P1P2', 'T1_P2P1', 'T1_P2'))
##' ct.parseGeneSymbol(fakeann, 'cellecta')
##' ct.parseGeneSymbol(fakeann, 'underscore')
##' @export
ct.parseGeneSymbol <- function(ann, format = c('cellecta', 'underscore')){
  
  format <- match.arg(format)
  stopifnot(all(c('ID', 'geneSymbol') %in% names(ann)))
  
  reagents <- as.character(ann$ID)
  symbols <- as.character(ann$geneSymbol)
  
  if(format %in% 'cellecta'){
    out <- lapply(symbols, 
                  function(x){
                    splut <- strsplit(x, split = '_')[[1]]
                    prom <- strsplit(splut[2], split = 'P')[[1]]
                    return(paste0(splut[1], '_P', prom[2:length(prom)]))
                  })
  }
  if(format %in% 'underscore'){
    out <- lapply(symbols, 
                  function(x){
                    return(strsplit(x, split = '_')[[1]])
                  })
  }
  
  names(out) <- reagents
  return(out)
}



