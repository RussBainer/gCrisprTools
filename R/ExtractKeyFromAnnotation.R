##' @title Check and optionally subset an annotation file for use in a Crispr Screen
##' @description This function processes a supplied annotation object for use in a pooled screening experiment. 
##' Originally this was processed into something special, but now it essentially returns
##' the original annotation object in which the geneSymbol column has been factorized. This is primarily used 
##' internally during a call to the \code{ct.generateResults()} function. Also performs some minor functionality checking.
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
          
          if((sum(is.na(ann$geneSymbol)) > 0) && (sum(ann$geneID == "no_gid") > 0)){
              message('NA and "no_gid" elements are both present in the supplied annotation, so I am using the "no_gid" elements. If you wish to select another set of gRNAs, please change geneSymb.')  
            } 
        
          if(sum(ann$geneID == "no_gid") > 0){
            ntc <- row.names(ann)[ann$geneID == "no_gid"]
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


