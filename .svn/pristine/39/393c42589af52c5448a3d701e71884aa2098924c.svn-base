
# helper functions --------------------------------------------------------

# These are internal functions to re-use for report generation so I am am not
# sticking to the naming convention with a 'ct' prefix, which is used for API.


#' Add formatted timestamp and extension to a file name
#'
#' @param name character
#' @param ext character - extension including a ".".
#'
#' @keywords internal
#' @return character
#'
appendDateAndExt <- function(name, ext) {
  paste0(name,
         format(Sys.time(), "_%b_%d_%Y_%H.%M.%S"),
         ext)
}

#' Internal wrapper to generate html markdown reports from existing templates
#'
#' @param reportNameBase character - name of the report's file.
#' @param templateName character - name of the rmarkdown template file in the
#'   standard location. Passed through to \code{template} argument of the
#'   \code{\link[rmarkdown]{draft}}.
#' @param rmdParamList list of named report parameters. Passed through to
#'   \code{params} argument of the \code{\link[rmarkdown]{render}}.
#' @param outdir An optional character string indicating the directory in which to 
#' generate the report. If \code{NULL}, a temporary directory will be automatically generated. 
#'
#' @keywords internal
#' @return  character with a path to html report in the temporary directory.
renderReport <- function(reportNameBase,
                         templateName,
                         rmdParamList, 
                         outdir = NULL) {
  
  if(is.null(outdir)){
      # use a temp directory to make a draft
      dir.create(outdir <- tempfile())
    } else {
        if (!is.character(outdir)) {
          stop("`outdir` must be path to a directory used ot save results")
        }
        outdir.created <- initOutDir(outdir)   #True/false indicating whether the directory was created
        outdir <- normalizePath(outdir)
      }

  render(
    draft(
      file = file.path(outdir, reportNameBase),
      create_dir = TRUE,
      template = templateName,
      package = "gCrisprTools",
      edit = FALSE
    ),
    params = rmdParamList
  )
}


# ## ----------------------------------------------------------------------


##' @title Generate a full experimental report from a pooled CRISPR screen
##' @description This is a function to generate an html report for a CRISPR screen, incorporating information about a specified contrast. 
##' The report contains a combination of experiment-level and contrast-specific analyses, largely collected from other functions in 
##' \code{gCrisprTools}. It is designed to be used 'as-is', and analysts interested in using different functionalities of the various 
##' functions should do that outside of this wrapper script.
##' @param fit An object of class \code{MArrayLM} containing, at minimum, a \code{coefficents} slot with coefficients from the comparison,
##' and a \code{stdev.unscaled} slot with the corresponding standard deviation of the coefficent estimates. The \code{row.names} attribute
##' should ideally match that which is found in \code{annotation}, but this will be checked internally.
##' @param eset An ExpressionSet object containing, at minimum, a matrix of gRNA abundances extractable with the \code{exprs()} function and some named
##' phenodata extractable with \code{pData()}.
##' @param sampleKey A sample key, supplied as an ordered factor linking the samples to experimental
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control set is assumed to be
##' the first \code{level}.
##' @param annotation An annotation object for the experiment. See the man page for \code{ct.prepareAnnotation()} for details and example format. 
##' @param results A data.frame summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}.
##' @param aln A numeric alignment matrix, where rows correspond to "targets", "nomatch", "rejections", and "double_match", and where columns 
##' correspond to experimentasl samples.
##' @param outdir A directory in which to generate the report; if \code{NULL}, a temporary directory will be automatically generated. 
##' The report will be located in a subdirectory whose name is internally generated (see below). The path to the report itself is returned by the function.
##' @param contrast.term A parameter passed to \code{ct.preprocessFit} in the event that the fit object contains data from multiple contrasts. See
##' that man page for further details.
##' @param identifier A character string to name the report and corresponding subdirectories. If provided, the final report will be called
##' '\code{identifier}.html' and will be located in a directory called \code{identifier} in the \code{outdir}. If \code{NULL}, a generic name 
##' including the timestamp will be generated.
##' @return The path to the generated html report.
##' @author Russell Bainer
##' @examples 
##' data('fit')
##' data('es')
##' 
##' ##' #Build the sample key
##' library(Biobase)
##' sk <- relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference")
##' names(sk) <- row.names(pData(es))
##' 
##' data('ann')
##' data('resultsDF')
##' data('aln')
##' path2report <- ct.makeReport(fit, es, sk, ann, resultsDF, aln, outdir = ".") 
##' @export

ct.makeReport <- function(fit, eset, sampleKey, annotation, results, aln, outdir = NULL, contrast.term = NULL, identifier = NULL){

  if(class(fit) != "MArrayLM"){
    stop("the provided fit does not appear to be a MArrayLM object.")
    }

  if(ncol(fit$coefficients) > 1){
    if(is.null(contrast.term)){
      stop("The fit object contains multiple coefficients. Please specify a contrast.term.")
    }
    fit <- ct.preprocessFit(fit, contrast.term)
  }

  #filter the annotation file as necessary for downstream processes
  annotation <- ct.prepareAnnotation(annotation, fit, throw.error = FALSE)

  if(!is.matrix(aln) | !setequal(row.names(aln), c("targets", "nomatch", "rejections", "double_match"))){
    stop("I don't think that the provided alignment matrix is actually an alignment matrix.")
  }
  if(!ct.resultCheck(results)){
    stop("Execution halted.")
  }
  
  
#Set up the path and folder for writing (heavily based on code from multiGSEA):
if (is.null(outdir)) {
  outdir.created <- FALSE
  noutdir <- character()
  } else {
  if (!is.character(outdir)) {
    stop("`outdir` must be path to a directory used ot save results")
  }
  outdir.created <- initOutDir(outdir)   #True/false indicating whether the directory was created
  noutdir <- normalizePath(outdir)
}

#make the sample names and build the parameter list for the rmd.
if(ct.inputCheck(sampleKey, eset)){
  sampleKey <- sampleKey[order(sampleKey)]
  }

rmdParamList <- list(fit= fit,
                    eset= eset,
                    sampleKey= sampleKey,
                    results = results,
                    annotation = annotation,
                    aln = aln)

#make the name of the rmd & output, and render it in the location of interest.
if(!is.null(identifier)){
  if(!is.character(identifier) | length(identifier) != 1){
    stop("identifier must be a character string of length 1.")
    }
  namebase <- identifier
  }else{
    namebase <- paste0("CrisprContrastReport", format(Sys.time(), "_%b_%d_%Y_%H.%M.%S"))
    }

rmdname <- file.path(noutdir, paste0(namebase, '.Rmd'))
outname <- file.path(noutdir, namebase, paste0(namebase, '.html'))

print(names(rmdParamList))

render(draft(rmdname, template = 'CRISPR_report', package = "gCrisprTools", edit = FALSE), params = rmdParamList)

return(outname)

}

##' @title Generate a QC report from a pooled CRISPR screen
##' @description This is a function to generate an html QC report for a CRISPR 
##'   screen, focusing on experiment-level and library-level analyses collected 
##'   from other functions in \code{gCrisprTools}. It is designed to be used 
##'   'as-is', and analysts interested in using different functionalities of the
##'   various functions should do that outside of this wrapper script.
##'   
##' @param eset An ExpressionSet object containing, at minimum, a matrix of gRNA
##'   abundances extractable with the \code{exprs()} function and some named 
##'   phenodata extractable with \code{pData()}.
##' @param trim The number of gRNAs to be trimmed from the top of the 
##'   distribution before estimating the abundance range. Empirically, this 
##'   usually should be equal to about 2 to 5 percent of the guides in the 
##'   library.
##' @param log2.ratio Maximum abundance of contaminant gRNAs, expressed on the
##'   log2 scale from the top of the trimmed range of each sample. That is,
##'   \code{log2.ratio = 4} means to discard all gRNAs whose abundance is
##'   (1/2)^4 of the trimmed maximum.
##' @param sampleKey A sample key, supplied as an ordered factor linking the 
##'   samples to experimental variables. The \code{names} attribute should 
##'   exactly match those present in \code{eset}, and the control set is assumed
##'   to be the first \code{level}.
##' @param annotation An annotation object for the experiment. See the man page 
##'   for \code{ct.prepareAnnotation} for details and example format.
##' @param aln A numeric alignment matrix, where rows correspond to "targets", 
##'   "nomatch", "rejections", and "double_match", and where columns correspond 
##'   to experimentasl samples.
##' @param identifier A character string to name the report and corresponding 
##'   subdirectories. If provided, the final report will be called 
##'   '\code{identifier}.html' and will be located in a directory called 
##'   \code{identifier}. If \code{NULL}, a generic name including the timestamp 
##'   will be generated.
##' @param lib.size Optional numeric vector of limma-type size factors
##' @param geneSymb The \code{geneSymbol} identifier(s) in \code{annotation} 
##'   that corresponds to gRNAs to be plotted on the curves. Passed through to 
##'   \code{\link{ct.gRNARankByReplicate}}, \code{\link{ct.viewControls}} and 
##'   \code{\link{ct.prepareAnnotation}} (as \code{controls} argument if it's
##'   not \code{NULL}). Default \code{NULL}.
##' @param outdir An optional character string indicating the directory in which
##'   to generate the report. If \code{NULL}, a temporary directory will be
##'   automatically generated.
##' @return The path to the generated html report.
##' @author Russell Bainer, Dariusz Ratman
##' @examples 
##' data('es')
##' data('ann')
##' data('aln')
##' 
##' ##' #Build the sample key
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' 
##' path2report <- ct.makeQCReport(es, trim = 1000, log2.ratio = 0.0625, sk, ann, aln, identifier = NULL, lib.size = NULL, outdir = ".") 
##' @export
ct.makeQCReport <-
  function(eset,
           trim,
           log2.ratio,
           sampleKey,
           annotation,
           aln,
           identifier = NULL,
           lib.size,
           geneSymb = NULL, 
           outdir = NULL) {
    # check and preprocess inputs
    if (!is.null(sampleKey)) {
      ct.inputCheck(sampleKey, eset)
      sampleKey <- sampleKey[order(sampleKey)]
    }
    
    annotation <- ct.prepareAnnotation(
      ann = annotation,
      object = eset,
      controls = ifelse(is.null(geneSymb), TRUE, geneSymb),
      throw.error = FALSE
    )
    
    if (!is.matrix(aln) |
        !setequal(row.names(aln),
                  c("targets", "nomatch", "rejections", "double_match"))) {
      stop("I don't think that the provided alignment matrix is actually an alignment matrix.")
    }
    

    # prepare params
    if (missing(identifier) | is.null(identifier)) {
      reportNameBase <- appendDateAndExt(name = "CRISPR_QC_report",
                                         ext = ".Rmd")
    } else {
      (is.character(identifier) & length(identifier) == 1) ||
        stop("Identifier must be a character string of length 1.")
      reportNameBase <- identifier
    }

    rmdParamList <- list(
      eset = eset,
      sampleKey = sampleKey,
      trim = trim,
      log2.ratio = log2.ratio,
      annotation = annotation,
      aln = aln,
      lib.size = lib.size,
      geneSymb = geneSymb
    )

    # render and return the path to report
    renderReport(
      reportNameBase = reportNameBase,
      templateName = 'CRISPR_QC_report',
      rmdParamList = rmdParamList, 
      outdir = outdir
    )
  }

##' @title Generate a Contrast report from a pooled CRISPR screen
##' @description This is a function to generate an html Contrast report for a 
##'   CRISPR screen, focusing on contrast-level analyses collected from other 
##'   functions in \code{gCrisprTools}. It is designed to be used 'as-is', and 
##'   analysts interested in using different functionalities of the various 
##'   functions should do that outside of this wrapper script.
##'   
##' @param eset An ExpressionSet object containing, at minimum, a matrix of gRNA
##'   abundances extractable with the \code{exprs()} function and some named 
##'   phenodata extractable with \code{pData()}.
##' @param fit A fit object for the contrast of interest, usually generated with
##'   \code{lmFit}.
##' @param sampleKey A sample key, supplied as an ordered factor linking the 
##'   samples to experimental variables. The \code{names} attribute should 
##'   exactly match those present in \code{eset}, and the control set is assumed
##'   to be the first \code{level}.
##' @param results A data.frame summarizing the results of the screen, returned 
##'   by the function \code{\link{ct.generateResults}}.
##' @param annotation An annotation object for the experiment. See the man page 
##'   for \code{ct.prepareAnnotation()} for details and example format.
##' @param comparison.id character with a name of the comparison.
##' @param identifier A character string to name the report and corresponding 
##'   subdirectories. If provided, the final report will be called 
##'   '\code{identifier}.html' and will be located in a directory called 
##'   \code{identifier} in the \code{outdir}. If \code{NULL}, a generic name
##' @param contrast.subset character vector containing the sample
##'   labels to be used in the analysis; all elements must be contained in the
##'   \code{colnames} of the specified \code{eset}. including the timestamp will
##'   be generated. Default: \code{colnames(eset)}.
##' @param outdir An optional character string indicating the directory in which to 
##'   generate the report. If \code{NULL}, a temporary directory will be automatically 
##'   generated. 
##' @return The path to the generated html report.
##' @author Russell Bainer, Dariusz Ratman
##' @examples 
##' data('es')
##' data('fit')
##' data('ann')
##' data('resultsDF')
##' 
##' ##' #Build the sample key
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' 
##' path2report <- ct.makeContrastReport(es, fit, sk, resultsDF, ann, comparison.id = NULL, outdir = ".") 
##' @export
ct.makeContrastReport <-
  function(eset,
           fit,
           sampleKey,
           results,
           annotation,
           comparison.id,
           identifier,
           contrast.subset = colnames(eset), 
           outdir = NULL
           ) {
    # check and preprocess inputs
    annotation <- ct.prepareAnnotation(annotation, throw.error = FALSE)
    if (!is.null(sampleKey)) {
      ct.inputCheck(sampleKey, eset)
      sampleKey <- sampleKey[order(sampleKey)]
    }
    
    if(!ct.resultCheck(results)){
      stop("Execution halted.")
    }
    
    # prepare params
    if (missing(identifier)) {
      reportNameBase <- appendDateAndExt(name = "CRISPR_contrast_report",
                                         ext = ".Rmd")
    } else {
      (is.character(identifier) & length(identifier) == 1) ||
        stop("Identifier must be a character string of length 1.")
      reportNameBase <- identifier
    }
    
    rmdParamList <- list(
      eset = eset,
      fit = fit,
      sampleKey = sampleKey,
      results = results,
      contrast.subset = contrast.subset,
      comparison.id = comparison.id,
      annotation = annotation
    )
    
    # render and return the path to report
    renderReport(
      reportNameBase = reportNameBase,
      templateName = 'CRISPR_contrast_report',
      rmdParamList = rmdParamList, 
      outdir = outdir
    )
  }

