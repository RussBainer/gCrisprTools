##' @title Normalize sample abundance estimates by the slope of the values in the central range
##' @description This function normalizes Crispr gRNA abundance estimates by equalizing the slopes of the middle (logged) values of the
##' distribution across samples. Specifically, the algorithm ranks the gRNA abundance estimates within each sample and determines a relationship between
##' rank change and gRNA within a trimmed region of the distribution via a linear fit. It then adjusts each sample such that the center of the logged
##' abundance distribution is strictly horizontal and returns these values as median-scaled counts in the appropriate slot of the input ExpressionObject.
##' @param ExpressionObject An ExpressionSet containing, at minimum, count data accessible by \code{exprs}, or an EList object with count data in the $E 
##' slot (usually returned by \link[limma]{voom}).
##' @param trim The proportion to be trimmed from each end of the distributionbefore performing the linear fit; algorithm defaults to 25% such that the
##' fit is performed on the interquartile range.
##' @param ... Other arguments to be passed to \code{ct.normalizeMedians()}, if desired.
##' @param lib.size An optional vector of size factor adjusted library size. 
##'   Default: \code{NULL} means to use sum of column counts as a lib.size.
##' @return A renormalized object of the same type as the provided object.
##' @author Russell Bainer
##' @import limma
##' @examples data('es')
##' data('ann')
##' 
##' #Build the sample key and library sizes for visualization
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' ls <- colSums(exprs(es))
##' 
##' es.norm <- ct.normalizeBySlope(es, lib.size= ls)
##' ct.gRNARankByReplicate(es, sk, lib.size= ls)
##' ct.gRNARankByReplicate(es.norm, sk, lib.size= ls)
##' @export
ct.normalizeBySlope <-
  function(ExpressionObject,
           trim = 0.25,
           lib.size = NULL,
           ...) {
    
  if(!(class(ExpressionObject) %in% c("ExpressionSet", "EList"))){stop(paste(deparse(substitute(ExpressionObject)), "is not an ExpressionSet or Elist."))}

  #log them
  if(class(ExpressionObject) == "EList"){
    e.dat <- ExpressionObject$E
  }else{
    counts <- exprs(ExpressionObject)
    if (is.null(lib.size)){
      lib.size <- colSums(counts)
     }else if(!is.numeric(lib.size) | length(lib.size) != ncol(counts)){
        stop('If specified, lib.size must be a numeric vector of the same length as the number of samples in the eset.')
      }  
    e.dat <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
    }

  #extract the inner 50% and fit a lm
  shuffle <- apply(e.dat, 2, sort, decreasing = TRUE, index.return = TRUE)
  locs <- round((nrow(e.dat)*trim)): round(nrow(e.dat)*(1-trim))
  slopes <- lapply(shuffle, function(x){lm(unlist(x[[1]][locs]) ~ locs)$coefficients[2]})

  
  #correct the values
  corrected <- vapply(names(slopes), function(x){
    outlist <- shuffle[[x]][[1]] + ((length(shuffle[[x]][[1]]):1) * slopes[[x]])
    return(outlist[row.names(e.dat)])}, 
    numeric(nrow(e.dat)))
  
  #Back to counts: 
  correctedCounts <- 2^corrected
  correctedCounts <- round(t(t(correctedCounts) * ((lib.size + 1) / 1e+06)) - 0.5)
  
  #update and return the object
  if(class(ExpressionObject) == "ExpressionSet"){
    exprs(ExpressionObject) <- correctedCounts
    ExpressionObject <- ct.normalizeMedians(ExpressionObject, lib.size)
    return(ExpressionObject)
    }else{
      co <- 2^corrected
      new.v <- voom(correctedCounts, design = ExpressionObject$design, normalize.method="scale")
      return(new.v)
      }
}


##' @title Normalize an ExpressionSet Containing a Crispr Screen
##' @description This function normalizes Crispr gRNA abundance estimates contained in an \code{ExpressionSet} object.
##' Currently four normalization methods are implemented: median scaling (via \code{normalizeMedianValues}), slope-based
##' normalization (via \code{ct.normalizeBySlope()}), scaling to the median of the nontargeting control values (via 
##' \code{ct.normalizeNTC()}), and spline fitting to the distribution of the nontargeting gRNAs (via \code{ct.normalizeSpline()}). 
##' Because of the peculiarities of pooled Crispr screening data, these implementations may be more stable than the endogenous methods 
##' used downstream by \link[limma]{voom}. See the respective man pages for further details about specific normalization approaches.
##' @param eset An ExpressionSet object with integer count data extractable with \code{exprs()}.
##' @param method The normalization method to use.
##' @param annotation The annotation object for the library, required for the methods employing nontargeting controls.
##' @param sampleKey An (optional) sample key, supplied as an ordered factor linking the samples to experimental
##' variables. The \code{names} attribute should exactly match those present in \code{eset}, and the control set is assumed to be
##' the first \code{level}. If `method` = `FQ`, the sampleKey is taken as the `sets` argument (and its format requirements are similarly 
##' relaxed; see `?ct.normalizeFC`).
##' @param lib.size An optional vector of voom-appropriate library size adjustment factors, usually calculated with \code{\link[edgeR]{calcNormFactors}} 
##' and transformed to reflect the appropriate library size. These adjustment factors are interpreted as the total library sizes for each sample, 
##' and if absent will be extrapolated from the columnwise count sums of the \code{exprs} slot of the \code{eset}.
##' @param plot.it Logical indicating whether to plot the ranked log2 gRNA count distributions before and after normalization. 
##' @param ... Other parameters to be passed to the individual normalization methods.
##' @return A renormalized ExpressionSet. If specified, the sample level counts will be scaled so as to maintain the validity 
##' of the specified \code{lib.size} values. 
##' @author Russell Bainer
##' @import limma
##' @seealso \code{\link{ct.normalizeMedians}}, \code{\link{ct.normalizeBySlope}}, \code{\link{ct.normalizeNTC}}, \code{\link{ct.normalizeSpline}}
##' @examples data('es')
##' data('ann')
##' 
##' #Build the sample key as needed
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' 
##' es.norm <- ct.normalizeGuides(es, 'scale', annotation = ann, sampleKey = sk, plot.it = TRUE)
##' es.norm <- ct.normalizeGuides(es, 'slope', annotation = ann, sampleKey = sk, plot.it = TRUE)
##' es.norm <- ct.normalizeGuides(es, 'controlScale', annotation = ann, sampleKey = sk, plot.it = TRUE, geneSymb = 'NoTarget')
##' es.norm <- ct.normalizeGuides(es, 'controlSpline', annotation = ann, sampleKey = sk, plot.it = TRUE, geneSymb = 'NoTarget')
##' @export
ct.normalizeGuides <- function(eset, method = c("scale", 'FQ', "slope", "controlScale", "controlSpline"), annotation = NULL, sampleKey = NULL, lib.size = NULL, plot.it = FALSE, ...){
  if(class(eset) != "ExpressionSet"){
    stop(paste(deparse(substitute(eset)), "is not an ExpressionSet."))
  }

  choices <- c("scale", 'FQ', "slope", "controlScale", "controlSpline")
  method <- match.arg(method, choices)

  if(method %in% c("controlScale", "controlSpline")){
    if(is.null(annotation)){
      stop("An annotation object must be provided to perform", deparse(substitute(method)), "normalization.")
    }
    annotation <- ct.prepareAnnotation(
      ann = annotation,
      object =  eset,
      controls = TRUE,
      throw.error = FALSE
    )
  }
  
  if(is.null(lib.size)){
    lib.size <- colSums(exprs(eset))
  } else if(!is.numeric(lib.size) | (length(lib.size) != ncol(eset))){
    stop('If specified, lib.size must be a numeric vector of length equal to the number of samples.')
  }
  
  
  new.eset <- switch(method, 
                     scale = ct.normalizeMedians(eset, lib.size = lib.size), 
                     FQ = ct.normalizeFQ(eset, sets = sampleKey, lib.size = lib.size),
                     slope = ct.normalizeBySlope(eset, lib.size = lib.size, ...), 
                     controlScale = ct.normalizeNTC(eset, annotation, lib.size = lib.size, ...), 
                     controlSpline = ct.normalizeSpline(eset, annotation, lib.size = lib.size, ...)
                    )
  # set negative counts to 0's if they happen to be present after normalization
  exprs(new.eset) <- apply(
    X = exprs(new.eset), 
    MARGIN = 2, 
    FUN = function(col) {
      col[col < 0] <- 0
      col
    })
    
  if(plot.it){
    par(mfrow = c(2, 1))
    ct.gRNARankByReplicate(eset, sampleKey, lib.size = lib.size) 
    ct.gRNARankByReplicate(new.eset, sampleKey, lib.size = lib.size)
   }
   return(new.eset)  
  }


##' @title Normalize sample abundance estimates by median gRNA counts
##' @description This function normalizes Crispr gRNA abundance estimates by equalizing the median gRNA abundance values after
##' correcting for library size. It does this by converting raw count values to log2 counts per million and optionally adjusting further in 
##' the usual way by dividing these values by user-specified library size factors. THis method should be more stable than the endogenous 
##' scaling functions used in \code{voom} in th especific case of Crispr screens or other cases where the median number of observed counts may be low. 
##' @param eset An \code{ExpressionSet} containing, at minimum, count data accessible by \code{exprs}.
##' @param lib.size An optional vector of voom-appropriate library size adjustment factors, usually calculated with \code{\link[edgeR]{calcNormFactors}} 
##' and transformed to reflect the appropriate library size. These adjustment factors are interpreted as the total library sizes for each sample, 
##' and if absent will be extrapolated from the columnwise count sums of the \code{exprs} slot of the \code{eset}.
##' @return A renormalized ExpressionSet object of the same type as the provided object.
##' @author Russell Bainer
##' @import limma
##' @examples data('es')
##' 
##' #Build the sample key and library sizes for visualization
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' ls <- colSums(exprs(es))
##' 
##' es.norm <- ct.normalizeMedians(es, lib.size= ls)
##' ct.gRNARankByReplicate(es, sampleKey = sk, lib.size= ls)
##' ct.gRNARankByReplicate(es.norm, sampleKey = sk, lib.size= ls)
##' @export
ct.normalizeMedians <- function(eset, lib.size = NULL){
  if(!class(eset) == 'ExpressionSet'){
    stop("Please provide an ExpressionSet object for normalization.")
  }
  
  counts <- exprs(eset)
  
  if (is.null(lib.size)){
    lib.size <- colSums(counts)
  }else if(!is.numeric(lib.size) | length(lib.size) != ncol(counts)){
      stop('If specified, lib.size must be a numeric vector of the same length as the number of samples in the eset.')
    } 
  
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  cmed <- apply(y, 2, median, na.rm = TRUE)
  cmed <- cmed - mean(cmed)
  correctedCounts <- 2^t(t(y) - cmed)
  correctedCounts <- (t(t(correctedCounts) * ((lib.size + 1) / 1e+06)) - 0.5)
  
  exprs(eset) <- round(correctedCounts)
  return(eset)
  }

##' @title Factored Quantile Normalization
##' @description This function applies quantile normalization to subsets of samples defined by a provided factor, correcting for library size. 
##' It does this by converting raw count values to log2 counts per million and optionally adjusting further in 
##' the usual way by dividing these values by user-specified library size factors; then this matrix is split into groups according to the 
##' provided factor that are quantile normalized, and then the groups are median scaled to each other before conversion back into 
##' raw counts. This method is best used in comparisons for long timecourse screens, where groupwise differences in growth rate 
##' cause uneven intrinsic dialation of construct distributions.
##' 
##' Note that this normalization strategy is not appropriate for experiments where significant distortion of the libraries is expected as a 
##' consequence of the screening strategy (e.g., strong selection screens). 
##' @param eset An \code{ExpressionSet} containing, at minimum, count data accessible by \code{exprs}.
##' @param sets A character or factor object delineating which samples shoudl be grouped together during the normalization step. Must 
##' be the same length as the number of columns in the provided eset, and cannot contain `NA` or `NULL` values. 
##' @param lib.size An optional vector of voom-appropriate library size adjustment factors, usually calculated with \code{\link[edgeR]{calcNormFactors}} 
##' and transformed to reflect the appropriate library size. These adjustment factors are interpreted as the total library sizes for each sample, 
##' and if absent will be extrapolated from the columnwise count sums of the \code{exprs} slot of the \code{eset}.
##' @return A renormalized ExpressionSet object of the same type as the provided object.
##' @author Russell Bainer
##' @import limma
##' @examples data('es')
##' 
##' #Build the sample key and library sizes for visualization
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' ls <- colSums(exprs(es))
##' 
##' es.norm <- ct.normalizeFQ(es, sets = gsub('(Death|Control)', '', pData(es)$TREATMENT_NAME), lib.size= ls)
##' ct.gRNARankByReplicate(es, sampleKey = sk, lib.size= ls)
##' ct.gRNARankByReplicate(es.norm, sampleKey = sk, lib.size= ls)
##' @export
ct.normalizeFQ <- function(eset, sets, lib.size = NULL){
  if(!class(eset) == 'ExpressionSet'){
    stop("Please provide an ExpressionSet object for normalization.")
  }
  
  stopifnot(length(na.omit(sets)) == ncol(eset), !any(vapply(sets, is.null, logical(1))))
  
  counts <- exprs(eset)
  
  if (is.null(lib.size)){
    lib.size <- colSums(counts)
  }else if(!is.numeric(lib.size) | length(lib.size) != ncol(counts)){
    stop('If specified, lib.size must be a numeric vector of the same length as the number of samples in the eset.')
  } 
  
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  
  #Apply Factored Quantiles
  quant <- lapply(unique(sets), function(x){limma::normalizeQuantiles(y[,(sets == x)])})
  y <- do.call('cbind', quant)
  
  cmed <- apply(y, 2, median, na.rm = TRUE)
  cmed <- cmed - mean(cmed)
  correctedCounts <- 2^t(t(y) - cmed)
  correctedCounts <- (t(t(correctedCounts) * ((lib.size + 1) / 1e+06)) - 0.5)
  
  exprs(eset) <- round(correctedCounts)
  return(eset)
}

##' @title Apply Factored Quantile Normalization to gRNA counts
##' @description This function normalizes Crispr gRNA abundance estimates by equalizing the median gRNA abundance values after
##' correcting for library size. It does this by converting raw count values to log2 counts per million and optionally adjusting further in 
##' the usual way by dividing these values by user-specified library size factors. THis method should be more stable than the endogenous 
##' scaling functions used in \code{voom} in th especific case of Crispr screens or other cases where the median number of observed counts may be low. 
##' @param eset An \code{ExpressionSet} containing, at minimum, count data accessible by \code{exprs}.
##' @param lib.size An optional vector of voom-appropriate library size adjustment factors, usually calculated with \code{\link[edgeR]{calcNormFactors}} 
##' and transformed to reflect the appropriate library size. These adjustment factors are interpreted as the total library sizes for each sample, 
##' and if absent will be extrapolated from the columnwise count sums of the \code{exprs} slot of the \code{eset}.
##' @return A renormalized ExpressionSet object of the same type as the provided object.
##' @author Russell Bainer
##' @import limma
##' @examples data('es')
##' 
##' #Build the sample key and library sizes for visualization
##' library(Biobase)
##' sk <- ordered(relevel(as.factor(pData(es)$TREATMENT_NAME), "ControlReference"))
##' names(sk) <- row.names(pData(es))
##' ls <- colSums(exprs(es))
##' 
##' es.norm <- ct.normalizeMedians(es, lib.size= ls)
##' ct.gRNARankByReplicate(es, sampleKey = sk, lib.size= ls)
##' ct.gRNARankByReplicate(es.norm, sampleKey = sk, lib.size= ls)
##' @export
ct.normalizeFactoredQuantiles <- function(eset, lib.size = NULL){
  if(!class(eset) == 'ExpressionSet'){
    stop("Please provide an ExpressionSet object for normalization.")
  }
  
  counts <- exprs(eset)
  
  if (is.null(lib.size)){
    lib.size <- colSums(counts)
  }else if(!is.numeric(lib.size) | length(lib.size) != ncol(counts)){
    stop('If specified, lib.size must be a numeric vector of the same length as the number of samples in the eset.')
  } 
  
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  cmed <- apply(y, 2, median, na.rm = TRUE)
  cmed <- cmed - mean(cmed)
  correctedCounts <- 2^t(t(y) - cmed)
  correctedCounts <- (t(t(correctedCounts) * ((lib.size + 1) / 1e+06)) - 0.5)
  
  exprs(eset) <- round(correctedCounts)
  return(eset)
}





