#This document contains functions for the RRAa algorithm. 

##' @title Checks and Possibly Sets the Number of Cores to be Used in Parallel Processing
##' @description This function determines the number of cores that the user is expecting to 
##' use during parallel processing operations, and if absent, sets the \code{mc.cores} option
##' to the maximum value. Users who do not wish to use all available cores during parallel 
##' processing should do so by invoking \code{options()} from the command line prior to analysis. 
##' 
##' @return Nothing, but invisibly sets \code{options(mc.cores)} if currently \code{NULL}. 
##' @author Russell Bainer, Pete Haverty
ct.numcores <- function()  {
    if(is.null(getOption('mc.cores'))){
       numcores = parallel::detectCores()
       options(mc.cores = numcores)
   }
    invisible()
  }

##' @title Aggregation of P-value Ranks using a Beta Distribution and Alpha Cutoff  
##' @description This function calculates an alpha-modified rho statistic from a set of normalized ranks by comparing them to a uniform distribution. 
##' Specifically, the ranks are ordered and p-values calculated at each position in the ordered vector by comparison to a Beta distribution. The rho value 
##' returned is the smallest p-value identified in this way across all scores. Should not be used by end users.
##' @param p.in A single column matrix of rank scores, with row.names indicating the gRNA labels.
##' @return A numeric rho value corresponding to the minimum rank order P. 
##' @author Russell Bainer, modified from code from the \code{RobustRankAggreg} package. 
##' 
##' Citation: 
##' Kolde, R. et al, Bioinformatics. 2012 Feb 15;28(4):573-80. doi: 10.1093/bioinformatics/btr709.
##' @import RobustRankAggreg
##' @examples
##' testp <- runif(20)
##'
##' RobustRankAggreg::betaScores((testp))
##' ct.alphaBeta(testp)  
##' @export

ct.alphaBeta <- function(p.in){ 
  p.in <- na.omit(p.in)
  n <- length(p.in)  
  if(n == 0){
    return(1)
  } else {
    p.in <- sort(p.in)
    return(min(pbeta(p.in, 1:n, n - 1:n + 1)))
  }
}


##' @title Aggregation of P-value Ranks using a Beta Distribution and Alpha Cutoff  
##' @description This function is called internally as a single instance of the beta aggregation step in RRAa. Users should not interact with it directly. 
##' The expected input is a list of rank statistics, and a paired \code{alpha} argument defining which values to consider in downstream analyses (see below). 
##' @param p A single column matrix of rank statistics, with \code{row.names} indicating the gRNA labels. 
##' @param g.key data.frame with guide and gene names
##' @param alpha The alpha cutoff parameter, corresponding to the P-value threshold or fold change proportion at which gRNAs should no longer be considered to be 
##' differentially expressed. Alternatively, this can be provided as a logical vector of the same length as the number of rows in \code{p}, 
##' containing only \code{TRUE} and \code{FALSE} elements indicating whether the element should be included during the aggeregation step. 
##' @param shuffle Logical indicating whether to shuffle the rank statistics prior to calculating the rho statistics (useful for permutation).
##' @param return.obj Name of the environment to record results in, or \code{TRUE} to return the RRAa results directly as a numeric vector of genewise rho statistics. 
##' If an environment is supplied, the function will directly increment the \code{target.positive.iterations} variable within it on the basis of the \code{obs} 
##' variable in the specified \code{environment}. 
##' @return Nothing, or a named list of target-level P-values, which are treated as a rho statistic in the permutation step. 
##' @author Russell Bainer
##' @examples 
##' data('fit')
##' data('ann')
##' geneScores <- ct.RRAalpha(fit$p.value, ann, alpha = 0.1, shuffle = FALSE, return.obj = TRUE)
##' @export

ct.RRAalpha <- function(p, g.key, alpha, shuffle = FALSE, return.obj = TRUE){
  
  #determine the input format and convert to a logical as needed. 
  if(length(alpha) == 1){

    ##tally the significant genes
    pass <- p[,1] <= alpha

    }else{
      if((sum(alpha %in% c(TRUE, FALSE)) != nrow(p))){stop('When provided as a vector, alpha values must only be TRUE or FALSE 
                                                           and exactly match the rows of the provided rank statistic.')}
      pass <- alpha  
    }
  
  #nonsignificant gRNAs are set to 1. Implicitly, this means that the significance 
  #relative to alpha is treated as an inherent property of the gRNA. 
  p[!pass,1] <- 1

  if(shuffle){
    p <- sample(p, nrow(p))
  }

  p.collect <- split(p, g.key$geneSymbol)
  rhoscores <- vapply(p.collect, ct.alphaBeta, numeric(1))
  
  if(!is.environment(return.obj) | deparse(substitute(return.obj)) %in% 'none'){
    return(rhoscores)
  } else {
     observations <- return.obj[['obs']]
     passcount <- return.obj[['target.positive.iterations']]
     passcount[rhoscores <= observations] <- passcount[rhoscores <= observations] + 1
     return.obj[['target.positive.iterations']] <- passcount
     invisible()
  }
}  



##' @title gRNA signal aggregation via RRAa, optionally using multiple cores. 
##' @description This is a wrapper function implementing the RRAalpha p-value aggregation algorithm. Takes in a set of gRNA rank scores (formatted as a single-column 
##' numeric matrix with row.names indicating the guide names) and a list object of gRNA annotations (names are the gene targets, and each element of the list contains 
##' a vector of the corresponding guide names). The rank scores are converted to gene-level statistics that are thenm transformed into empirical p-values by permutation. 
##' @param p A single column matrix of ranking scores, with row.names indicating the gRNA labels
##' @param g.key An annotation data frame of gRNAs, minimally containing a factorized "geneSymbol" column indicating the target names. This is typically generated by calling the \code{ct.buildKeyFromAnnotation()} function.  
##' @param alpha The alpha cutoff parameter, corresponding to the P-value threshold or fold change proportion at which gRNAs should no longer be considered to be 
##' differentially expressed. Alternatively, this can be provided as a logical vector of the same length as the number of rows in \code{p}, 
##' containing only \code{TRUE} and \code{FALSE} elements indicating whether the element should be included during the aggeregation step. 
##' @param permute Number of permutations to be used during empirical p-value estimation. In a multicore context the exact number of permutations may vary somewhat
##' to accomodate the corresponding system archetecture but should be close to the specified permutation number. 
##' @param multicore Logical indicating whether to use multiple cores to calculate p-values.
##' @param core.perm Maximum number of permutations to run on each core (only relevant when \code{multicore} is \code{TRUE}). 
##' @param permutation.seed numeric seed for permutation reproducibility.
##'   Default: \code{NULL} means to not set any seed.
##' @return A named list of target-level empirical P-values. 
##' @author Russell Bainer
##' @examples 
##' data('fit')
##' data('ann')
##' genePvals <- ct.RRAaPvals(fit$p.value, ann, alpha = 0.1, permute = 100, multicore = FALSE)
##' genePvals <- ct.RRAaPvals(fit$p.value, ann, alpha = 0.1, permute = 100, multicore = TRUE, core.perm = 10)
##' @export
ct.RRAaPvals <- function(p,
                         g.key,
                         alpha,
                         permute,
                         multicore = TRUE,
                         core.perm = 100,
                         permutation.seed = NULL) {
  
  #Input checks
  if(class(p) != "matrix" | is.null(ncol(p))){
    stop("P-values should be input as a single-column matrix with row names contained in the gs.list")}
  if(ncol(p) > 1){
    warning(paste('Multiple columns are present in the p-value matrix. Using the first column:', colnames(p)[1]))}
  if(!is.numeric(core.perm) | length(core.perm) > 1){
    stop('core.perm must be supplied as a single numeric value.')
  }
  
  if(length(alpha) > 1){
    
    if((length(alpha) != nrow(p)) | (sum(alpha %in% c(TRUE, FALSE)) != nrow(p))){
      stop("When provided as a vector, alpha must contain exactly one element for each row of the p object, and only contain TRUE or FALSE values.")
    }

  }else{      
    if(!is.numeric(alpha) | (length(alpha) != 1)){
      stop("If not provided as a list, the alpha parameter must be a single numeric value.")
    }
  }
  
  if(!is.data.frame(g.key)){stop("The annotation provided must be a data frame.")}
  if(!("geneSymbol" %in% names(g.key))){stop("The provided annotation does not contain a geneSymbol column.")}
  if(!setequal(row.names(g.key), row.names(p))){stop("Provided p-value list and annotation object contain different elements.")}
  
  is.null(permutation.seed) ||
    is.numeric(permutation.seed) ||
    stop("'permutation.seed' must be numeric or NULL.")
  
  #Everything apparently ok, generate P-values. 
  result.environment <- new.env()

  result.environment$obs <- ct.RRAalpha(p, g.key, alpha, return.obj = 'none')
  nguides <- nrow(p)
  resgenes <- levels(g.key$geneSymbol)
  result.environment$ngenes <- length(resgenes)
  result.environment$target.positive.iterations <- rep(0, result.environment$ngenes)

  if(multicore){
    #figure out the right number of jobs to send to each core
    ct.numcores()
    cores <- getOption('mc.cores')

    jpc <- min(floor(permute/cores), core.perm)
    if(jpc == 0){jpc <- 1}
    njobs <- ceiling(permute/jpc)
    permute <- jpc * njobs 
    
    if (is.null(permutation.seed)) {
      batch.perm.seeds <- NULL
    } else {
      # generate a seed for each permutation batch
      batch.perm.seeds <-
        seq(from = permutation.seed, to = permutation.seed + njobs - 1)
    } 
    
    message(paste("Permuting", permute, 'times, using', cores, 'cores.'))
    out <-
      mclapply(1:njobs, function(x) {
        ct.RRAalphaBatch(
          p,
          g.key,
          alpha,
          result.environment,
          batch.size = jpc,
          permutation.seed = batch.perm.seeds[x] 
        )
      }, mc.preschedule = FALSE, mc.cores = cores) 
    result.environment$target.positive.iterations <- rowSums(as.data.frame(out))
    } else {
      set.seed(permutation.seed) # default NULL will have no effect
      message(paste("Permuting", permute, 'times, this may take a while...'))
      iter <- (replicate(permute, ct.RRAalpha(p, g.key, alpha, shuffle = TRUE, return.obj = 'none'))) 
      result.environment$target.positive.iterations <- unlist(lapply(1:nrow(iter), function(x){sum(iter[x,] <= result.environment$obs[x])}))
      }

  out <- result.environment$target.positive.iterations/permute
  names(out) <- names(result.environment$obs)
  return(out)
}

##' @title Create Batches of Null Permutations for a Crispr Screen
##' @description This is a wrapper function to partition batches of calls to \code{ct.RRAalpha()} for multicore processing. 
##' It is called internally as a single instance of the beta aggregation step in RRAa. Users should not interact with it directly. 
##' @param p A single column matrix of rank statistics, with row.names indicating the gRNA labels. 
##' @param g.key data.frame with guide and gene names
##' @param alpha The alpha cutoff parameter, corresponding to the P-value threshold or fold change proportion at which gRNAs should no longer be considered to be 
##' differentially expressed. Alternatively, this can be provided as a logical vector of the same length as the number of rows in \code{p}, 
##' containing only \code{TRUE} and \code{FALSE} elements indicating whether the element should be included during the aggeregation step. 
##' @param result.environment The target environment containing the quasi-global variables incremented during the permutations in the child functions. 
##' @param batch.size Number of iterations to deploy to each daughter process. 
##' @param permutation.seed numeric seed for permutation reproducibility. Default is \code{NULL}, in which case no seed is set.
##' @return An integer vector indicating the number of iterations in which each gene's score was better than those indicated in \code{result.environment$obs}. 
##' @author Russell Bainer
##' @examples 
##' data('fit')
##' data('ann')
##' batch.size <- 50
##' env <- new.env()
##' env$obs <- ct.RRAalpha(fit$p.value, ann, alpha = 0.1, return.obj = 'none')
##' env$ngenes <- length(levels(ann$geneSymbol))
##' passed <- ct.RRAalphaBatch(fit$p.value, ann, alpha = 0.1, result.environment = env, batch.size = batch.size)
##' hist(passed/batch.size, main = 'Empirical P-value', xlab = 'P')
##' @export

ct.RRAalphaBatch <- function(p, g.key, alpha, result.environment, batch.size = 100, 
                             permutation.seed = NULL){
  
  #make a batch environment
  batch.env <- new.env()
  batch.env$obs <- result.environment$obs
  batch.env$target.positive.iterations <- rep.int(0, result.environment$ngenes)
  
  #run permutations and increment as needed
  set.seed(permutation.seed) # default NULL will have no effect
  invisible(replicate(batch.size, ct.RRAalpha(p, g.key, alpha, shuffle = TRUE, return.obj = batch.env))) 
  
  return(batch.env$target.positive.iterations)
}
