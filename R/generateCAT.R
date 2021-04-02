##' @title Compare Two CRISPR Screens via a CAT plot 
##' @description This is a function for comparing the results of two screening experiments. Given two \code{summaryDF}, 
##' the function places them in register with one another, generates a Concordance At The Top (CAT) plot, and returns an 
##' invisible dataframe containing the relevant gene-level signals. Signals are aggregated by P-value threshold, such 
##' that the concordance is represented as the portion of shared values meeting or exceeding that significance threshold.
##' 
##' This function is conceptually similar to the `ct.ROC` and `ct.PRC()` functions, but is appropriate when considering 
##' consistency of ranked values rather than an interchangeable set; the most common use case is for comparing primary 
##' and replication screens, where the underlying technology and selection criteria are expected to be highly similar. 
##' CAT plots are fundamentally about comparing rankings, and so only targets in common between the two provided 
##' screens are considered. If the totality of list overlap is important, consider using `ct.PRC()` or `ct.ROC()`.
##'
##' @param dflist A list of results dataframes. Names will be preserved, and the enrichment calculation is conditioned on 
##' the first element of the list.
##' @param targets Column of the provided \code{summaryDF} to consider. Must be \code{geneID} or \code{geneSymbol}.
##' @param switch.dir Logical indicating whether to test overlap of signals in the same direction, or whether the 
##' directionality is expected to reverse. `same.dir = FALSE` looks at the consistency between depleted signals in `df1` and 
##' enriched signals in `df2`.
##' @param plot.it Logical indicating whether to compose the plots on the default device. Two CAT plots summarizing overlap in 
##' both enrichment directions are drawn. 
##' @return Invisibly, a data.frame containing the relevant summary stats for each target in both screens. 
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' cat <- ct.CAT(list('first' = resultsDF, 'second' = resultsDF[1:2000,]))
##' head(cat)
##' @export
ct.CAT <-
  function(dflist,
           targets = c('geneSymbol', 'geneID'),
           switch.dir = FALSE, 
           plot.it = TRUE) {

    #Check the input:
    stopifnot(is(switch.dir, 'logical'), is(plot.it, 'logical'), is(dflist, 'list'), length(dflist) == 2)
    targets <- match.arg(targets)
    dflist <- ct.regularizeContrasts(dflist, collapse = targets)
    titles <- paste0(names(dflist)[2], ' signal by ', names(dflist)[1], ', ', c('Enrichment', 'Depletion'))
    cuts <- c(min(sum(dflist[[1]]$direction %in% 'enrich'), sum(dflist[[2]]$direction %in% 'enrich')), 
              min(sum(dflist[[2]]$direction %in% 'deplete'), sum(dflist[[2]]$direction %in% 'deplete')))
    cuts <- cuts/nrow(dflist[[1]])

    #Order is out, rank is in
    r1 <- ct.rankSimple(dflist[[1]])
    
    if(switch.dir){
      titles <- paste0(names(dflist)[2], c(' Depletion', ' Enrichment'), ' by ', names(dflist)[1], c(' Enrichment', ' Depletion'))
      r2 <- ct.rankSimple(dflist[[2]], top = 'deplete')
    } else {
      r2 <- ct.rankSimple(dflist[[2]])
    }

    #Calculate overlaps by rank
    enrich.over <- vapply(sort(unique(r1), decreasing = FALSE)[],
                          function(x){
                            return(sum(r2[r1 <= x] <= x)/sum(r1 <= x))
                            },
                          numeric(1))

    #standardize the backend
    r2 <- r2 + (max(r1) - max(r2))
    
    deplete.over <- vapply(sort(unique(r1), decreasing = TRUE),
                          function(x){
                            return(sum(r2[r1 >= x] >= x)/sum(r1 >= x))
                          },
                          numeric(1))
    #Plot it
    if(plot.it){
      
      plot(sort(unique(r1), decreasing = FALSE)/max(r1), 
           enrich.over, 
           main = titles[1], 
           xlab = paste0(names(dflist)[1], ' Signal Rank'), ylab = 'Concordance', 
           xlim = c(0,cuts[1]), ylim = c(0,1),
           pch = 19, col = rgb(0,0,0.7), 
           type = 'l')
      abline(0,1,col='red')
      
      plot(sort(unique(r1), decreasing = FALSE)/max(r1), 
           deplete.over, 
           main = titles[2], 
           xlab = paste0(names(dflist)[1], ' Signal Rank'), ylab = 'Concordance', 
           xlim = c(0,cuts[2]), ylim = c(0,1),
           pch = 19, col = rgb(0,0,0.7), 
           type = 'l')
      abline(0,1,col='red')
    }
    return(invisible(data.frame('Enrich'=enrich.over, 'Deplete'=deplete.over)))
  }
