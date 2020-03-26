##' @title Compare Two CRISPR Screens via a CAT plot 
##' @description This is a function for comparing the results of two screening experiments. Given two \code{summaryDF}, 
##' the function places them in register with one another, generates a Concordance At The Top (CAT) plot, and returns an 
##' invisible dataframe containing the relevant gene-level signals.
##' 
##' This function is conceptually similar to the `ct.ROC` and `ct.PRC()` functions, but is appropriate when considering 
##' consistency of ranked values rather than an interchangeable set; the most common use case is for comparing primary 
##' and replication screens, where the underlying technology and selection criteria are expected to be highly similar. 
##' CAT plots are fundamentally about comparing rankings, and so only targets in common between the two provided 
##' screens are considered. If the totality of list overlap is important, consider using `ct.PRC()` or `ct.ROC()`.
##' 
##' Note that ranking statistics in CRISPR screens are (usually) permutation-based, and so some granularity in the 
##' rankings is expected. This function does a little extra work to limit the influence of this granularity and to 
##' ensure that hits are counted as soon as the requisite value of the ranking statistic is reached regardless of 
##' where the target is located within the block of equally-significant hits. Functionally, this means that the 
##' drawn curve is somewhat anticonservative in cases where the target ranks are not well differentiated.  
##'
##' @param df1 A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param df2 A dataframe summarizing the results of the screen, returned by the function \code{\link{ct.generateResults}}. 
##' @param targets Column of the provided \code{summaryDF} to consider. Must be \code{geneID} or \code{geneSymbol}.
##' @param enrich Logical indicating whether to test for enrichment or depletion.
##' @param plot.rho Logical indicating whether to plot the Rho values in addition to the P-values, which sometimes 
##' have better ranking properties. 
##' @return Invisibly, a data.frame containing the relevant summary stats for each target in both screens. 
##' @author Russell Bainer
##' @examples 
##' data('resultsDF')
##' cat <- ct.CAT(resultsDF, resultsDF[1:2000,], enrich = TRUE)
##' head(cat)
##' @export
ct.CAT <-
  function(df1, df2,
           targets = c('geneSymbol', 'geneID'),
           enrich, 
           plot.rho = TRUE) {

    #Check the input: 
    if(any(!ct.resultCheck(df1), !ct.resultCheck(df2))){
      stop("Invalid summary dataframe provided. Execution halted.")
    }
    targets <- match.arg(targets)
    stopifnot(is.logical(enrich), is.logical(plot.rho))
    
    #Convert to gene-level stats & place in register. 
    df1 <- df1[!duplicated(df1[,targets]),]
    row.names(df1) <- df1[,targets]
    df2 <- df2[!duplicated(df2[,targets]),]
    row.names(df2) <- df2[,targets]
    
    out <- data.frame('Target' = union(df1[,targets], df2[,targets]))
    row.names(out) <- out$Target
    out$screen1_enrich.rho <- df1[out$Target,'Rho_enrich']
    out$screen1_enrich.p <- df1[out$Target,'Target-level Enrichment P']
    out$screen1_deplete.rho <- df1[out$Target,'Rho_deplete']
    out$screen1_deplete.p <- df1[out$Target,'Target-level Depletion P']
    out$screen2_enrich.rho <- df2[out$Target,'Rho_enrich']
    out$screen2_enrich.p <- df2[out$Target,'Target-level Enrichment P']
    out$screen2_deplete.rho <- df2[out$Target,'Rho_deplete']
    out$screen2_deplete.p <- df2[out$Target,'Target-level Depletion P']
    
    common <- intersect(row.names(df1), row.names(df2))
    
    message(paste0(length(common), ' ', targets, ' elements in common.'))
    
    ranks <- apply(out[common,2:9], 2, base::rank, na.last = TRUE, ties.method = 'min')

    if(enrich){
      message('Sorting on enrichment.')
      ranks <- ranks[,c(1,5,2,6)]
      out <- out[order((-log10(out$screen1_enrich.rho) + -log10(out$screen2_enrich.rho))/2, decreasing = TRUE),]
    } else {
      message('Sorting on depletion.')
      ranks <- ranks[,c(3,7,4,8)]
      out <- out[order((-log10(out$screen1_deplete.rho) + -log10(out$screen2_deplete.rho))/2, decreasing = TRUE),]
    }
    
    message(paste0(min(apply(ranks,2, function(x){length(unique(x))})), ' distinct ranks identified.'))
    
    #reorder to pick the least granular statistic as the dependent variable.
    if(length(unique(ranks[,1])) < length(unique(ranks[,2]))){ ranks[,1:2] <- ranks[,2:1]}
    if(length(unique(ranks[,3])) < length(unique(ranks[,4]))){ ranks[,3:4] <- ranks[,4:3]}
    
    #Compose statistics. more significant = Bigger numbers. 
    rho.cat <- t(vapply(unique(ranks[,1]), 
                        function(x){
                          c(sum(ranks[(ranks[,1] <= x),2] <= x), sum(ranks[,1] <= x)) 
                        }, 
                        numeric(2)))
    p.cat <- t(vapply(unique(ranks[,3]), 
                        function(x){
                          c(sum(ranks[(ranks[,3] <= x),4] <= x), sum(ranks[,3] <= x)) 
                        }, 
                        numeric(2)))
    p.cat <- p.cat[order(p.cat[,2], decreasing = TRUE),] 
    rho.cat <- rho.cat[order(rho.cat[,2], decreasing = TRUE),] 
    
    #Plot it
    plot((p.cat[,2]/nrow(ranks)), p.cat[,1]/p.cat[,2], 
         main = paste0(ifelse(enrich, 'Enrichment ', 'Depletion '), 'CAT'), 
         xlab = 'Signal Rank', ylab = 'Concordance', 
         xlim = c(0,1), ylim = c(0,1),
         pch = 19, col = rgb(0,0,0.7))
    abline(0,1,col='red')
    if(plot.rho){lines(rho.cat[,2]/nrow(rho.cat), rho.cat[,1]/rho.cat[,2], 
                       lty = 2, col = 'darkgrey', lwd = 2)
      legend('bottomright', c('P', 'Rho'), fill = c(rgb(0,0,0.7), 'darkgrey'))
      }
    return(invisible(out))
  }













