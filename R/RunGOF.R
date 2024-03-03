#' Run Goodness of Fit (GOF) methods for feature selection.
#'
#' @description Take a raw count matrix (gene in the row and cell in the column) and the NB or ANB fitted output, runs one of the 4 GOF variants for feature selection.
#' The output of this function is a list that contains the goodness of fit measures for all genes in the count matrix and top.n features selected based on the goodness of fit
#' measures. Note: both Wdist.mean and Wdist.med methods further filter the genes based on their sample and theoretical variance and only keep the genes with their sample variance
#' greater than their theoretical variance.
#'
#' @param counts A raw count matrix with genes in the row and cells in the column.
#' @param countsFit A list generated after running FitDist().
#' @param top.n Number of features to select.
#'     (default: 2000)
#' @return Depending on the choice of method, RunGOF returns a list that contains
#' \itemize {
#' \item PP.ABC: a vector of ABC values for each gene using the P-P method.
#' \item topPP: top features selected from P-P method.
#' \item QQ.ABC: a vector ABC values for each gene using the Q-Q method.
#' \item topQQ: top features selected from Q-Q method.
#' \item Wdist.mn: a vector of 1-Wasserstein distance for each gene adjusted by gene mean expression.
#' \item topWdist.mn: top features selected from 1-Wasserstein distance adjusted by gene mean expression method.
#' \item Wdist.med: a vector of 1-Wasserstein distance for each gene adjusted by gene median expression.
#' \item topWdist.med: top features selected from 1-Wasserstein distance adjusted by gene median expression method.
#' }
#'
#' @export
RunGOF <- function(counts,
                   countsFit,
                   method = c("PP", "QQ", "Wdist.mean", "Wdist.med"),
                   top.n = 2000) {

  if (method=="PP") {
    ### run PP ####
    pp.abc <- unlist(lapply(countsFit, function(x) {
      calcPPabc(x$vpi, x$vqi) }))
    toppp <- names(sort(pp.abc, decreasing=T)[1:top.n])
    return(list("PP.abc"=sort(pp.abc, decreasing=TRUE), "topPP"=toppp))
  }
  if (method=="QQ") {
    ### run QQ ###
    qq.abc <- unlist(lapply(countsFit, function(x) {
      calcQQabc(x$Data, x$FittedData, "Observed", "Expected") }))
    topqq <- names(sort(qq.abc, decreasing=T)[1:top.n])
    return(list("QQ.abc"=sort(qq.abc, decreasing=TRUE), "topQQ"=topqq))
  }
  if (method=="Wdist.mean") {
    ### run 1-Wasserstein ###
    wasser.out <- sapply(names(countsFit), function(x)
      calcDiscWasser1(vdata=countsFit[[x]]$Data,
                      vpi=countsFit[[x]]$vpi,
                      vqi=countsFit[[x]]$vqi,
                      nmax=countsFit[[x]]$nmax,
                      method="mean"),
      simplify=FALSE, USE.NAMES = TRUE)
    qdata.list <- lapply(countsFit, '[[', 2)
    Wdist.mean <- unlist(lapply(wasser.out, '[[', 2))

    # elimiate genes with sample variance < theoretical variance
    gene.obs.var <- apply(counts, 1, var)
    gene.exp.var <- unlist(lapply(qdata.list, function(x) var(x)))
    filt <- which(gene.obs.var > gene.exp.var)
    Wdist.mean <- Wdist.mean[filt]
    topWdist.mn <- names(sort(Wdist.mean, decreasing=TRUE))[1:top.n]
    return(list("Wdist.mn"=sort(Wdist.mean, decreasing=TRUE), "topWdist.mn"=topWdist.mn))
  }
  if (method=="Wdist.med") {
    ### run 1-Wasserstein ###
    wasser.out <- sapply(names(countsFit), function(x)
      calcDiscWasser1(vdata=countsFit[[x]]$Data,
                      vpi=countsFit[[x]]$vpi,
                      vqi=countsFit[[x]]$vqi,
                      nmax=countsFit[[x]]$nmax,
                      method="median"),
      simplify=FALSE, USE.NAMES = TRUE)
    qdata.list <- lapply(countsFit, '[[', 2)
    Wdist.med <- unlist(lapply(wasser.out, '[[', 2))
    # elimiate genes with sample variance < theoretical variance
    gene.obs.var <- apply(counts, 1, var)
    gene.exp.var <- unlist(lapply(qdata.list, function(x) var(x)))
    filt <- which(gene.obs.var > gene.exp.var)
    Wdist.med <- Wdist.med[filt]
    topWdist.med <- names(sort(Wdist.med, decreasing=TRUE))[1:top.n]
    return(list("Wdist.med"=sort(Wdist.med, decreasing=TRUE), "topWdist.med"=topWdist.med))
  }

}

