#' Compute 1-Wasserstein distance between the fitted gene count distribution and the observed gene count distribution
#' @description This function computes the 1-Wasserstein distance between two discrete distributions.
#'
#' @param vdata A numeric vector of the empirical gene count data.
#' @param vpi A numeric vector of the empirical probabilities.
#' @param vqi A numeric vector of the theoretical probabilities.
#' @param nmax The maximum of the biggest count in vp_i and the biggest count in vq_i.
#' @param method The method used to adjust the 1-Wasserstein distance metric.
#' @return calcDiscWasser1 returns a list that contains
#' \itemize {
#' \item Wdist: a vector of the unadjusted 1-Wasserstein distance for each gene.
#' \item Wdist.adj: a vector of adjusted 1-Wasserstein distance for each gene.
#' \item vDj: a vector of the difference between the two quantile functions.
#' }
#' @export
calcDiscWasser1 <- function(vdata,
                            vpi, # empirical probability
                            vqi,  # theoretical probability
                            nmax,
                            method=c("mean", "median")) {
  # initialization
  ptilde <- vpi[1]
  qtilde <- vqi[1]
  Wdist <- 0
  vDj <- c()

  # iteration: looping through integers 1,2,...,nmax
  for (j in 1:nmax) {
    Dj <- ptilde - qtilde
    Cj <- abs(Dj)

    if (Dj >= 0) {
      ptilde <- vpi[j+1] + Cj
      qtilde <- vqi[j+1]
    } else {
      ptilde <- vpi[j+1]
      qtilde <- vqi[j+1] + Cj
    }

    Wdist <- Wdist + Cj
    vDj[j] <- Dj

    if (method=="mean") {
      Wdist.adj <- Wdist/mean(vdata)
    }
    if (method=="median") {
      Wdist.adj <- Wdist/(median(vdata) + 1)
    }
  }

  return(list("Wdist"=Wdist,
              "Wdist.adj"=Wdist.adj,
              "vDj"=vDj))
}

