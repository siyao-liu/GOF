#' Fit a distribution to gene count distribution
#'
#' @description This function fits either a Negative Binonimal (NB) or Average Negative Binomial (ANB) distribution to the observed gene count distribution.
#'
#' @param vdata A vector of data points - gene count data.
#' @param family A string specifying the model type to be used for fitting the observed distribution. By default, we fit an Average Negative Binomial distribution to the observed gene count distribution.
#' @param p.threshold Probability threshold used to find the big quantile of the NB or ANB fit.
#' @param common.sf The common spread factor parameter. By default, it is 0.7440919 which was estimated from the 3 cell line mixture data.
#' @return FitDist returns a list with containing the following components:
#' \itemize {
#' \item Data: a vector of empirical gene count data.
#' \item FittedData: a vector of fitted gene count data.
#' \item vigrid: a vector of grid from 0,1,2,..., maximum of data point.
#' \item vpi: a vector of empirical p_i probabilities.
#' \item vqi: a vector of fitted q_i probabilities.
#' \item nmax: the maximum of the biggest count in p_i and the biggest count in q_i.
#' }
#' @export
FitDist <- function(vdata,
                    family="average negative binomial",
                    p.threshold=1-10^(-5),
                    common.sf=0.7440919 # threshold derived from the 3 cell line original data (may subject to change)
) {

  # turn into empirical probabilities
  # generate grid 0,1,2,..., max data point
  vpigrid <- seq(0, max(vdata), 1)

  vcounts <- hist(vdata, breaks=c(vpigrid, (max(vpigrid)+1)), plot=FALSE, right=FALSE)$counts

  # vector of empirical p_i probabilities
  vpi <- vcounts / length(vdata)

  # Fit Negative Binomial
  # fit Negative Binomial parameters
  require(fitdistrplus)
  nbfit <- fitdistrplus::fitdist(data = vdata, distr = "nbinom", method = "mle",
                                 control=list(trace=1, REPORT=1), lower=c(0, 0))
  nbParm <- nbfit$estimate

  if (family=="negative binomial") {
    # Big quantile of the NB fit
    qbig <- qnbinom(p=p.threshold, size = nbParm[1], mu = nbParm[2])
    # Theoretical quantile of the NB fit
    qdata <- qnbinom(p=ppoints(length(vdata)), size = nbParm[1], mu = nbParm[2])

    if (qbig > max(vdata)) { # q igrid is bigger than p igrid, then pad vpi with 0s
      nmax <- qbig
      vigrid <- seq(0, nmax, 1)
      vpi <- c(vpi, rep(0, (nmax+1-length(vpigrid))))
    } else { # q grid no bigger than p grid, then just use that
      nmax <- length(vpigrid)
      vigrid <- vpigrid
    }
    # vector of NB q_i probabilities
    vqi <- dnbinom(x=vigrid, size = nbParm[1], mu = nbParm[2])

  }

  if (family=="average negative binomial") {

    # Fit Average Negative Binomial
    # find big quantile of the ANB fit
    qbig <- qnbinom(p=p.threshold, size = (1 / (common.sf)^2), mu = nbParm[2])
    # Theoretical quantile of the ANB fit
    qdata <- qnbinom(p=ppoints(length(vdata)), size = (1 / (common.sf)^2), mu = nbParm[2])

    if (qbig > max(vdata)) { # q igrid is bigger than p igrid, then pad vpi with 0s
      nmax <- qbig
      vigrid <- seq(0, nmax, 1)
      vpi <- c(vpi, rep(0, (nmax+1-length(vpigrid))))
    } else { # q grid no bigger than p grid, then just use that
      nmax <- length(vpigrid)
      vigrid <- vpigrid
    }

    # vector of ANB q_i probabilities
    vqi <- dnbinom(x=vigrid, size = (1 / (common.sf)^2), mu = nbParm[2])
  }

  return(list("Data"=vdata,
              "FittedData"=qdata,
              "vigrid"=vigrid,
              "vpi"=vpi,
              "vqi"=vqi,
              "nmax"=nmax))
}
