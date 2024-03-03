#' Q-Q plot comparing two samples with a few integer values.
#' @description A function for generating a Q-Q plot to assess whether two samples with a few integer values come from the same distribution.
#'
#'
#' @param P A numeric vector of sampled data points from the first distribution.
#' @param Q A numeric vector of sampled data points from the second distribution.
#' @param sample1 A character to denote the first distribution.
#' @param sample2 A character to denote the second distribution.
#' @return qqplot_small_test returns a Q-Q plot using ggplot().
#'
#' @examples
#' P <- sample(0:2,1000, replace=TRUE, prob=c(1/3, 1/2, 1/6))
#' Q <- sample(0:1,1000, replace=TRUE, prob=c(2/3, 1/3))
#' qqplot_small_test(P, Q, "P", "Q")
#'
#' @references {
#' Pan, Y., Landis, J.T., Moorad, R. et al. The Poisson distribution model fits UMI-based single-cell RNA-sequencing data.
#' \emph{BMC Bioinformatics 24, 256 (2023).}
#' \url{https://doi.org/10.1186/s12859-023-05349-2}
#' }
#' @export
qqplot_small_test <- function(P, Q, sample1, sample2) {
  require(ggplot2)
  require(dplyr)
  new_quantile <- function(data, sample){

    if(!is.numeric(data)){
      warning("need numeric values for input data")
    }

    dfp <- data.frame(value = data)
    dfp <- dfp %>%
      do(data.frame(., fval = ecdf(.$value)(.$value)))
    dfp <- dfp %>% distinct()
    dfp <- dfp[with(dfp, order(fval)),]
    dfp$cat <- sample
    dfp <- rbind(data.frame(dfp), data.frame(value = c(min(data)-1), cat = sample, fval = c(0)))
    dfp$value_new <- dfp$value + 1/2
    dfp <- dfp[with(dfp, order(fval)),]
    return(dfp)
  }

  qq_interpolation <- function(dfp, dfq, sample1, sample2){

    `%notin%` <- Negate(`%in%`)

    if (!all(colnames(dfp) == c("value", "fval", "cat", "value_new"))){
      warning("column names for first data not match, run new_quantile() function first")
    }

    if (!all(colnames(dfq) == c("value", "fval", "cat", "value_new"))){
      warning("column names for second data not match, run new_quantile() function first")
    }

    if (dfp$fval[1] != 0){
      warning("fvalue for first data not start from 0, run new_quantile() function first")
    }

    if (dfq$fval[1] != 0){
      warning("fvalue for second data not start from 0, run new_quantile() function first")
    }

    # interpolation
    for (i in 1:nrow(dfp)){
      if (dfq$fval[1] == 0 &dfp$fval[i] < 1 & dfp$fval[i] %notin% dfq$fval){
        cat_add <- sample2
        fval_add <- dfp$fval[i]
        index <- max(which(dfq$fval < fval_add))
        if (index < nrow(dfq))
          value_new_add <- approx(dfq$fval[index : (index+1)], dfq$value_new[index : (index+1)], fval_add)$y
        dfq <- rbind(data.frame(dfq), data.frame(value = (value_new_add-1/2), fval = fval_add, cat = cat_add, value_new = value_new_add))
        dfq <- dfq[with(dfq, order(fval)),]}
    }

    for (i in 1:nrow(dfq)){
      if (dfp$fval[1] == 0 & dfq$fval[i] < 1 & dfq$fval[i] %notin% dfp$fval){
        cat_add <- sample1
        fval_add <- dfq$fval[i]
        index <- max(which(dfp$fval < fval_add))
        if (index < nrow(dfp))
          value_new_add <- approx(dfp$fval[index : (index+1)], dfp$value_new[index : (index+1)], fval_add)$y
        dfp <- rbind(data.frame(dfp), data.frame(value = (value_new_add-1/2), fval = fval_add, cat = cat_add, value_new = value_new_add))
        dfp <- dfp[with(dfp, order(fval)),]}
    }

    colnames(dfp) <- c("value", "fval", "cat", "value_new_p")
    dfp1 <- dfp[,c("fval", "value_new_p")]
    colnames(dfq) <- c("value", "fval", "cat", "value_new_q")
    dfq1 <- dfq[,c("fval", "value_new_q")]
    df_pq <- merge(dfp1, dfq1, by = "fval")
    #df_tq <- df_pq
    return(df_pq)
  }

  dfp <- new_quantile(P, sample1)
  dfq <- new_quantile(Q, sample2)

  df_tq <- qq_interpolation(dfp, dfq, sample1, sample2)

  # scaling: scale the quantiles by taking the sample mean of the theoretical quantile for computing QQ ABC
  f <- mean(df_tq$value_new_q)

  p <- ggplot(data = df_tq) +
    geom_line(data = df_tq, aes(x=value_new_q, y=value_new_p), color="red") +
    geom_point(data = df_tq, aes(x=value_new_q, y=value_new_p)) +
    geom_segment(data = df_tq, aes(x = value_new_q,
                                   y = value_new_q,
                                   xend = value_new_q,
                                   yend = value_new_p),
                 color="black") +
    geom_abline(intercept = 0, slope = 1, color = "goldenrod", size=1) +
    xlab(sample2) + ylab(sample1) +
    coord_cartesian(xlim =c(-1, max(df_tq$value_new_q + 1)),
                    ylim = c(-1, max(df_tq$value_new_p + 1))) +
    ggtitle(paste("Q-Q Plot, scaled ABC: ",
                  round(calcQQabc(P, Q, sample1="Sample quantiles",
                                  sample2="Theoretical quantiles"), 3), sep="") ) +
    theme_GOF(base_size=14, angle=0, strip.text_size=12)
  print(p)
}



