#' Compute the area between curves (abc) in the Q-Q plot
#' @description This function computes the area between curves in the Q-Q plot.
#'              It takes the empirical and theoretical quantiles as input and calculates the area between the 45 degree line and the the quantile of the theoretical distribution against the quantile of observed distribution.
#'
#' @param P A numeric vector of sampled data points from the first distribution.
#' @param Q A numeric vector of sampled data points from the second distribution.
#' @param sample1 A character to denote the first distribution.
#' @param sample2 A character to denote the second distribution.
#' @return calcQQabc returns a vector of the area between curves (abc) values for each gene using the Q-Q method.
#' @export
calcQQabc <- function(P, Q, sample1, sample2) {
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

  tog <- data.frame("x"=df_tq$value_new_q/f, "y"=df_tq$value_new_p/f)

  # the height (vertical), which is the difference between x and y
  v <- abs(tog$y - tog$x)

  # delta x which is the difference between xi
  h <- diff(tog$x)

  if ( all(tog[-c(1, nrow(tog)), ]$x < tog[-c(1, nrow(tog)), ]$y) |
       all(tog[-c(1, nrow(tog)), ]$x > tog[-c(1, nrow(tog)), ]$y) )  {
    # the red curve is on one side of the 45 degree line

    area <- c()
    # trapezoid calculation
    for (i in 1: length(h)) {
      area[i] <- (v[i] + v[i+1])*h[i]/2
    }
    # sum up all the trapezoids
    abc <- sum(area)
  }
  else { # the red curve crosses with the 45 degree line

    # calculate the area of each trapezoid
    area <- c()
    # trapezoid calculation
    for (i in 1: length(h)) {
      area[i] <- (v[i] + v[i+1])*h[i]/2
    }

    # calculate the difference between x and y
    delta_xy <- tog$x - tog$y

    if ( all(delta_xy==0) ) {
      abc <- 0

    } else {
      # find the first sign change
      pos <- which(diff(sign(delta_xy))!=0)

      xc_area <- c()
      for (i in 1: length(pos)) {
        # find indices of the neighboring points on the red ROC curve that intersect the 45 degree line
        ind1 <- pos[i]
        ind2 <- pos[i] + 1
        # get the coordinates of these two points
        pt1 <- tog[ind1, ]
        pt2 <- tog[ind2, ]

        # calculate m
        m <- (pt2$y - pt1$y) / (pt2$x - pt1$x)

        # then calculate Xc
        xc <- (pt1$y - m*pt1$x) / (1-m)

        # calculate the triangles around Xc (crossing where the 45 degree line and ROC curve crosses)
        xc_area[i] <- 0.5*(abs(pt1$y-pt1$x)*abs(xc-pt1$x)) + 0.5*(abs(pt2$y-pt2$x)*abs(pt2$x-xc))
      }

      # sum up the crossing triangles
      xc_area_sum <- sum(xc_area)

      # calculate the shaded area by summing up each trapezoids and substract the crossing trapezoids and
      # then add back the crossing triangles
      abc <- sum(area[-pos]) + xc_area_sum
    }
  }

  return(abc)
}
