#' Compute the area between curves (abc) in the P-P plot
#' @description This function computes the area between curves in the P-P plot.
#'              It takes the empirical and theoretical probabilities as input and calculates the area between the 45 degree line and the CDF of the theoretical distribution against the CDF of the observed distribution.
#' @param P A numeric vector of the empirical probabilities.
#' @param Q A numeric vector of the theoretical probabilities.
#' @return calcPPabc returns a vector of the area between curves (abc) values for each gene using the P-P method.
#' @export
calcPPabc <- function(P, Q) {

  tog <- data.frame("x"=P, "y"=Q)
  tog <- rbind(c(0,0), tog, c(1,1)) # add point(0,0) & point(1,1)

  # the height (vertical), which is the difference between x and y
  v <- abs(tog$y - tog$x)

  # delta x which is the difference between xi
  h <- diff(tog$x)

  if ( all(tog[-c(1, nrow(tog)), ]$x < tog[-c(1, nrow(tog)), ]$y) |
       all(tog[-c(1, nrow(tog)), ]$x > tog[-c(1, nrow(tog)), ]$y) )  { # the red curve is on one side of the 45 degree line

    area <- c()
    # trapezoid calculation
    for (i in 1: length(h)) {
      area[i] <- (v[i] + v[i+1])*h[i]/2
    }
    # sum up all the trapezoids
    abc <- sum(area)

  } else { # the red curve crosses with the 45 degree line

    # calculate the area of each trapezoid
    area <- c()
    # trapezoid calculation
    for (i in 1: length(h)) {
      area[i] <- (v[i] + v[i+1])*h[i]/2
    }

    # calculate the difference between x and y
    delta_xy <- tog$x - tog$y

    # find the first sign change
    pos <- which(diff(sign(delta_xy))!=0)
    pos <- pos[-c(1, length(pos))] # remove the (0,0) point and (1,1) point

    xc_area <- c()
    for (i in 1: length(pos)) {
      # find indices of the neighboring points on the red ROC curve that intersect the 45 degree line
      ind1 <- pos[i]
      ind2 <- pos[i] + 1
      # get the coordinates of these two points
      pt1 <- tog[ind1, ]
      pt2 <- tog[ind2, ]

      if (pt2$x != pt1$x) {
        # calculate m
        m <- (pt2$y - pt1$y) / (pt2$x - pt1$x)
        # then calculate Xc
        xc <- (pt1$y - m*pt1$x) / (1-m)
      } else {
        xc <- pt1$x
      }

      # calculate the triangles around Xc (crossing where the 45 degree line and ROC curve crosses)
      xc_area[i] <- 0.5*(abs(pt1$y-pt1$x)*abs(xc-pt1$x)) + 0.5*(abs(pt2$y-pt2$x)*abs(pt2$x-xc))

    }

    # sum up the crossing triangles
    xc_area_sum <- sum(xc_area)

    # calculate the shaded area by summing up each trapezoids and substract the crossing trapezoids and
    # then add back the crossing triangles
    abc <- sum(area[-pos]) + xc_area_sum
  }
  return(abc)
}
