#' P-P plot comparing two probabilities.
#'
#' @description A function for generating a P-P plot to compare whether two samples come from the same distribution.
#'
#' @param P A numeric vector of the empirical probabilities.
#' @param Q A numeric vector of the theoretical probabilities.
#' @param prob1 A character to denote the first probability distribution.
#' @param prob2 A character to denote the first probaility distribution.
#' @return ppplot returns a P-P plot using ggplot().
#' @export
ppplot <- function(P, # empirical probabilities
                   Q, # theoretical probabilities
                   prob1,
                   prob2) {
  require(ggplot2)
  df <- data.frame("x"=P, "y"=Q)
  tog <- rbind(c(0,0), df, c(1,1)) # add point(0,0) & point(1,1)

  p <- ggplot(data=tog, aes(y, x)) +
    geom_point(color = "black", shape=1) +
    xlab(prob2) + ylab(prob1) +
    geom_line(data=tog, aes(y, x), color="red", size=1) +
    geom_abline(intercept=0, slope=1, color="goldenrod", size=1) +
    geom_point(x=0,y=0, color="black") +
    geom_point(x=1,y=1, color="black") +
    geom_segment(aes(x = y, y = y, xend = y, yend = x), color="black") +
    ggtitle(paste("P-P Plot, ABC: ", round(calcPPabc(tog$y, tog$x), 3), sep="") ) +
    scale_x_continuous(limits=c(0, 1)) + scale_y_continuous(limits=c(0, 1)) +
    theme_GOF(base_size=14, angle=0, strip.text_size=12)
  print(p)
}
