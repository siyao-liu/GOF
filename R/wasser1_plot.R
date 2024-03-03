#' Make diagnostic plot for discrete 1-Wasserstein distance method
#' @description This function generates diagnostic plots for assessing the 1-Wasserstein distance method.
#'
#' @param vigrid A vector of grid from 0,1,2,..., maximum of data point.
#' @param vpi A numeric vector of the empirical probabilities.
#' @param vqi A numeric vector of the theoretical probabilities.
#' @param vDj A vector of the difference between the two quantile functions.
#' @param Wdist A string of the 1-Wasserstein distance metric.
#' @param nmax The maximum of the biggest count in vp_i and the biggest count in vq_i.
#' @param GeneName The name of the gene.
#' @return wasser1_plot returns a diagnostic plot using ggplot().
#'         Left panel is a bar plot of the empirical gene count distribution and the theroretical gene count distribution.
#'         Right panel is a bar plot of the difference between the two distributions.
#' @export
wasser1_plot <- function(vigrid, vpi, vqi, vDj, Wdist, nmax, GeneName) {
  require(ggplot2)
  require(cowplot)
  plot.df <- data.frame("vigrid"=vigrid,
                        "vpi"=vpi,
                        "vqi"=vqi)
  # make probability distributions
  p1 <- ggplot(plot.df) +
    # first plot observed bars - black
    geom_bar(stat="identity", aes(x = vigrid, y=vpi, fill="vpi")) +
    # then plot expected bars - green
    geom_bar(stat="identity", aes(x = vigrid, y=vqi, fill="vqi"), width=0.6, alpha=0.6) +
    ggtitle(GeneName) +
    scale_fill_manual("", values=c("vpi"="black", "vqi"="green4"),
                      labels=c("Observed", "Expected")) +
    xlab("Gene count")  + ylab("Probability") +
    theme_GOF(base_size=14, angle=0, strip.text_size=12) +
    theme(legend.position = c(0.75, 0.85),
          legend.background = element_rect(fill = NA, color = NA, size=1))

  # plot vDj
  plot.df2 <- data.frame("index"=seq(0, nmax, 1), "D"=c(0, vDj))
  plot.df2$sign <- as.character(sign(plot.df2$D))
  p2 <- ggplot(plot.df2, aes(x=index, y=D, fill=sign)) +
    geom_bar(stat="identity", width=0.6) +
    scale_fill_manual("", values=c("-1"="darkblue", "1"="firebrick")) +
    geom_hline(yintercept=0) +
    xlab("Gene count") +
    ylab("D") +
    theme_GOF(base_size=14, angle=0, strip.text_size=12) +
    theme(legend.position = "none") +
    ggtitle(paste0(GeneName, ", Wdist=", round(Wdist, 2)))

  # put two plots together
  cowplot::plot_grid(p1, p2, ncol=1)

}
