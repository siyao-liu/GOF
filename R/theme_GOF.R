#' GOF theme ggplots
#' @description This function generates ggplot object with theme elements that are used in the plots used in the GOF package.
#'
#' @param base_size Base font size in pts.
#' @param angle Angle (in [0, 45, 90]) for the axis label.
#' @param strip.text_size Text size in pts.
#' @return theme_GOF returns a list that can be added to a ggplot object.
#' @export
theme_GOF <- function(base_size=14, angle=0, strip.text_size=12) {
  if (angle==90) {
    theme_bw(base_size=base_size) +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(color="black"),
        panel.background = element_rect(color="black"),
        strip.text = element_text(size=strip.text_size ),
        strip.background = element_rect(fill="white"))
  } else if (angle==45) {
    theme_bw(base_size=base_size) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0),
            axis.text = element_text(color="black"),
            panel.background = element_rect(color="black"),
            strip.text = element_text(size=strip.text_size ),
            strip.background = element_rect(fill="white"))
  } else {
    theme_bw(base_size=base_size) +
      theme(
        axis.text = element_text(color="black"),
        panel.background = element_rect(color="black"),
        strip.text = element_text(size=strip.text_size ),
        strip.background = element_rect(fill="white"))

  }

}
