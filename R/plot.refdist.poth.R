#' Plot reference distribution for given network
#'
#' @description
#' Plot a histogram of simulated reference distribution for POTH
#'
#' @param x R object of class \code{refdist}.
#' @param observed A logical indicating whether a line with the observed POTH
#'   should be drawn in the plot.
#' @param probability A logical indicating whether the probability of a POTH
#'   greater than the observed POTH should be printed in the legend. This
#'   argument is ignored if \code{observed = FALSE}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Plots a histogram of the simulated POTH values from the reference
#' distribution. (Wigle et al., 2024).
#'
#' @return
#' A ggplot2 object.
#'
#' @author Augustine Wigle \email{amhwigle@uwaterloo.ca},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @references
#' Wigle A, Béliveau A, Salanti G, Rücker G, Schwarzer G, Mavridis D,
#' Nikolakopoulou A (2024):
#' Precision of treatment hierarchy: A metric for quantifying uncertainty in
#' treatment hierarchies in network meta-analysis.
#' Preprint on arXiv, \url{https://doi.org/10.48550/arXiv.2501.11596}
#'
#' @examples
#' library("netmeta")
#' data(smokingcessation)
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1, random = FALSE)
#'
#' # Reference distribution
#' dist1 <- refdist(net1, pooled = "common", nsim = 2)
#' plot(dist1, observed = FALSE)
#' 
#' \donttest{
#' data(Senn2013)
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", random = FALSE)
#'
#' # Reference distribution
#' dist2 <- refdist(net2, pooled = "common", nsim = 25)
#' plot(dist2, probability = TRUE)
#' }
#' 
#' @method plot refdist
#' @export

plot.refdist <- function(x, observed = TRUE, probability = FALSE, ...) {
  
  chkclass(x, "refdist")
  
  # Get rid of warning "no visible binding for global variable"
  #
  dist <- obs <- name2 <- NULL
  
  dat <- data.frame(dist = x$dist)
  
  if (observed) {
    p <- poth(x = x$netmeta, pooled = x$pooled)
    #
    if (probability)
      dat2 <-
        data.frame(obs = c(p$poth, 0),
                   name2 = c(paste0("Observed\nPOTH = ", round(p$poth, 3)),
                             paste0("Pr(POTH > ", round(p$poth, 3), ") = ",
                                    round(mean(x$dist > p$poth), 3))))
    
    else
      dat2 <-
        data.frame(obs = p$poth,
                   name2 = paste0("Observed\nPOTH = ", round(p$poth, 3)))
  }
  
  g <- ggplot(data = dat,
              aes(x = dist)) +
    geom_histogram(col = "black", fill = "lightgrey",
                   binwidth = (max(x$dist) - min(x$dist)) / 30) +
    theme_bw() +
    # scale_x_continuous(limits = c(0, 1)) +
    coord_cartesian(xlim = c(0, 1)) +
    geom_hline(yintercept = 0) +
    labs(x = "POTH", y = "Frequency", color = "")
  #
  if (observed)
    g <- g +
    geom_vline(data = dat2, aes(xintercept = obs, color = name2),
               show.legend = TRUE, linetype = "dashed", linewidth = 1) +
    scale_color_manual(values = c("cornflowerblue", "transparent")) +
    labs(x = "POTH", y = "Frequency", color = "") +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  #
  g
}
