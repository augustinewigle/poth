#' Plot results of leave-one-out method
#' 
#' @description
#' Plot results of leave-one-out method for separation in ranking (SIR) metric
#' 
#' @param x R object of class \code{sir}.
#' @param plot_trts A logical indicating whether treatment names be included
#'   in the plot.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Plot results of leave-one-out method for separation in ranking (SIR) metric
#' (Wiggle et al., 2024).
#'
#' @return
#' A ggplot2 object.
#'
#' @author Augustine Wigle \email{amhwigle@uwaterloo.ca},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @references
#' Wiggle A, ... (2024):
#' Separation In Ranking: A Metric for Quantifying Uncertainty in Treatment
#' Hierarchies in Network Meta-Analysis
#'
#' @examples
#' \dontrun{
#' # R package netmeta must be available
#' library("netmeta")
#' data(smokingcessation)
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1, random = FALSE)
#' 
#' # Use P-scores to calculate SIR
#' loo1 <- loo(sir(net1))
#' loo1
#' plot(loo1)
#' }
#' 
#' @method plot loo.sir
#' @export

plot.loo.sir <- function(x, plot_trts = TRUE, ...) {
  
  chkclass(x, "loo.sir")
  
  res_pos <- x$residuals > 0
  
  # Get rid of warning "no visible binding for global variable"
  
  n <- residuals <- score_type <- plot_trts <- NULL
  n.seq <- seq_len(x$n)
  
  g <- ggplot(x,
              aes(x = n.seq, y = residuals,
                  fill = res_pos, alpha = abs(residuals))) +
    geom_col(col = "black") +
    coord_cartesian(ylim = c(-max(abs(x$residuals)), max(abs(x$residuals)))) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(breaks = c(F, T), values = c("red4", "darkgreen"),
                      labels = c("Reduces certainty", "Increases certainty")) +
    labs(y = expression(SIR-SIR^"-i"),
         x = paste0("Global Rank Using ", score_type, "s"),
         fill = str_wrap(paste0("Contribution to SIR\n(Global SIR ",
                                ifelse(sir < 0.001, "< 0.001",
                                       paste0("= ", round(sir, digits = 3))),
                                ")"), width = 20)) +
    guides(alpha = "none") +
    theme_bw()
  
  if (plot_trts) {
    g <- g + scale_x_continuous(
      sec.axis = sec_axis(~.,
                          breaks = n.seq,
                          labels = str_trunc(names(x$residuals), width = 13),
                          name = "Treatment Name"),
      breaks = n.seq, minor_breaks = NULL) +
      theme(axis.text.x.top =
              element_text(angle = 90, hjust = 0.01, vjust = 0.4),
            axis.title.x.top = element_blank())
  }
  else
    g <- g + scale_x_continuous(breaks = n.seq, minor_breaks = NULL)
  
  g  
}
