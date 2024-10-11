#' Plot results of leave-one-out method
#'
#' @description
#' Plot results of leave-one-out method for precision of treatment hierarchy (POTH) metric
#'
#' @param x R object of class \code{poth}.
#' @param labels A logical indicating whether treatment names should be
#'   shown in the plot.
#' @param digits Minimal number of significant digits for global POTH, see
#'   \code{\link{print.default}}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Plot results of leave-one-out method for precision of treatment hierarchy (POTH) metric
#' (Wigle et al., 2024).
#'
#' @return
#' A ggplot2 object.
#'
#' @author Augustine Wigle \email{amhwigle@uwaterloo.ca},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @references
#' Wigle, A., Béliveau, A., Salanti, G., Rücker, G., Schwarzer, G., Mavridis, D.,
#' Nikolakopoulou, A. (2024):
#' Precision of Treatment Hierarchy: A Metric for Quantifying Uncertainty in Treatment
#' Hierarchies in Network Meta-Analysis
#'
#' @examples
#' \dontrun{
#' library("netmeta")
#' data(smokingcessation)
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1, random = FALSE)
#'
#' # Leave-one-out method
#' loo1 <- loo(poth(net1))
#' loo1
#' plot(loo1)
#'
#' data(Senn2013)
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", random = FALSE)
#'
#' # Leave-one-out method (without sorting by ranking metric)
#' loo2 <- loo(poth(net2), sort = FALSE)
#' loo2
#' plot(loo2)
#' }
#'
#' @method plot loo.poth
#' @export

plot.loo.poth <- function(x, labels = TRUE, digits = 3, ...) {

  chkclass(x, "loo.poth")

  res_pos <- x$resid > 0
  poth <- attr(x, "poth")
  score_type <- attr(x, "score_type")

  # Get rid of warning "no visible binding for global variable"

  n <- resid <- score_type <- NULL
  n.seq <- seq_len(nrow(x))

  g <- ggplot(x,
              aes(x = n.seq, y = resid,
                  fill = res_pos, alpha = abs(resid))) +
    geom_col(col = "black") +
    coord_cartesian(ylim = c(-max(abs(x$resid)), max(abs(x$resid)))) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(breaks = c(FALSE, TRUE), values = c("red4", "darkgreen"),
                      labels = c("Reduces certainty", "Increases certainty")) +
    labs(y = expression(r["POTH"]~"(j)"),
         x = paste0("Place in Hierarchy"),
         fill = str_wrap(paste0("Contribution to POTH (",
                                formatPT(poth, digits = digits, lab = TRUE,
                                         labval = "POTH"),
                                ")"), width = 20)) +
    guides(alpha = "none") +
    theme_bw()

  if (labels) {
    g <- g + scale_x_continuous(
      sec.axis = sec_axis(~.,
                          breaks = n.seq,
                          labels = str_trunc(x$trt, width = 13),
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
