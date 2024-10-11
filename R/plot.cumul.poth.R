#' Plot cumulative precision of treatment hierarchy (POTH) values
#'
#' @description
#' Plot cumulative precision of treatment hierarchy (POTH) values
#'
#' @param x R object of class \code{cumul.poth}.
#' @param labels A logical indicating whether treatment names should be
#'   shown in the plot.
#' @param nchar.trts Number of characters to keep for each treatment name if labels = TRUE.
#' @param digits Minimal number of significant digits for cumulative POTH, see
#'   \code{\link{print.default}}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Plot results of cumulative method for precision of treatment hierarchy (POTH) metric
#' (Wigle et al., 2024).
#'
#' @return
#' A ggplot2 object.
#'
#' @author Augustine Wigle \email{amhwigle@@uwaterloo.ca},
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
#' data(Senn2013)
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE)
#'
#' # Cumulative method
#' c1 <- cumul(poth(net1))
#' c1
#' plot(c1)
#' plot(c1, labels = TRUE)
#'
#' c2 <- cumul(poth(net1), sort = FALSE)
#' c2
#' plot(c2)
#' plot(c2, labels = TRUE)
#' }
#'
#' @method plot cumul.poth
#' @export

plot.cumul.poth <- function(x, labels = FALSE, nchar.trts = 4, digits = 3, ...) {

  chkclass(x, "cumul.poth")
  chknumeric(nchar.trts, min = 1, length = 1)

  poth_cum <- x$poth_cum[-1]
  df <- data.frame(poth = poth_cum, id = seq_along(poth_cum) + 1,
                   trt = x$trt[-1])
  #
  if (labels) {
    trts.abbr <- treats(x$trt, nchar.trts)
    df$labels <- paste("+", trts.abbr[-1])
    df$labels[1] <- paste(trts.abbr[1], df$labels[1])
    df$labels <- factor(df$labels, levels = df$labels)
    #
    xlab <- "Treatments"
  }
  else {
    df$labels <- df$id
    #
    if (attr(x, "sort"))
      xlab <- "Best k Treatments"
    else
      xlab <- "Number of treatments"
  }

  p <- ggplot(df, aes(x = labels, y = poth)) +
    geom_col(col = "black") +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = xlab, y = "cPOTH")
  #
  if (!labels)
    p <- p + scale_x_continuous(breaks = 2:max(df$id), minor_breaks = NULL)

  p
}
