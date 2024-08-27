#' Plot results of leave-one-out method
#'
#' @description
#' Plot results of cumulative method for separation in ranking (SIR) metric
#'
#' @param x R object of class \code{sir}.
#' @param labels A logical indicating whether treatment names should be
#'   shown in the plot.
#' @param trt.trunc Number of characters to keep for each treatment name if labels = T
#' @param digits Minimal number of significant digits for global SIR, see
#'   \code{\link{print.default}}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' Plot results of cumulative method for separation in ranking (SIR) metric
#' (Wigle et al., 2024).
#'
#' @return
#' A ggplot2 object.
#'
#' @author Augustine Wigle \email{amhwigle@@uwaterloo.ca},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @references
#' Wigle A, ... (2024):
#' Separation In Ranking: A Metric for Quantifying Uncertainty in Treatment
#' Hierarchies in Network Meta-Analysis
#'
#' @method plot cumul.sir
#' @export

plot.cumul.sir <- function(x, labels = FALSE, trt.trunc = 3, digits = 3, ...) {

  chkclass(x, "cumul.sir")

  df <- data.frame(sir = x$cumul.sir,
                   grp = 2:(length(x$groups)+1))

  if(labels) {

    xlab <- "Treatments"

    if(trt.trunc > 0) {

      short <- str_trunc(names(x$ranking), width = trt.trunc, ellipsis = "")
      ordered_short <- short[order(x$ranking, decreasing = TRUE)]
      labs <- sapply(2:(length(df$sir)+1), function(x) paste(ordered_short[1:x], collapse = ", "))
      df$labels = str_wrap(labs, width = 10)


    } else {

      ordered <- names(x$ranking)[order(x$ranking, decreasing = TRUE)]
      labs <- sapply(2:(length(df$sir)+1), function(x) paste(ordered[1:x], collapse = ", "))
      df$labels = str_wrap(labs, width = 10)


    }

    p <- ggplot(df, aes(x = labels, y = sir)) +
      geom_col(col = "black") +
      geom_hline(yintercept = 0) +
      scale_y_continuous(limits = c(0,1)) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = xlab, y = "cPOTH")

  } else {

    xlab <- "Best k Treatments"

    df$labels <- df$grp

    p <- ggplot(df, aes(x = labels, y = sir)) +
      geom_col(col = "black") +
      geom_hline(yintercept = 0) +
      scale_y_continuous(limits = c(0,1)) +
      scale_x_continuous(breaks = 2:max(df$grp), minor_breaks = NULL) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = xlab, y = "cPOTH")

  }



  return(p)

}
