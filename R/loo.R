#' Leave-one-out method for separation in ranking (SIR) metric
#'
#' @param x An R object of class \code{sir}.
#' @param sort A logical indicating whether results should be sorted
#'   by decreasing ranking metric.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param \dots Additional arguments.
#'
#' @return A data frame with additional class \code{loo.sir} and the following
#'   variables:
#' \item{trt}{Treatment names.}
#' \item{score}{Ranking metric (global).}
#' \item{rank}{Treatment rank (global).}
#' \item{residual}{Residual (global SIR minus leave-one-out SIR.}
#' \item{ratio}{Ratio of residual devided by absolute sum of residuals.}
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
#' # Leave-one-out method
#' loo1 <- loo(sir(net1))
#' loo1
#'
#' data(Senn2013)
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", random = FALSE)
#'
#' # Leave-one-out method (without sorting by ranking metric)
#' loo2 <- loo(sir(net2), sort = FALSE)
#' loo2
#' }
#'
#' @rdname loo
#' @method loo sir
#' @export

loo.sir <- function(x, sort = TRUE, ...) {

  chkclass(x, "sir")

  n <- x$n
  trts <- x$trts

  if (x$input == "mcmc.samples") {
    score_type <- "SUCRA"

    ranking <- x$ranking
    samples <- x$x
    colnames(samples) <- trts
    small.values <- x$small.values
    #
    if (sort)
      seq <- order(ranking, decreasing = FALSE)
    else
      seq <- seq_along(ranking)
    #
    loo_rps <-
      lapply(seq_len(n),
             function(x)
               rankMCMC(samples[, -x], small.values))
    #
    loo_sirs <- sapply(loo_rps, function(x) sir(x)$sir)
    names(loo_sirs) <- colnames(samples)
  }
  else if (x$input %in% c("effects.se", "netmeta")) {
    score_type <- "P-score"
    #
    ranking <- x$ranking
    #
    if (x$input == "effects.se") {
      TE <- x$x
      seTE <- x$se
    }
    else {
      TE <- x$TE
      seTE <- x$se
    }
    #
    small.values <- x$small.values
    #
    if (sort)
      seq <- order(ranking, decreasing = TRUE)
    else
      seq <- seq_along(ranking)
    #
    loo_pscores <-
      lapply(seq_len(n),
             function(drp)
               pscores(TE[-drp, -drp], seTE[-drp, -drp], small.values))
    names(loo_pscores) <- colnames(TE)
    #
    loo_sirs <- sapply(loo_pscores, function(x) sir(x)$sir)
    names(loo_sirs) <- names(loo_pscores)
  }
  else
    stop("Leave-one-out method not available for input type '", x$input, "'.")

  # Put everything together
  #
  residuals <- x$sir - loo_sirs
  #

  res <- data.frame(trt = names(residuals),
                    rank = rank(-ranking),
                    score = ranking,
                    residuals = residuals,
                    ratio = residuals / sum(abs(residuals)))[seq, ]
  #
  attr(res, "sir") <- x$sir
  attr(res, "score_type") <- score_type
  #
  class(res) <- c("loo.sir", class(res))
  #
  res
}


#' @rdname loo
#' @export loo

loo <- function(x, ...)
  UseMethod("loo")


#' @rdname loo
#' @keywords print
#' @method print loo.sir
#' @export

print.loo.sir <- function(x, digits = 3, ...) {

  chkclass(x, "loo.sir")
  #
  sir <- attr(x, "sir")
  score_type <- attr(x, "score_type")
  #
  rownames(x) <- x$trt
  x$trt <- NULL
  #
  x$score <- round(x$score, digits)
  x$residuals <- round(x$residuals, digits)
  x$ratio <- round(x$ratio, digits)
  #
  class(x) <- "data.frame"
  #
  print(x)
  #
  invisible(NULL)
}
