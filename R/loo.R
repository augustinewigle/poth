#' Leave-one-out method for precision of treatment hierarchy (POTH) metric
#'
#' @param x An R object of class \code{poth}.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param \dots Additional arguments.
#'
#' @return A data frame with additional class \code{loo.poth} and the following
#'   variables:
#' \item{trt}{Treatment names.}
#' \item{score}{Ranking metric (global).}
#' \item{rank}{Treatment rank (global).}
#' \item{residual}{Residual (global POTH minus leave-one-out POTH.}
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
#' loo1 <- loo(poth(net1))
#' loo1
#'
#' data(Senn2013)
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", random = FALSE)
#'
#' }
#'
#' @rdname loo
#' @method loo poth
#' @export

loo.poth <- function(x, ...) {

  chkclass(x, "poth")

  n <- x$n
  trts <- x$trts

  if (x$input == "mcmc.samples") {
    score_type <- "SUCRA"

    ranking <- x$ranking
    samples <- x$x
    colnames(samples) <- trts
    small.values <- x$small.values
    #
    # if (sort)
      seq <- order(ranking, decreasing = TRUE)
    # else
    #   seq <- seq_along(ranking)
    #
    loo_rps <-
      lapply(seq_len(n),
             function(x)
               rankMCMC(samples[, -x], small.values))
    #
    loo_poths <- sapply(loo_rps, function(x) poth(x)$poth)
    names(loo_poths) <- colnames(samples)
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
    # if (sort)
      seq <- order(ranking, decreasing = TRUE)
    # else
    #   seq <- seq_along(ranking)
    #
    loo_pscores <-
      lapply(seq_len(n),
             function(drp)
               pscores(TE[-drp, -drp], seTE[-drp, -drp], small.values))
    names(loo_pscores) <- colnames(TE)
    #
    loo_poths <- sapply(loo_pscores, function(x) poth(x)$poth)
    names(loo_poths) <- names(loo_pscores)
  }
  else
    stop("Leave-one-out method not available for input type '", x$input, "'.")

  # Put everything together
  #
  residuals <- x$poth - loo_poths
  #

  res <- data.frame(trt = names(residuals),
                    rank = rank(-ranking),
                    score = ranking,
                    residuals = residuals,
                    ratio = residuals / sum(abs(residuals)))[seq, ]
  #
  attr(res, "poth") <- x$poth
  attr(res, "score_type") <- score_type
  #
  class(res) <- c("loo.poth", class(res))
  #
  res
}


#' @rdname loo
#' @export loo

loo <- function(x, ...)
  UseMethod("loo")


#' @rdname loo
#' @keywords print
#' @method print loo.poth
#' @export

print.loo.poth <- function(x, digits = 3, ...) {

  chkclass(x, "loo.poth")
  #
  poth <- attr(x, "poth")
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
