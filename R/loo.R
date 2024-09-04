#' Leave-one-out method for precision of treatment hierarchy (POTH) metric
#'
#' @param x An R object of class \code{poth}.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments.
#'
#' @return A data frame with additional class \code{loo.poth} and the following
#'   variables:
#' \item{trt}{Treatment names.}
#' \item{rank}{Treatment rank (global).}
#' \item{score}{Ranking metric (global).}
#' \item{poth_loo}{Leave-one-out POTH.}
#' \item{resid}{Residuals (global POTH minus leave-one-out POTH.}
#' \item{ratio}{Ratio of residual devided by absolute sum of residuals.}
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
    poth_loo <- sapply(loo_rps, function(x) poth(x)$poth)
    names(poth_loo) <- colnames(samples)
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
    poth_loo <- sapply(loo_pscores, function(x) poth(x)$poth)
    names(poth_loo) <- names(loo_pscores)
  }
  else
    stop("Leave-one-out method not available for input type '", x$input, "'.")

  # Put everything together
  #
  resid <- x$poth - poth_loo
  #

  res <- data.frame(trt = names(resid),
                    rank = rank(-ranking),
                    score = ranking,
                    poth_loo = poth_loo,
                    resid = resid,
                    ratio = resid / sum(abs(resid)))[seq, ]
  #
  attr(res, "poth") <- x$poth
  attr(res, "score_type") <- score_type
  attr(res, "pooled") <- x$pooled
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

print.loo.poth <- function(x, digits = 3, legend = TRUE, ...) {

  chkclass(x, "loo.poth")
  #
  chknumeric(digits, min = 0, length = 1)
  chklogical(legend)
  #
  poth <- attr(x, "poth")
  score_type <- attr(x, "score_type")
  pooled <- attr(x, "pooled")
  
  txt <- "Leave-one-out method"
  #
  if (pooled != "")
    txt <- paste0(txt,
                  " (",
                  if (pooled == "common") "common" else "random",
                  " effects model)")
  #
  txt <- paste0(txt, "\n\n")
  #
  cat(txt)
  #
  cat(paste0("Precision of treatment hierarchy (global POTH) = ",
             round(poth, digits = digits), "\n\n"))
  
  rownames(x) <- x$trt
  x$trt <- NULL
  #
  x$score <- round(x$score, digits)
  x$poth_loo <- round(x$poth_loo, digits)
  x$resid <- round(x$resid, digits)
  x$ratio <- round(x$ratio, digits)
  #
  class(x) <- "data.frame"
  #
  print(x)
  
  if (legend) {
    cat("\nLegend:\n")
    cat(" rank     - Treatment rank (global)\n")
    cat(" score    - Ranking metric (global)\n")
    cat(" poth_loo - Leave-one-out POTH\n")
    cat(" resid    - Residual (global POTH minus leave-one-out POTH)\n")
    cat(" ratio    - Ratio of residual devided by absolute sum of residuals\n")
  }
  
  invisible(NULL)
}
