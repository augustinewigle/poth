#' Calculate the local POTH for a subset of treatments
#'
#' @param x An object of class \code{poth}.
#' @param subset A character vector of treatment names to consider as the set
#'   of competing treatments.
#' @param top A single integer to define the number of treatments with the
#'   largest ranking metric to consider in subset.
#' @param bottom A single integer to define the number of treatments with the
#'   smallest ranking metric to consider in subset.
#' @param \dots Additional arguments (ignored).
#'
#' @return An R object of class \code{poth}.
#'
#' @examples
#' library("netmeta")
#' data(smokingcessation)
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'   event = list(event1, event2, event3), n = list(n1, n2, n3),
#'   data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1, random = FALSE)
#'
#' # Use P-scores to calculate local POTH for treatments "A" and "C"
#' subset(poth(net1), subset = c("A", "C"))
#'
#' # Use P-scores to calculate local POTH for first three treatments
#' subset(poth(net1), top = 3)
#'
#' # Use P-scores to calculate local POTH for first three treatments
#' subset(poth(net1), bottom = 3)
#'
#' @method subset poth
#' @export

subset.poth <- function(x, subset, top, bottom, ...) {

  chkclass(x, "poth")

  if (x$input == "mcmc.samples") {
    score_type <- "SUCRA"
    #
    scores <- x$ranking
    #
    samples <- x$x
    #
    trts <- x$trts
    small.values <- x$small.values
    n <- ncol(samples)
  }
  else if (x$input %in% c("effects.se", "netmeta")) {
    score_type <- "P-score"
    #
    scores <- x$ranking
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
    trts <- colnames(TE)
    small.values <- x$small.values
    n <- ncol(TE)
  }
  else
    stop("Local POTH for a subset of treatments not available for input type '",
         x$input, "'.")

  if (as.numeric(!missing(subset)) +
      as.numeric(!missing(top) | !missing(bottom)) != 1)
    stop("Please provide either argument 'subset' ",
         "or argument(s) 'top' and / or 'bottom'.")
  #
  trts <- x$trts
  if (!missing(subset))
    seq <- setchar(subset, trts)
  #
  if (!missing(top)) {
    chknumeric(top, min = 1, max = n - 1, integer = TRUE)
    #
    seq.top <- rank(scores, ties.method = "first") > n - top
    #
    if (missing(bottom))
      seq <- seq.top
  }
  #
  if (!missing(bottom)) {
    chknumeric(bottom, min = 1, max = n - 1, integer = TRUE)
    #
    seq.bottom <- rank(scores, ties.method = "first") <= bottom
    #
    if (missing(top))
      seq <- seq.bottom
  }
  #
  if (!missing(top) & !missing(bottom))
    seq <- seq.top | seq.bottom

  if (x$input == "mcmc.samples") {
    if (x$pooled != "")
      return(poth(samples[, seq, drop = FALSE], small.values = small.values,
                  pooled = x$pooled))
    else
      return(poth(samples[, seq, drop = FALSE], small.values = small.values))
  }
  else if (x$input %in% c("effects.se", "netmeta")) {
    if (x$pooled != "")
      return(poth(TE[seq, seq, drop = FALSE],
                  seTE[seq, seq, drop = FALSE],
                  small.values = small.values,
                  pooled = x$pooled))
    else
      return(poth(TE[seq, seq, drop = FALSE],
                  seTE[seq, seq, drop = FALSE],
                  small.values = small.values))
  }
}
