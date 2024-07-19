#' Calculate the local SIR for a subset of treatments
#'
#' @param x An object of class \code{sir}.
#' @param subset A character vector of treatment names to consider as the set
#'   of competing treatments.
#' @param top A single integer to define the number of treatments with the
#'   largest ranking metric to consider in subset.
#' @param bottom A single integer to define the number of treatments with the
#'   smallest ranking metric to consider in subset.
#' @param \dots Additional arguments (ignored).
#' 
#' @return An R object of class \code{sir}.
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
#' # Use P-scores to calculate local SIR for treatments "A" and "C"
#' subset(sir(net1), subset = c("A", "C"))
#' 
#' # Use P-scores to calculate local SIR for first three treatments
#' subset(sir(net1), top = 3)
#' 
#' # Use P-scores to calculate local SIR for first three treatments
#' subset(sir(net1), bottom = 3)
#' }
#' 
#' @method subset sir
#' @export

subset.sir <- function(x, subset, top, bottom, ...) {
  
  chkclass(x, "sir")
  
  if (x$input == "mcmc.samples") {
    score_type <- "SUCRA"
    scores <- x$ranking
    samples <- x$x
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
    stop("Local SIR for a subset of treatments not available for input type '",
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
  
  if (x$input == "mcmc.samples")
    return(sir(samples[, seq, drop = FALSE], small.values = small.values))
  else if (x$input %in% c("effects.se", "netmeta"))
    return(sir(pscores(TE[seq, seq, drop = FALSE],
                       seTE[seq, seq, drop = FALSE],
                       small.values)))
}
