#' Leave-one-out method for separation in ranking (SIR) metric
#'
#' @param x An R object of class \code{sir}.
#' @param sort A logical indicating whether results should be sorted
#'   by decreasing ranking metric.
#' @param \dots Additional arguments.
#' 
#' @return A data frame with the following variables:
#' \item{residuals}{...}
#' \item{sir}{...}
#' \item{trt_name}{...}
#' \item{orig_score}{...}
#' \item{score_type}{...}
#' \item{rank}{...}
#' \item{ratio}{...}
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
    scores <- x$ranking
    samples <- x$x
    small.values <- x$small.values
    #
    if (sort) {
      seq <- order(scores, decreasing = FALSE)
      samples <- samples[ , seq]
    }
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
    small.values <- x$small.values
    #
    if (sort) {
      seq <- order(scores, decreasing = TRUE)
      TE <- TE[seq, seq]
      seTE <- seTE[seq, seq]
    }
    #
    loo_pscores <-
      lapply(seq_len(n),
             function(x)
               pscores(TE[-x, -x], seTE[-x, -x], small.values))
    names(loo_pscores) <- colnames(TE)
    #
    loo_sirs <- sapply(loo_pscores, function(x) sir(x)$sir)
    names(loo_sirs) <- names(loo_pscores)
  }
  else
    stop("Leave-one-out method not available for input type '", x$input, "'.")
  
  #
  # Put everything together
  #
  sir_res <- x$sir - loo_sirs
  #
  res <- data.frame(residuals = sir_res,
                    sir = x$sir,
                    #
                    trt_name = names(sir_res),
                    orig_score = scores[seq],
                    score_type = score_type,
                    rank = rank(-scores[seq]),
                    ratio = sir_res / sum(abs(sir_res)))
  #
  if (sort)
    res <- res[order(res$rank), ]
  
  res
}


#' @rdname loo
#' @export loo

loo <- function(x, ...)
  UseMethod("loo")
