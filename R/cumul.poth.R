#' Best k treatments method for precision of treatment hierarchy (POTH)
#' metric
#'
#' @param x An R object of class \code{poth}.
#' @param sort A logical indicating whether results should be sorted
#'   by decreasing ranking metric.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param \dots Additional arguments.
#'
#' @return A data frame with additional class \code{bestk.poth} and the
#'   following variables:
#' \item{trt}{Name of added treatment.}
#' \item{rank}{Treatment rank (global).}
#' \item{score}{Ranking metric (global).}
#' \item{poth_bestk}{Best k treatments POTH.}
#'
#' @examples
#' \donttest{
#' library("netmeta")
#' data(Senn2013)
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013, sm = "MD", random = FALSE)
#'
#' # Best k treatments method
#' bk1 <- bestk(poth(net1))
#' bk1
#' plot(bk1)
#' plot(bk1, labels = TRUE)
#' bk2 <- bestk(poth(net1), sort = FALSE)
#' bk2
#' plot(bk2)
#' plot(bk2, labels = TRUE)
#' }
#'
#' @rdname bestk
#' @method bestk poth
#' @export

bestk.poth <- function(x, sort = TRUE, ...) {

  chkclass(x, "poth")

  n <- x$n
  trts <- x$trts

  if (x$input == "mcmc.samples") {
    score_type <- "SUCRA"
    ranking <- x$ranking
    small.values <- x$small.values
    #
    samples <- x$x
    # gives order from best treatment to worst
    if (sort)
      seq <- order(ranking, decreasing = TRUE)
    else
      seq <- seq_along(ranking)
    #
    cum_rps <-
      lapply(2:n,
             function(x)
               rankMCMC(samples[, seq[1:x]], small.values))
    #
    cum_poths <- c(NA, sapply(cum_rps, function(x) poth(x)$poth))
  }
  else if (x$input %in% c("effects.se", "netmeta")) {
    score_type <- "P-score"
    ranking <- x$ranking
    small.values <- x$small.values
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
    # gives order from best treatment to worst
    if (sort)
      seq <- order(ranking, decreasing = TRUE)
    else
      seq <- seq_along(ranking)
    #
    cum_pscores <-
      lapply(2:n,
             function(x)
               pscores(TE[seq[1:x], seq[1:x]],
                       seTE[seq[1:x], seq[1:x]],
                       small.values))
    #
    cum_poths <- c(NA, sapply(cum_pscores, function(x) poth(x)$poth))
  }
  else
    stop("Best k treatments method not available for input type '", x$input, "'.")
  #
  ranking <- ranking[seq]
  trts <- trts[seq]

  res <- data.frame(trt = trts,
                    rank = rank(-ranking),
                    score = ranking,
                    poth_bestk = cum_poths)
  rownames(res) <- trts
  #
  attr(res, "poth") <- x$poth
  attr(res, "pooled") <- x$pooled
  #
  attr(res, "sort") <- sort
  attr(res, "score_type") <- score_type
  attr(res, "small.values") <- small.values
  #
  class(res) <- c("bestk.poth", class(res))
  #
  res
}


#' @rdname bestk
#' @keywords print
#' @method print bestk.poth
#' @export

print.bestk.poth <- function(x, digits = 3, legend = TRUE, ...) {

  chkclass(x, "bestk.poth")
  #
  chknumeric(digits, min = 0, length = 1)
  chklogical(legend)
  #
  poth <- attr(x, "poth")
  score_type <- attr(x, "score_type")
  pooled <- attr(x, "pooled")

  txt <- "Best k treatments precision of treatment hierarchy (POTH)"
  #
  if (pooled != "")
    txt <- paste0(txt,
                  " based on ",
                  if (pooled == "common") "common" else "random",
                  " effects model")
  #
  txt <- paste0(txt, "\n\n")
  #
  cat(txt)

  x$trt <- NULL
  #
  x$score <- round(x$score, digits)
  #
  x$poth_bestk <- round(x$poth_bestk, digits)
  x$poth_bestk <- ifelse(is.na(x$poth_bestk), rep_len("", nrow(x)), x$poth_bestk)
  #
  rownames(x) <- paste("Adding", rownames(x))
  names(x)[names(x) == "poth_bestk"] <- "POTH"
  #
  class(x) <- "data.frame"
  #
  print(x)

  if (legend) {
    cat("\nLegend:\n")
    cat(" rank  - Treatment rank (global)\n")
    cat(" score - Ranking metric (global)\n")
    cat(" POTH - Best k Treatments POTH\n")
  }

  invisible(NULL)
}


#' @rdname bestk
#' @export bestk

bestk <- function(x, ...)
  UseMethod("bestk")
