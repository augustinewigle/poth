#' Cumulative method for precision of treatment hierarchy (POTH) metric
#'
#' @param x An R object of class \code{poth}.
#' @param sort A logical indicating whether results should be sorted
#'   by decreasing ranking metric.
#' @param \dots Additional arguments.
#'
#' @return A vector of POTH values for the top 2 to n treatments with additional class \code{cumul.poth}
#'
#' @rdname cumul
#' @method cumul poth
#' @export

cumul.poth <- function(x, sort = TRUE, ...) {

  chkclass(x, "poth")

  n <- x$n
  trts <- x$trts

  if (x$input == "mcmc.samples") {
    score_type <- "SUCRA"
    ranking <- x$ranking
    samples <- x$x
    small.values <- x$small.values
    #

    # gives order from best treatment to worst
    seq <- order(ranking, decreasing = T)

    #
    cum_rps <-
      lapply(2:n,
             function(x)
               rankMCMC(samples[, seq[1:x]], small.values))

    #
    cum_poths <- sapply(cum_rps, function(x) poth(x)$poth)
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
    seq <- order(ranking, decreasing = TRUE)

    cum_pscores <-
      lapply(2:n,
             function(x)
               pscores(TE[seq[1:x], seq[1:x]], seTE[seq[1:x], seq[1:x]], small.values))
    #
    cum_poths <- sapply(cum_pscores, function(x) poth(x)$poth)

  }
  else
    stop("Leave-one-out method not available for input type '", x$input, "'.")

  names(cum_poths) <- paste0("top ",2:n)

  res <- cum_poths

  grps <- sapply(2:n, function(x) paste(trts[seq[1:x]], collapse = ", "))

  res <- list(cumul.poth = cum_poths,
              groups = grps,
              ranking = x$ranking,
              poth = x$poth)

  # attr(res, "score_type") <- score_type
  # attr(res, "ranking") <- x$ranking
  # attr(res, "groups") <- grps
  #
  class(res) <- c("cumul.poth", class(res))
  #

  return(res)
}

#' @rdname cumul
#' @export cumul

cumul <- function(x, ...)
  UseMethod("cumul")
