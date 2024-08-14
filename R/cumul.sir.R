#' Cumulative method for separation in ranking (SIR) metric
#'
#' @param x An R object of class \code{sir}.
#' @param sort A logical indicating whether results should be sorted
#'   by decreasing ranking metric.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param \dots Additional arguments.
#'
#' @return A vector of sir values for the top 2 to n treatments with additional class \code{cumul.sir}
#'
#' @rdname cumul
#' @method cumul sir
#' @export

cumul.sir <- function(x, sort = TRUE, ...) {

  chkclass(x, "sir")

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
    cum_sirs <- sapply(cum_rps, function(x) sir(x)$sir)
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
    cum_sirs <- sapply(cum_pscores, function(x) sir(x)$sir)

  }
  else
    stop("Leave-one-out method not available for input type '", x$input, "'.")

  names(cum_sirs) <- paste0("top ",2:n)

  res <- cum_sirs

  grps <- sapply(2:n, function(x) paste(trts[seq[1:x]], collapse = ", "))

  res <- list(cumul.sir = cum_sirs,
              groups = grps,
              ranking = x$ranking,
              sir = x$sir)

  # attr(res, "score_type") <- score_type
  # attr(res, "ranking") <- x$ranking
  # attr(res, "groups") <- grps
  #
  class(res) <- c("cumul.sir", class(res))
  #

  return(res)
}

#' @rdname cumul
#' @export cumul

cumul <- function(x, ...)
  UseMethod("cumul")
