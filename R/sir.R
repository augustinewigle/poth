#' Calculate separation in ranking (SIR) metric
#' 
#' @description
#' Separation in ranking (SIR) is a metric to quantify the uncertainty in
#' a treatment hierarchy in network meta-analysis
#' 
#' @param x Mandatory argument with suitable information on the treatment
#'   hierarchy (see Details).
#' @param pooled A character string indicating whether the treatment hierarchy
#'   is based on a common or random effects model. Either \code{"common"} or
#'   \code{"random"}, can be abbreviated.
#' @param trts An optional vector with treatment names. Must match the
#'   order of treatments provided for argument \code{x}.
#' @param sort A logical indicating whether printout should be sorted
#'   by decreasing ranking metric.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param object An object of class \code{summary.sir}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' This function calculates the separation in ranking (SIR) metric to quantify 
#' the uncertainty in a treatment hierarchy in network meta-analysis
#' (Wiggle et al., 2024).
#'   
#' Argument \code{x} providing information on the treatment hierarchy is the
#' only mandatory argument. The following input formats can be provided:
#' \itemize{
#'  \item vector representing a ranking metric, i.e., SUCRAs or P-scores,
#'  \item square matrix with the probabilities for each possible rank
#'  (with treatments in rows and ranks in columns).
#' }
#' It is also possible to provide an R object created with
#' \code{\link[netmeta]{netrank}} (ranking metric) or
#' \code{\link[netmeta]{rankogram}} (probabilities for each possible rank)
#' from R package \bold{netmeta}.
#'
#' @return
#' An object of class \code{sir} with corresponding \code{print}
#' function. The object is a list containing the following components:
#' \item{sir}{Separation in ranking metric.}
#' \item{ranking}{A named numeric vector with rankings, i.e.,
#'   SUCRAs or P-scores.}
#' \item{ranking.matrix}{A square matrix with the probabilities
#'   for each possible rank (if information is available).}
#' \item{pooled}{As defined above.}
#'
#' @author Augustine Wigle \email{amhwigle@uwaterloo.ca},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @references
#' Wiggle A, ... (2024):
#' Separation In Ranking: A Metric for Quantifying Uncertainty in Treatment
#' Hierarchies in Network Meta-Analysis
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
#' # Calculate probabilities for each possible rank
#' set.seed(1909) # make results reproducible
#' rg1 <- rankogram(net1)
#' rg1
#' 
#' # Calculate SIR
#' s1 <- sir(rg1)
#' s1
#' 
#' # Also print probabilities for each possible rank
#' summary(s1)
#' 
#' # Use SUCRAs to calculate SIR
#' nr1 <- netrank(rg1)
#' nr1
#' sir(nr1)
#' sir(nr1$ranking.common)
#' 
#' data(Senn2013)
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", random = FALSE)
#' 
#' # Use P-scores to calculate SIR
#' nr2 <- netrank(net2)
#' nr2
#' sir(nr2)
#' }
#' 
#' @export sir

sir <- function(x, pooled, trts = NULL) {
  
  #
  # Set argument pooled
  #
  if (!missing(pooled)) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    pooled[pooled == "fixed"] <- "common"
  }
  else if (inherits(x, c("netrank", "rankogram"))) {
    if (!x$common & x$random)
      pooled <- "random"
    else
      pooled <- "common"
  }
  else
    pooled <- ""
  
  
  #
  # Calculate SIR
  #
  
  if (is.matrix(x)) {
    # Error checking
    if (nrow(x) != ncol(x))
      stop("Argument 'x' must be a square ranking matrix.")
    else if (any(abs(apply(x, 1, sum) - 1) > 1e-7))
      warning("The rows of a ranking matrix should sum to 1.")
    else if (any(abs(apply(x, 2, sum) - 1) > 1e-7))
      warning("The columns of a ranking matrix should sum to 1.")
    #
    ranking.matrix <- x
    #
    n <- nrow(x)
    n.seq <- seq_len(n)
    rank_e <- apply(n.seq * x, 1, sum)
    #
    rank_vars <- vector("numeric", n)
    #
    for (i in 1:n)
      rank_vars[i] <- sum((n.seq - rank_e[i])^2 * x[i, ])
    #
    ranking <- (n - rank_e) / (n - 1)
    #
    sir <- 1 - sum(rank_vars) * 12 / n / (n + 1) / (n - 1)
  }
  else if ((is.vector(x) & !is.list(x)) ||
           inherits(x, c("netrank", "rankogram"))) {
    #
    ranking.matrix <- NULL
    #
    if (inherits(x, "rankogram")) {
      n <- x$x$n
      #
      if (pooled == "common") {
        ranking <- x$ranking.common
        ranking.matrix <- x$ranking.matrix.common
        #
        x <- x$ranking.matrix.common
      }
      else {
        ranking <- x$ranking.random
        ranking.matrix <- x$ranking.matrix.random
        #
        x <- x$ranking.matrix.random
      }
    }
    else if (inherits(x, "netrank")) {
      n <- x$x$n
      #
      if (pooled == "common") {
        ranking <- x$ranking.common
        #
        x <- x$ranking.matrix.common
      }
      else {
        ranking <- x$ranking.random
        ranking.matrix <- x$ranking.matrix.random
        #
        x <- x$ranking.matrix.random
      }
    }
    else {
      # error checking
      if (abs(mean(x)) - 0.5 > 1e-7)
        warning("The mean of the ranking should be 0.5. Check your values.")
      #
      n <- length(x)
      ranking <- x
    }
    #
    sir <- 3 * ((1 - n) / (n + 1) + 4 * (n - 1) / (n * (n + 1)) * sum(ranking^2))
  }
  else
    stop("Argument 'x' must be a ranking matrix or vector.")
  
  if (!is.null(trts)) {
    if (length(trts) != length(ranking))
      stop("Different number of treatment names and rankings.")
    #
    names(ranking) <- trts
  }
  
  res <- list(sir = sir, ranking = ranking, ranking.matrix = ranking.matrix,
              pooled = pooled)
  class(res) <- "sir"
  #
  res
}

#' @rdname sir
#' @keywords print
#' @method print sir
#' @export

print.sir <- function(x, sort = TRUE, digits = 3, ...) {
  class(x) <- "list"
  #
  if (sort)
    seq <- rev(order(x$ranking))
  else
    seq <- seq_along(x$ranking)
  #
  x$sir <- round(x$sir, digits)
  x$ranking <- round(x$ranking[seq], digits)
  #
  x$ranking.matrix <- NULL
  #
  if (x$pooled == "")
    x$pooled <- NULL
  #
  print(x)
  #
  invisible(NULL)
}

#' @rdname sir
#' @keywords summary
#' @method summary sir
#' @export

summary.sir <- function(object, ...) {
  res <- object
  class(res) <- "summary.sir"
  res
}


#' @rdname sir
#' @keywords print
#' @method print summary.sir
#' @export

print.summary.sir <- function(x, sort = TRUE, digits = 3, ...) {
  class(x) <- "list"
  #
  if (sort)
    seq <- rev(order(x$ranking))
  else
    seq <- seq_along(x$ranking)
  #
  x$sir <- round(x$sir, digits)
  x$ranking <- round(x$ranking[seq], digits)
  #
  if (!is.null(x$ranking.matrix))
    x$ranking.matrix <- round(x$ranking.matrix[seq, ], digits)
  #
  if (x$pooled == "")
    x$pooled <- NULL
  #
  print(x)
  #
  invisible(NULL)
}
