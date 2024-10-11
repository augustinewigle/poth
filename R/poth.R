#' Calculate precision of treatment hierarchy (POTH) metric
#'
#' @description
#' Precision of treatment hierarchy (POTH) is a metric to quantify the uncertainty in
#' a treatment hierarchy in network meta-analysis
#'
#' @param x Mandatory argument with suitable information on the treatment
#'   hierarchy (see Details).
#' @param se Matrix of estimated standard errors for relative effects.
#' @param small.values A character string specifying whether small
#'   outcome values indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated.
#' @param pooled A character string indicating whether the treatment hierarchy
#'   is based on a common or random effects model. Either \code{"common"} or
#'   \code{"random"}, can be abbreviated.
#' @param trts An optional vector with treatment names. Must match the
#'   order of treatments provided for argument \code{x}.
#' @param sort A logical indicating whether printout should be sorted
#'   by decreasing ranking metric.
#' @param digits Minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param object An object of class \code{summary.poth}.
#' @param \dots Additional arguments (ignored).
#'
#' @details
#' This function calculates the precision of treatment hierarchy (POTH) metric to quantify
#' the uncertainty in a treatment hierarchy in network meta-analysis
#' (Wigle et al., 2024).
#'
#' Argument \code{x} providing information on the treatment hierarchy is the
#' only mandatory argument. The following input formats can be provided:
#' \enumerate{
#'  \item vector representing a ranking metric, i.e., SUCRAs or P-scores,
#'  \item square matrix with the probabilities for each possible rank
#'  (with treatments in rows and ranks in columns),
#'  \item MCMC samples (with samples in rows and treatments in columns),
#'  \item relative effect matrix,
#'  \item R object created with \code{\link[netmeta]{netmeta}},
#'    \code{\link[netmeta]{netrank}}, or \code{\link[netmeta]{rankogram}}
#'    object from R package \bold{netmeta}.
#' }
#'
#' Argument \code{se} must be provided if argument \code{x} is a matrix with
#' relative effects. Otherwise, argument \code{se} is ignored.
#'
#' Argument \code{small.values} must be provided if argument \code{x} contains
#' MCMC samples, relative effects, or is an object created with
#' \code{\link[netmeta]{netmeta}}. This argument can be provided for an R
#' object created with \code{\link[netmeta]{netrank}} or
#' \code{\link[netmeta]{rankogram}} and is ignored otherwise.
#'
#' Argument \code{trts} is ignored for \code{\link[netmeta]{netmeta}},
#' \code{\link[netmeta]{netrank}}, and \code{\link[netmeta]{rankogram}} objects.
#'
#' @return
#' An object of class \code{poth} with corresponding \code{print}
#' function. The object is a list containing the following components:
#' \item{poth}{Separation in ranking metric.}
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
#' Wigle, A., Béliveau, A., Salanti, G., Rücker, G., Schwarzer, G., Mavridis, D.,
#' Nikolakopoulou, A. (2024):
#' Precision of Treatment Hierarchy: A Metric for Quantifying Uncertainty in Treatment
#' Hierarchies in Network Meta-Analysis
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
#' # Calculate probabilities for each possible rank
#' set.seed(1909) # make results reproducible
#' rg1 <- rankogram(net1)
#' rg1
#'
#' # Calculate POTH
#' s1 <- poth(rg1)
#' s1
#'
#' # Also print probabilities for each possible rank
#' summary(s1)
#'
#' # Use SUCRAs to calculate POTH
#' nr1 <- netrank(rg1)
#' nr1
#' poth(nr1)
#' poth(nr1$ranking.common)
#'
#' data(Senn2013)
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", random = FALSE)
#'
#' # Use P-scores to calculate POTH
#' nr2 <- netrank(net2)
#' nr2
#' poth(nr2)
#' }
#'
#' @export poth

poth <- function(x, se = NULL, small.values, pooled, trts = NULL) {

  #
  # Set arguments
  #

  if (missing(small.values)) {
    if (inherits(x, c("netmeta", "netrank", "rankogram")))
      small.values <- x$small.values
    else
      small.values <- "desirable"
  }
  #
  small.values <- setsv(small.values)
  #
  if (!missing(pooled)) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    pooled[pooled == "fixed"] <- "common"
  }
  else if (inherits(x, c("netmeta", "netrank", "rankogram"))) {
    if (!x$common & x$random)
      pooled <- "random"
    else
      pooled <- "common"
  }
  else
    pooled <- ""


  #
  # Calculate POTH
  #

  if ((is.matrix(x) & !missing(se)) || inherits(x, "netmeta")) {
    #
    # Input: matrices with relative effects and standard errors
    #
    if (inherits(x, "netmeta")) {
      if (pooled == "common") {
        TE <- x$TE.common
        se <- x$seTE.common
      }
      else {
        TE <- x$TE.random
        se <- x$seTE.random
      }
      #
      input <- "netmeta"
      #
      trts <- NULL
    }
    else {
      if (!is.matrix(se))
        stop("Argument 'se' must be a matrix with standard errors.")
      #
      TE <- x
      #
      input <- "effects.se"
    }
    #
    ranking.matrix <- NULL
    #
    n <- ncol(TE)
    ranking <- pscores(TE, se, small.values, trts)
    #
    poth <-
      3 * ((1 - n) / (n + 1) + 4 * (n - 1) / (n * (n + 1)) * sum(ranking^2))
  }
  else if ((is.matrix(x) | is.data.frame(x)) & missing(se)) {
    if (ncol(x) == nrow(x)) {
      #
      # Input: ranking matrix
      #
      if (any(abs(apply(x, 1, sum) - 1) > 1e-7))
        warning("The rows of a ranking matrix should sum to 1.")
      else if (any(abs(apply(x, 2, sum) - 1) > 1e-7))
        warning("The columns of a ranking matrix should sum to 1.")
      #
      ranking.matrix <- x
      #
      input <- "ranking.matrix"
    }
    else {
      #
      # Input: MCMC samples
      #
      ranking.matrix <- rankMCMC(x, small.values, trts)
      #
      input <- "mcmc.samples"
    }
    #
    n <- ncol(x)
    n.seq <- seq_len(n)
    #
    rank_e <- apply(col(ranking.matrix) * ranking.matrix, 1, sum)
    #
    rank_vars <- vector("numeric", n)
    #
    for (i in n.seq)
      rank_vars[i] <- sum((n.seq - rank_e[i])^2 * ranking.matrix[i, ])
    #
    ranking <- (n - rank_e) / (n - 1)
    #
    poth <- 1 - sum(rank_vars) * 12 / n / (n + 1) / (n - 1)
  }
  else if ((is.vector(x) & !is.list(x)) ||
           inherits(x, c("netrank", "rankogram"))) {
    #
    # Input: ranking metrix or R object created with netrank() or rankogram()
    #
    if (inherits(x, "rankogram")) {
      n <- x$x$n
      #
      if (pooled == "common") {
        ranking <- x$ranking.common
        ranking.matrix <- x$ranking.matrix.common
      }
      else {
        ranking <- x$ranking.random
        ranking.matrix <- x$ranking.matrix.random
      }
      #
      input <- "rankogram"
      #
      trts <- NULL
    }
    else if (inherits(x, "netrank")) {
      n <- x$x$n
      #
      if (pooled == "common") {
        ranking <- x$ranking.common
        ranking.matrix <- x$ranking.matrix.common
      }
      else {
        ranking <- x$ranking.random
        ranking.matrix <- x$ranking.matrix.random
      }
      #
      input <- "netrank"
    }
    else {
      if (abs(mean(x)) - 0.5 > 1e-7)
        warning("The mean of the ranking should be 0.5. Check your values.")
      #
      n <- length(x)
      #
      ranking <- x
      ranking.matrix <- NULL
      #
      input <- "ranking"
      #
      trts <- NULL
    }
    #
    poth <-
      3 * ((1 - n) / (n + 1) + 4 * (n - 1) / (n * (n + 1)) * sum(ranking^2))
  }
  else
    stop("Argument 'x' must be a ranking vector, matrix or vector.")


  if (!is.null(trts)) {
    if (length(trts) != length(ranking))
      stop("Different number of treatment names and rankings.")
    #
    names(ranking) <- trts
    #
    if (!is.null(ranking.matrix))
      rownames(ranking.matrix) <- trts
    #
    if (input == "mcmc.samples")
      colnames(x) <- trts
    #
    if (input == "effects.se") {
      rownames(TE) <- colnames(TE) <- trts
      rownames(se) <- colnames(se) <- trts
    }
  }
  else
    trts <- names(ranking)

  res <- list(poth = poth, ranking = ranking, ranking.matrix = ranking.matrix,
              small.values = small.values, pooled = pooled,
              n = n,
              x = x, se = se, input = input, trts = trts)
  #
  if (inherits(x, "netmeta"))
    res$TE <- TE
  #
  class(res) <- "poth"
  #
  res
}


#' @rdname poth
#' @keywords print
#' @method print poth
#' @export

print.poth <- function(x, sort = TRUE, digits = 3, ...) {

  chkclass(x, "poth")

  class(x) <- "list"
  #
  if (sort)
    seq <- rev(order(x$ranking))
  else
    seq <- seq_along(x$ranking)
  #
  x$poth <- round(x$poth, digits)
  x$ranking <- round(x$ranking[seq], digits)

  txt <- "Precision of treatment hierarchy (POTH)"
  #
  if (x$pooled != "")
    txt <- paste0(txt,
                  " based on ",
                  if (x$pooled == "common") "common" else "random",
                  " effects model")
  #
  txt <- paste0(txt, "\n\n")
  #
  cat(txt)
  #
  cat(paste0("POTH = ", x$poth, "\n\n"))
  #
  cat("Ranking:\n")
  print(x$ranking)

  invisible(NULL)
}


#' @rdname poth
#' @keywords summary
#' @method summary poth
#' @export

summary.poth <- function(object, ...) {
  res <- object
  class(res) <- "summary.poth"
  res
}


#' @rdname poth
#' @keywords print
#' @method print summary.poth
#' @export

print.summary.poth <- function(x, sort = TRUE, digits = 3, ...) {

  chkclass(x, "summary.poth")

  if (sort)
    seq <- rev(order(x$ranking))
  else
    seq <- seq_along(x$ranking)
  #
  x$poth <- round(x$poth, digits)
  x$ranking <- round(x$ranking[seq], digits)
  #
  if (!is.null(x$ranking.matrix))
    x$ranking.matrix <- round(x$ranking.matrix[seq, ], digits)

  txt <- "Precision of treatment hierarchy (POTH)"
  #
  if (x$pooled != "")
    txt <- paste0(txt,
                  " based on ",
                  if (x$pooled == "common") "common" else "random",
                  " effects model")
  #
  txt <- paste0(txt, "\n\n")
  #
  cat(txt)
  #
  cat(paste0("POTH = ", x$poth, "\n\n"))
  #
  cat("Ranking:\n")
  print(x$ranking)
  #
  if (!is.null(x$ranking.matrix)) {
    cat("Ranking matrix:\n")
    print(x$ranking.matrix)
  }

  invisible(NULL)
}
