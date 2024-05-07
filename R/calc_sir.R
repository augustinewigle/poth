#' Calculate separation in ranking (SIR) metric
#' 
#' @description
#' Separation in ranking (SIR) is a metric to quantify the uncertainty in
#' a treatment hierarchy in network neta-analysis
#' 
#' @param x A square matrix with rankings or a vector of SUCRA values
#'   or P-scores (see Details).
#' @param trts An optional vector with treatment names. Must match the
#'   order of treatments provided for argument \code{x}.
#'
#' @details
#' Argument \code{x} can be
#' \itemize{
#'  \item a square matrix where the rows represent ranks and the columns
#'    treatments,
#'  \item a vector with SUCRAs or P-scores.
#' }
#'
#' @return A list with elements 'sucras' and 'sir'.
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
#' # Make results reproducible
#' set.seed(1909)
#' 
#' # R package netmeta must be available
#' library("netmeta")
#' data(smokingcessation)
#' p1 <- pairwise(list(treat1, treat2, treat3),
#'                event = list(event1, event2, event3), n = list(n1, n2, n3),
#'                data = smokingcessation, sm = "OR")
#' net1 <- netmeta(p1)
#' 
#' rg1 <- rankogram(net1)
#' rg1$ranking.matrix.common
#' t(rg1$ranking.matrix.common)
#' 
#' netrank(rg1, random = FALSE)
#' sir(t(rg1$ranking.matrix.common))
#' sir(netrank(rg1)$ranking.common)
#' 
#' data(Senn2013)
#' net2 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'                 data = Senn2013, sm = "MD", random = FALSE, nchar.trts = 4)
#' 
#' rg2 <- rankogram(net2)
#' rg2$ranking.matrix.common
#' t(rg2$ranking.matrix.common)
#' 
#' netrank(rg2, random = FALSE)
#' sir(t(rg2$ranking.matrix.common))
#' sir(netrank(rg2)$ranking.common)
#' }
#' 
#' @export sir

sir <- function(x, trts = NULL) {

  if (is.matrix(x)) {
    # error checking
    if (nrow(x) != ncol(x))
      stop("Argument 'x' must be a square matrix.")
    else if (any(abs(apply(x, 1, sum) - 1) > 1e-7))
      warning("The rows of a ranking matrix should sum to 1.")
    else if (any(abs(apply(x, 2, sum) - 1) > 1e-7))
      warning("The columns of a ranking matrix should sum to 1.")

    n <- ncol(x)

    rank_e <- apply(seq(1:n) * x, 2, sum)

    sucras <- (n - rank_e) / (n - 1)


    rank_vars <- numeric(n)
    
    for (i in 1:n)
      rank_vars[i] <- sum((seq(1:n) - rank_e[i])^2 * x[, i])
    
    sir <- 1 - sum(rank_vars) * 12 / n / (n + 1) / (n - 1)
  }
  else if (is.vector(x) & !is.list(x)) {
    # error checking
    if (abs(mean(x)) - 0.5 > 1e-7)
      warning("The mean of the sucras should be 0.5. Check your values.")
    
    n <- length(x)
    
    sucras <- x
    sir <- 3 * ((1 - n) / (n + 1) + 4 * (n - 1) / (n * (n + 1)) * sum(sucras^2))
  }
  else
    stop("Argument 'x' must be a ranking matrix or vector.")
  
  if (!is.null(trts)) {
    if (length(trts) != length(sucras))
      stop("Different number of treatment names and SUCRAs / P-scores.")
    #
    names(sucras) <- trts
  }
  
  return(list(sucras = sort(sucras, decreasing = TRUE), sir = sir))
}
