#' Calculate a ranking probabilities matrix from MCMC samples
#' 
#' @param x a matrix or data.frame of MCMC samples, where rows are MCMC samples and columns are relative effects (relative to anchor) for treatments.
#' must have column names that are the name of each treatment.
#' @param small.values A character string specifying whether small
#'   outcome values indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated.
#' @param trts character vector of treatment names, optional if samples has column names
#' 
#' @export

rankMCMC <- function(x, small.values = "desirable", trts = NULL) {
  
  small.values <- setsv(small.values)
  
  n <- ncol(x)
  n.seq <- seq_len(n)
  
  # name check
  if (length(trts) != n) {
    if (is.null(colnames(x))) {
      trts <- paste0("trt", n.seq)
      colnames(x) <- trts
    }
    #
    trts <- colnames(x)
  }
  else
    colnames(x) <- trts
  
  # Ranks for every row of the matrix:
  direction <- ifelse(small.values == "undesirable", -1, 1)
  ranks <- t(apply(x * direction, 1, rank))
  colnames(ranks) <- colnames(x)
  
  # rows of rank_mat are treatments, columns are ranks
  rank_mat <- matrix(nrow = n, ncol = n)
  #
  for (i in n.seq)
    for (j in n.seq)
      rank_mat[i, j] <- mean(ranks[, i] == j)
  #
  rownames(rank_mat) <- colnames(x)
  colnames(rank_mat) <- seq_len(nrow(rank_mat))
  
  rank_mat
}
