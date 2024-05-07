#' Calculate a ranking probabilities matrix from MCMC samples
#' 
#' @param samples a matrix or data.frame of MCMC samples, where rows are MCMC samples and columns are relative effects (relative to anchor) for treatments.
#' must have column names that are the name of each treatment.
#' @param largerbetter logical indicating if larger values indicate a more effective treatment
#' @param trts character vector of treatment names, optional if samples has column names
#' 
#' @export

rank_mat <- function(samples, largerbetter, trts = NA) {

  nt <- ncol(samples)

  # name check
  if (length(trts) != nt) {

    if(is.null(colnames(samples))) {

      warning("Please specify treatment names. Assigning generic treatment names")

      trts <- paste0("trt", 1:nt)

      colnames(samples) <- trts

    }

    trts <- colnames(samples)

  } else {

    colnames(samples) <- trts

  }

  nt <- ncol(samples)

  rank_mat <- matrix(nrow = nt, ncol = nt)

  # Ranks for every row of the matrix:
  sorted <- t(apply(samples * ifelse(largerbetter, -1, 1), 1, rank))
  colnames(sorted) <- colnames(samples)

  # columns of rank_mat are treatments, rows are ranks
  for (i in 1:nt)
    for (j in 1:nt)
      rank_mat[j, i] <- mean(sorted[, i] == j)
  #
  colnames(rank_mat) <- colnames(samples)
  rownames(rank_mat) <- 1:nrow(rank_mat)
  
  rank_mat
}
