#' Calculate a ranking probabilities matrix from mcmc samples
#' @param samples a matrix or data.frame of MCMC samples, where rows are MCMC samples and columns are relative effects (relative to anchor) for treatments.
#' must have column names that are the name of each treatment.
#' @export
rank_mat <- function(samples, positive_good) {

  # samples <- combine.mcmc(samples)
  #
  # samples <- samples[,startsWith(colnames(samples), "d[")] # get rid of sigma and other variables

  nt <- ncol(samples)

  rank_mat <- matrix(nrow = nt, ncol = nt)

  # Sort every row of the matrix:
  sorted <- t(apply(samples, 1, function(x) {sort(x, index.return = T, decreasing = positive_good)$ix}))
  colnames(sorted) <- colnames(samples)

  # columns of rank_mat are treatments, rows are ranks
  for(i in 1:nt) {

    for(j in 1:nt) {

      # check how often its the j'th rank (1/0 vector) and calculate mean of j'th ranks
      rank_mat[j, i] <- sum(sorted[,i] == j)/nrow(sorted)
      # sum(apply(samples, 1, function(x) {sort(x, index.return = T, decreasing = positive_good)$ix[j] == i}))/

    }

  }

  colnames(rank_mat) <- colnames(samples)
  rownames(rank_mat) <- 1:nrow(rank_mat)

  return(rank_mat)

}
