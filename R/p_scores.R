#' Calculate p-scores from a set of relative effects and standard errors
#'
#' @param diffs Matrix of relative effects
#' @param ses Matrix of estimated standard errors for relative effects
#' @param trts optional; vector of treatment names matching order in diffs and sds
#' @param largerbetter Logical, T if larger outcomes are good. Defauls is T
#'
#' @return named vector of P-scores
#'
#' @importFrom stats pnorm
#'
#' @export

pscores <- function(diffs,
                    ses,
                    largerbetter = TRUE,
                    trts = NULL) {

  n <- nrow(diffs)

  # name check
  if (length(trts) != n) {
    if (is.null(colnames(diffs))) {
      warning("Using generic treatment names.")
      trts <- paste0("trt", 1:n)
      colnames(diffs) <- colnames(ses) <- trts
    }
    #
    trts <- colnames(diffs)
  }
  
  a_mat <- matrix(NA, nrow = n, ncol = n)
  
  for (i in 1:n)
    for (j in 1:n)
      if (i != j)
        a_mat[i,j] <- diffs[i,j]/ses[i,j]

  dir <- ifelse(largerbetter, 1, -1)

  pscores <- numeric(n)

  # calculate p-scores
  for (i in 1:n)
    pscores[i] <- sum(pnorm(a_mat[i,]*dir), na.rm = T)/(n-1)
  
  names(pscores) <- trts

  pscores
}
