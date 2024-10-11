#' Calculate P-scores from a set of relative effects and standard errors
#'
#' @param TE Matrix of relative effects
#' @param seTE Matrix of estimated standard errors for relative effects
#' @param trts optional; vector of treatment names matching order in TE and sds
#' @param small.values A character string specifying whether small
#'   outcome values indicate a beneficial (\code{"desirable"}) or
#'   harmful (\code{"undesirable"}) effect, can be abbreviated.
#'
#' @return named vector of P-scores


pscores <- function(TE, seTE, small.values = "desirable", trts = NULL) {

  small.values <- setsv(small.values)

  n <- nrow(TE)
  n.seq <- seq_len(n)

  # name check
  if (length(trts) != n) {
    if (is.null(colnames(TE))) {
      trts <- paste0("trt", n.seq)
      colnames(TE) <- colnames(seTE) <- trts
    }
    #
    trts <- colnames(TE)
  }
  else
    colnames(TE) <- trts

  a_mat <- matrix(NA, nrow = n, ncol = n)
  #
  for (i in n.seq)
    for (j in n.seq)
      if (i != j)
        a_mat[i, j] <- TE[i, j] / seTE[i, j]


  # Calculate p-scores
  pscores <- numeric(n)
  direction <- ifelse(small.values == "undesirable", 1, -1)
  #
  if (n == 1)
    pscores <- 1
  else {
    for (i in n.seq)
      pscores[i] <- sum(pnorm(a_mat[i, ] * direction), na.rm = TRUE) / (n - 1)
  }
  #
  names(pscores) <- trts

  pscores
}
