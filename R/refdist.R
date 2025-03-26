#' Generate reference distribution for POTH for a given network structure
#'
#' @param x A \code{\link[netmeta]{netmeta}} object.
#' @param d A vector of the desired relative effects, must be in the same
#'   order as \code{x$trts}.
#' @param pooled A character string indicating whether the treatment hierarchy
#'   is based on a common or random effects model. Either \code{"common"} or
#'   \code{"random"}, can be abbreviated.
#' @param nsim Number of samples from reference distribution.
#' @param verbose A logical indicating whether progress information should
#'   be printed.
#'
#' @details
#' By default, argument \code{pooled} is equal to "random" if only the random
#' effects model was considered in the network meta-analysis \code{x}.
#' Otherwise, argument \code{pooled} is equal to "common".
#'
#' If argument \code{d} is missing, the respective relative effects are taken
#' to be all 0.
#'
#' @return A vector of POTH values.
#'
#' @seealso \code{\link[netmeta]{netmeta}}
#'
#' @examples
#' \donttest{
#' library("netmeta")
#' data(Senn2013)
#' net1 <- netmeta(TE, seTE, treat1.long, treat2.long, studlab,
#'   data = Senn2013,
#'   sm = "MD")
#'
#' # POTH (based on common effects model)
#' poth(net1)
#'
#' # Sample POTH values from reference distribution (common effects model)
#' set.seed(1909)
#' poths <- refdist(net1)
#' summary(poths)
#'
#' # POTH (based on random effects model)
#' poth(net1, pooled = "random")
#'
#' # Sample POTH values from reference distribution (common effect model)
#' poths.r <- refdist(net1, pooled = "random")
#' summary(poths.r)
#' }
#'
#' @export

refdist <- function(x, d, pooled, nsim = 25, verbose = TRUE) {

  chkclass(x, "netmeta")
  #
  if (!missing(pooled)) {
    pooled <- setchar(pooled, c("common", "random", "fixed"))
    pooled[pooled == "fixed"] <- "common"
  }
  else {
    if (!x$common & x$random)
      pooled <- "random"
    else
      pooled <- "common"
  }
  #
  if (missing(d)) {
    d <- rep(0, x$n)
  }
  else
    chknumeric(d, length = length(x$trts))
  #
  chknumeric(nsim, min = 1, length = 1)
  chklogical(verbose)

  # Simulate data with the desired relative effects and identical structure
  # and heterogeneity
  #
  meanvec <- x$X.matrix %*% d

  # Standard errors based on common or random effects model, ignoring multi-arm corrections
  #
  if (pooled == "random")
    sdvec <- x$seTE.adj.random
  else
    sdvec <- x$seTE.adj.common

  simdata <- replicate(nsim,
                       rnorm(length(meanvec), mean = meanvec, sd = sdvec))

  # Calculate POTHs
  #
  poths <- numeric(nsim)
  #
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  #
  for (i in seq_len(nsim)) {

    if(any(x$multiarm)) { # correct data to be consistent

      mstudies <- x$studlab[x$multiarm]

      # loop over studies
      for(study in mstudies) {

        ix <- which(x$studlab == study)
        narm <- unique(x$n.arms[ix])
        basicix  <- ix[1:(narm-1)] # indices for first ai-a contrasts for this study
        funcix <- ix[narm:length(ix)] # indices for remaining contrasts to be made internally consistent

        A <- t(x$X.matrix[basicix,])
        if(narm == 3) {

          B <- x$X.matrix[funcix,]

        } else {

          B <- t(x$X.matrix[funcix,])

        }

        simdata[funcix,i] <- matrix(simdata[basicix,i], nrow = 1) %*% MASS::ginv(A)%*%B

      }

    }
    poths[i] <-
      poth(netmeta(TE = simdata[, i],
                   seTE = x$seTE,
                   treat1 = x$treat1,
                   treat2 = x$treat2,
                   studlab = x$studlab,
                   sm = x$sm,
                   small.values = x$small.values,
                   keepdata = FALSE), pooled = pooled)$poth
    #
    if (verbose)
      setTxtProgressBar(pb, i)
  }
  #
  res <- list(dist = poths,
              d = d,
              pooled = pooled,
              netmeta = x)
  #
  class(res) <- c("refdist.poth", class(res))
  return(res)
}
