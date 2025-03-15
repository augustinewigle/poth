#' Generate reference distribution for POTH for a given network structure
#' (without multi-arm trials)
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
#' from the network meta-analysis, i.e., \code{x$TE.common[, 1]} if
#' \code{pooled = "common"} and \code{x$TE.random[, 1]} if
#' \code{pooled = "random"}.
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
#'   data = Senn2013, subset = studlab != "Willms1999",
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
   if (pooled == "common")
     d <- x$TE.common[, 1]
   else
     d <- x$TE.random[, 1]
  }
  else
    chknumeric(d, length = length(x$trts))
  #
  chknumeric(nsim, min = 1, length = 1)
  chklogical(verbose)
  #
  if (any(x$multiarm))
    stop("Method not implemented for networks with multi-arms.",
         call. = FALSE)
  
  # Simulate data with the desired relative effects and identical structure
  # and heterogeneity
  #
  meanvec <- x$X.matrix %*% d
  #
  # Standard errors based on common or random effects model
  #
  if (pooled == "random")
    sdvec <- 1 / sqrt(x$w.random)
  else
    sdvec <- 1 / sqrt(x$w.common)
  
  simdata <- replicate(nsim,
                       rnorm(length(meanvec), mean = meanvec, sd = sdvec))
  
  # Calculate POTHs
  #
  poths <- numeric(nsim)
  #
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  #
  for (i in seq_len(nsim)) {
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
  poths
}
