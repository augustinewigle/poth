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
#' @param object A \code{refdist} object.
#' @param \dots Additional arguments passed on to print or summary function.
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
  if (missing(d))
    d <- rep(0, x$n)
  else
    chknumeric(d, length = length(x$trts))
  #
  chknumeric(nsim, min = 1, length = 1)
  chklogical(verbose)
  
  # Simulate data with the desired relative effects and identical structure
  # and heterogeneity
  #
  meanvec <- x$X.matrix %*% d
  
  # Standard errors based on common or random effects model with multi-arm
  # corrections
  #
  if (pooled == "random")
    sdvec <- x$seTE.adj.random
  else
    sdvec <- x$seTE.adj.common
  
  
  # Calculate POTHs
  #
  poths <- vector("numeric", nsim)
  #
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  #
  for (i in seq_len(nsim)) {
    # Vector with treatment estimates
    TE.i <- rnorm(length(meanvec), mean = meanvec, sd = sdvec)
    #
    if (any(x$multiarm)) { # correct data to be consistent
      #
      # Loop over multi-arm studies
      #
      for (j in unique(x$studlab[x$multiarm])) {
        sel.j <- which(x$studlab == j)
        narm <- unique(x$n.arms[sel.j])
        # Indices for first ai - a contrasts for this study
        # (basic parameters)
        basic.j  <- sel.j[1:(narm - 1)]
        # Indices for remaining contrasts to be made internally consistent
        # (functional parameters)
        funct.j <- sel.j[narm:length(sel.j)]
        #
        A <- t(x$X.matrix[basic.j, ])
        B <- t(x$X.matrix[funct.j, , drop = FALSE])
        #
        TE.i[funct.j] <- matrix(TE.i[basic.j], nrow = 1) %*% MASS::ginv(A) %*% B
      }
    }
    #
    net.i <-  netmeta(TE = TE.i,
                      seTE = x$seTE,
                      treat1 = x$treat1,
                      treat2 = x$treat2,
                      studlab = x$studlab,
                      sm = x$sm,
                      small.values = x$small.values,
                      keepdata = FALSE)
    #
    poths[i] <- poth(net.i, pooled = pooled)$poth
    #
    if (verbose)
      setTxtProgressBar(pb, i)
  }
  #
  if (verbose)
    cat("\n")
  #
  res <- list(dist = poths,
              d = d,
              pooled = pooled,
              netmeta = x)
  #
  class(res) <- c("refdist", class(res))
  #
  res
}


#' @rdname refdist
#' @method print refdist
#' @export

print.refdist <- function(x, ...) {
  chkclass(x, "refdist")
  #
  print(x$dist, ...)
  #
  invisible(NULL)
}


#' @rdname refdist
#' @method summary refdist
#' @export

summary.refdist <- function(object, ...) {
  chkclass(object, "refdist")
  #
  summary(object$dist, ...)
}
