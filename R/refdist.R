#' Generate reference distribution for given network structure (with no multi-arm trials)
#' @param net A netmeta object which has the desired structure including tau2
#' @param d a vector of the desired relative effects, must be in the same order as the treatment ordering in net$X.matrix
#' @param nsim number of times to simulate from reference distribution
#' @returns A vector of the POTH values
#' @export
refdist <- function(net, d, nsim) {

  # simulate data with desired relative effects and identical structure and heterogeneity
  meanvec <- net$X.matrix%*%d
  sdvec <- 1/sqrt(net$w.random) # assuming random effects desired - TODO make flexible

  simdata <- replicate(nsim, rnorm(length(meanvec), mean = meanvec, sd = sdvec))

  # calculate POTH
  poths <- numeric(nsim)
  for(i in 1:nsim) {

    newnet <-  netmeta(TE = simdata[,i],
                       seTE = net$seTE,
                       treat1 = net$treat1,
                       treat2 = net$treat2,
                       studlab = net$studlab,
                       sm = net$sm,
                       small.values = net$small.values,
                       keepdata = F)

    poths[i] <- poth(newnet)$poth
  }

  return(poths)


}
