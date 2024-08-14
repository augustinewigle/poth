#' sir: Brief overview of methods and general hints
#' 
#' R package \bold{sir} allows to calculate the separation in ranking (SIR)
#' metric to quantify the uncertainty in a treatment hierarchy in network
#' meta-analysis (Wigle et al., 2024).
#'
#' @name sir-package
#'
#' @details
#' R package \bold{sir} provides the following methods:
#' \itemize{
#' \item Calculate the separation in ranking metric (\code{\link{sir}})
#' \item Conduct leave-one-out analysis (\code{\link{loo.sir}})
#' }
#' 
#' Type \code{help(package = "sir")} for a listing of R functions
#' available in \bold{sir}.
#' 
#' Type \code{citation("sir")} on how to cite \bold{sir} in
#' publications.
#' 
#' The development version of \bold{sir} is available on GitHub
#' \url{https://github.com/augustinewigle/sir}.
#' 
#' @author Augustine Wigle \email{amhwigle@uwaterloo.ca},
#'   Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @references
#' Wigle A, ... (2024):
#' Separation In Ranking: A Metric for Quantifying Uncertainty in Treatment
#' Hierarchies in Network Meta-Analysis
#'
#' @keywords package
#'
#' @import ggplot2
#' @importFrom stringr str_wrap str_trunc
#' @importFrom stats pnorm

"_PACKAGE"

NULL
