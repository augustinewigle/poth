#' poth: Brief overview of methods and general hints
#' 
#' R package \bold{poth} allows to calculate the precision of treatment hierarchy (POTH)
#' metric to quantify the uncertainty in a treatment hierarchy in network
#' meta-analysis (Wigle et al., 2024).
#'
#' @name poth-package
#'
#' @details
#' R package \bold{poth} provides the following methods:
#' \itemize{
#' \item Calculate the separation in ranking metric (\code{\link{poth}})
#' \item Conduct leave-one-out analysis (\code{\link{loo.poth}})
#' }
#' 
#' Type \code{help(package = "poth")} for a listing of R functions
#' available in \bold{poth}.
#' 
#' Type \code{citation("poth")} on how to cite \bold{poth} in
#' publications.
#' 
#' The development version of \bold{poth} is available on GitHub
#' \url{https://github.com/augustinewigle/poth}.
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
