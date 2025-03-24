#' poth: Brief overview of methods and general hints
#'
#' R package \bold{poth} allows to calculate the precision of treatment
#' hierarchy (POTH) metric to quantify the uncertainty in a treatment hierarchy
#' in network meta-analysis (Wigle et al., 2024).
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
#' Wigle A, Béliveau A, Salanti G, Rücker G, Schwarzer G, Mavridis D,
#' Nikolakopoulou A (2024):
#' Precision of treatment hierarchy: A metric for quantifying uncertainty in
#' treatment hierarchies in network meta-analysis.
#' Preprint on arXiv, \url{https://doi.org/10.48550/arXiv.2501.11596}
#'
#' @keywords package
#'
#' @import ggplot2
#' @importFrom stringr str_wrap str_trunc
#' @importFrom stats pnorm rnorm
#' @importFrom MASS ginv
#' @importFrom netmeta netmeta treats
#' @importFrom utils setTxtProgressBar txtProgressBar

"_PACKAGE"

NULL
