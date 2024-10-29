#'  Safety of Immune Checkpoint Inhibitors in Cancer
#'
#' @description
#' Data from a network meta-analysis on immune checkpoint inhibitors (ICIs)
#' to assess the safety of ICI drugs as cancer treatment (Xu et al., 2018).
#'
#' @details
#' Data were obtained from Rosenberger et al. (2021), who re-analysed the data.
#' There are seven treatments and 23 studies. The outcome of interest is the
#' number of treatment-related adverse events, so smaller values indicate a
#' safer treatment.
#'
#' @name Xu2018
#'
#' @docType data
#'
#' @format
#' A data frame with the following columns:
#' \tabular{rl}{
#' \bold{\emph{studyID}}\tab study id \cr
#' \bold{\emph{treatment}}\tab treatment name \cr
#' \bold{\emph{adverse}}\tab number of adverse events \cr
#' \bold{\emph{n}}\tab group sample size
#' }
#'
#' @source
#' Rosenberger, K.J., Duan, R., Chen, Y. et al. (2021):
#' Predictive P-score for treatment ranking in Bayesian network meta-analysis.
#' \emph{BMC Med Res Methodol} \bold{21}, 213.
#' \doi{10.1186/s12874-021-01397-5}
#'
#' @keywords datasets
#'
#' @seealso \code{\link[netmeta]{netmeta}}
#'
#' @examples
#' data(Xu2018)
#' head(Xu2018)
#'
#' library("netmeta")
#' pw <- pairwise(treat = treatment, event = adverse, n = n,
#'   studlab = studyID, data = Xu2018, sm = "OR")
#' #
#' net <- netmeta(pw, small.values = "desirable", method.tau = "REML",
#'   common = FALSE)
#' #
#' poth(net)

NULL
