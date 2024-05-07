#' Calculate the local SIR for a subset of treatments
#'
#' @param diffs Matrix of relative effects, must have column names that are the name of each treatment or trts must be specified
#' @param ses Matrix of estimated standard errors for relative effects, must have column names that are the name of each treatment or trts must be specified
#' @param trts optional; vector of treatment names matching order in diffs and sds
#' @param largerbetter Logical, T if larger outcomes are good. Defauls is T
#' @param samples a matrix or data.frame of MCMC samples, where rows are MCMC samples and columns are relative effects (relative to anchor) for treatments.
#' must have column names that are the name of each treatment or trts must be specified
#' @param subset A character vector of treatment names to consider as the set of competing treatments
#' @param topx Calculate the local SIR for the top x treatments, as ranked by global SUCRA/P-scores
#' @param bottomx Calculate the local SIR for the bottom x treatments, as ranked by global SUCRA/P-scores
#' @param do_loo Should a leave-one-out analysis of the subset be performed? Default is F
#' 
#' @return A named list with two elements:
#' df - a data frame with the information used to draw the plot
#' plot - a ggplot object
#' 
#' @import ggplot2
#' @importFrom stringr str_wrap
#' 
#' @export

local_sir <- function(diffs = NA,
                      ses = NA,
                      samples = NA,
                      trts = NULL,
                      largerbetter,
                      subset = NA,
                      topx = NA,
                      bottomx = NA,
                      do_loo = F) {

  #
  # Checks, define score type, check treatment names
  #
  
  if (!is.null(dim(diffs))) {
    bayes <- FALSE # we are in frequentist world
    
    if (is.null(dim(ses)))
      stop("Argument 'ses' mandatory if argument 'diffs' is provided.")
    
    score_type <- "P-score"
    n <- nrow(diffs)
    
    # name check
    if (length(trts) != n) {
      if (is.null(colnames(diffs)))
        stop("No treatment names specified. ",
          "Treatment names are required for local_sir either through the ",
          "trts argument or the column names of diffs.")
      #
      trts <- colnames(diffs)
    }
    else
      colnames(diffs) <- colnames(ses) <- trts
  }
  else {
    # Bayesian
    bayes <- TRUE
    
    n <- ncol(samples)
    score_type <- "SUCRA"
    
    # name check
    if (length(trts) != n) {
      if (is.null(colnames(samples)))
        stop("No treatment names specified. ",
             "Treatment names are required for local_sir either through the ",
             "trts argument or the column names of diffs.")
      #
      trts <- colnames(samples)
    }
    else
      colnames(samples) <- trts
  }
  
  
  #
  # Define subsets
  #
  
  if (length(subset) == 1) {
    if (is.na(subset)) { # subset is NA
      if (is.na(topx)&is.na(bottomx)) {
        stop("One of subset, topx, or bottomx must be specified")
      }
      else {
        # Calculate ranks
        if (bayes) {
          rmat <- rank_mat(samples = samples, largerbetter = largerbetter, trts = trts)
          scores <- sir(rmat, trts = trts)$sucras
        }
        else
          scores <- pscores(diffs, ses, largerbetter, trts)
        
        trt_order <- names(scores)[order(scores, decreasing = TRUE)]
        
        if (!is.na(topx)) { # define subset as the top x treatments
          if (topx > n - 1)
            stop('topx should be less than the number of treatments')
          #
          subset <- trt_order[1:topx]
        }
        else if (!is.na(bottomx)) { # define subset as the top x treatments
          if (bottomx > n - 1) {
            stop('bottomx should be less than the number of treatments')
          }
          #
          subset <- trt_order[(n - bottomx + 1):n]
        }
      }
    }
  }
  
  
  #
  # Check that treatment names exist
  #
  
  if (any(!(subset %in% trts))) {
    wrong_names <- subset[which(!(subset %in% trts))]
    #
    stop("The following treatments in subset were not found in the ",
         "treatment names:\n    ",
         paste(wrong_names, collapse = ",\n    "),
         "\n  All treatments in subset must be found in the treatment names")
  }
  
  
  #
  # Calculate scores and SIR for the subset
  #
  
  if (bayes) {
    # Calculate SUCRAs and global SIR
    rps <- rank_mat(samples[,subset], largerbetter, trts)
    both <- sir(rps, trts = colnames(rps))
    scores <- both$sucras
    sir <- both$sir
    
    if (do_loo)
      loo <- loo_sir(samples = samples[,subset],
                     largerbetter = largerbetter,
                     trts = subset)
    else
      loo <- NULL
  }
  else {
    # Calculate P-scores
    scores <- pscores(diffs[subset, subset], ses[subset, subset],
                      largerbetter, trts = subset)
    # Regular SIR
    sir <- sir(scores, trts = subset)$sir
    #
    if (do_loo)
      loo <- loo_sir(diffs = diffs[subset, subset],
                     ses = ses[subset, subset],
                     largerbetter = largerbetter,
                     trts = subset)
    else
      loo <- NULL
  }
  
  
  list(scores = scores,
       score_type = score_type,
       sir = sir,
       loo = loo)
}
