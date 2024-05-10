#' Given a matrix of relative effects and standard errors, calculate the LOO SIR contribution
#'
#' @param diffs Matrix of relative effects, must have column names that are the name of each treatment or trts must be specified
#' @param ses Matrix of estimated standard errors for relative effects, must have column names that are the name of each treatment or trts must be specified
#' @param trts optional; vector of treatment names matching order in diffs and sds
#' @param largerbetter Logical, T if larger outcomes are good. Defauls is T
#' @param samples a matrix or data.frame of MCMC samples, where rows are MCMC samples and columns are relative effects (relative to anchor) for treatments.
#' must have column names that are the name of each treatment or trts must be specified
#' @param plot_trts Logical; should treatment names be included in the plot? Default it TRUE
#' 
#' @return A named list with two elements:
#' df - a dataframe with the information used to draw the plot
#' plot - a ggplot object
#' 
#' @import ggplot2
#' @importFrom stringr str_trunc
#' 
#' @export

loo_sir <- function(diffs = NA,
                    ses = NA,
                    samples = NA,
                    trts = NULL,
                    largerbetter,
                    plot_trts = TRUE) {
  
  #
  # Checks and definitions
  #
  
  if (!is.null(dim(diffs))) {
    bayes <- FALSE # we are in frequentist world
    
    if (is.null(dim(ses)))
      stop("Argument 'ses' mandatory if argument 'diffs' is provided.")
    
    score_type <- "P-score"
    n <- nrow(diffs)
    
    #
    # name check
    #
    
    if (length(trts) != n) {
      if (is.null(colnames(diffs))) {
        warning("Assigning generic treatment names.")
        #
        trts <- paste0("trt", 1:n)
        colnames(diffs) <- colnames(ses) <- trts
      }
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
    
    #
    # name check
    #
    
    if (length(trts) != n) {
      
      if (is.null(colnames(samples))) {
        warning("Assigning generic treatment names.")
        #
        trts <- paste0("trt", 1:n)
        #
        colnames(samples) <- trts
      }
      #
      trts <- colnames(samples)
    }
    else
      colnames(samples) <- trts
  }
  
  #
  # Calculate scores and do LOO
  #
  
  if (bayes) {
    # Calculate SUCRAs and global SIR
    rps <- rank_mat(samples, largerbetter, trts)
    both <- sir(rps, trts = trts)
    scores <- both$sucras[trts]
    sir <- both$sir
    #
    id_order <- order(scores, decreasing = F)
    ordered_samples <- samples[,id_order]
    # LOO
    # recalculate rp matrix without i
    # recalculate sir
    loo_rps <- lapply(1:n, function(x) rank_mat(ordered_samples[,-x], largerbetter, trts = colnames(ordered_samples[-x])))
    loo_sirs <- sapply(loo_rps, function(x) sir(x)$sir)
    names(loo_sirs) <- colnames(ordered_samples)
  }
  else {
    # Calculate P-scores
    scores <- pscores(diffs, ses, largerbetter, trts)
    # Regular SIR
    sir <- sir(scores, trts = trts)$sir
    # Calculate loo
    # treatment order
    id_order <- order(scores, decreasing = TRUE)
    ordered_diffs <- diffs[id_order, id_order]
    ordered_ses <- ses[id_order, id_order]
    # loo_pscores[[1]] corresponds to the p-scores we calculate if we
    # leave out the number one ranked treatment based on the base p-scores
    loo_pscores <-
      lapply(1:n,
             function(x)
               pscores(ordered_diffs[-x, -x],
                        ordered_ses[-x,-x], largerbetter,
                        colnames(ordered_diffs)[-x]))
    names(loo_pscores) <- colnames(ordered_diffs)
    #
    loo_sirs <- sapply(loo_pscores, function(x) sir(x, trts = names(x))$sir)
    names(loo_sirs) <- names(loo_pscores)
  }
  
  
  #
  # Put everything together
  #
  sir_res <- sir - loo_sirs
  #
  df <- data.frame(trt_name = names(sir_res),
                   sir_res = sir_res,
                   orig_score = scores[id_order],
                   score_type = score_type,
                   rank = rank(-scores[id_order]),
                   sir_res_ratio = sir_res/sum(abs(sir_res)),
                   SIR= sir)
  #
  df$res_pos <- df$sir_res > 0
  #
  df <- df[order(df$rank), ]
  # # df$mag_order <- order(abs(df$sir_res))
  #
  # df$xplot <- paste0(interaction(df$rank, df$trt_name, sep = "\n("), ")")
  #
  # df$xplot <- factor(df$xplot, levels = df$xplot[order(df$rank)])
  #
  # # df$short_names <- unlist(lapply(as.list(df$trt_name),
  # #                          function(x) ifelse(strwidth(x) > 10, paste0(substr(x, start = 1, stop = 9), "."), x)))
  
  # Get rid of warning 'Undefined global functions or variables'
  res_pos <- NULL
  
  g <- ggplot(df,
              aes(x = 1:n, y = sir_res, fill = res_pos,alpha = abs(sir_res))) +
    geom_col(col = "black") +
    coord_cartesian(ylim = c(-max(abs(df$sir_res)), max(abs(df$sir_res)))) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(breaks = c(F, T), values = c("red4", "darkgreen"),
                      labels = c("Reduces certainty", "Increases certainty")) +
    labs(y = expression(SIR-SIR^"-i"),
         x = paste0("Global Rank Using ", score_type, "s"),
         fill = str_wrap(paste0("Contribution to SIR\n(Global SIR ",
                                ifelse(sir < 0.001, "< 0.001",
                                       paste0("= ", round(sir, digits = 3))),
                                ")"), width = 20)) +
    guides(alpha = "none") +
    theme_bw()
  
  if (plot_trts) {
    g <- g + scale_x_continuous(
      sec.axis = sec_axis(~.,
                          breaks = 1:n,
                          labels = str_trunc(df$trt_name, width = 13),
                          name = "Treatment Name"),
      breaks = 1:n, minor_breaks = NULL) +
      theme(axis.text.x.top =
              element_text(angle = 90, hjust = 0.01, vjust = 0.4),
            axis.title.x.top = element_blank())
  }
  else
    g <- g + scale_x_continuous(breaks = 1:n, minor_breaks = NULL)
  
  list(plot = g, info = df)
}
