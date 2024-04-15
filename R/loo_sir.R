#' Given a matrix of relative effects and standard errors, calculate the LOO SIR contribution
#'
#' @param diffs Matrix of relative effects, must have column names that are the name of each treatment or trt_names must be specified
#' @param ses Matrix of estimated standard errors for relative effects, must have column names that are the name of each treatment or trt_names must be specified
#' @param trt_names optional; vector of treatment names matching order in diffs and sds
#' @param largerbetter Logical, T if larger outcomes are good. Defauls is T
#' @param samples a matrix or data.frame of MCMC samples, where rows are MCMC samples and columns are relative effects (relative to anchor) for treatments.
#' must have column names that are the name of each treatment or trt_names must be specified
#' @returns A named list with two elements:
#' df - a dataframe with the information used to draw the plot
#' plot - a ggplot object
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @export
loo_sir <- function(diffs = NA,
                    ses = NA,
                    samples = NA,
                    trt_names = NA,
                    largerbetter) {

  if(!is.null(dim(diffs))) { # we are in Frequentist world

    if(is.null(dim(ses))) {

      stop("If diffs is specified then ses must also be specified")

    }

    score_type <- "P-score"
    n <- nrow(diffs)

    # name check
    if(length(trt_names) != n) {

      if(is.null(colnames(diffs))) {

        warning("No treatment names specified. Assigning generic treatment names")

        trt_names <- paste0("trt", 1:n)

        colnames(diffs) <- colnames(ses) <- trt_names

      }

      trt_names <- colnames(diffs)

    } else {

      colnames(diffs) <- colnames(ses) <- trt_names

    }

    # Calculate P-scores
    scores <- p_scores(diffs, ses, largerbetter, trt_names)

    # Regular SIR
    sir <-calc_sir(sucras = scores, trt_names = trt_names)$sir

    # Calculate loo
    # treatment order
    id_order <- order(scores, decreasing = T)
    ordered_diffs <- diffs[id_order, id_order]
    ordered_ses <- ses[id_order, id_order]
    # loo_pscores[[1]] corresponds to the p-scores we calculate if we leave out the number one ranked treatment based on the base p-scores
    loo_pscores <- lapply(1:n, function(x) p_scores(ordered_diffs[-x,-x], ordered_ses[-x,-x], largerbetter, colnames(ordered_diffs)[-x]))
    names(loo_pscores) <- colnames(ordered_diffs)

    loo_sirs <- sapply(loo_pscores, function(x) calc_sir(sucras = x, trt_names = names(x))$sir)
    names(loo_sirs) <- names(loo_pscores)

  } else { # we are in Bayesian world

    n <- ncol(samples)
    score_type <- "SUCRA"

    # name check
    if(length(trt_names) != n) {

      if(is.null(colnames(samples))) {

        warning("No treatment names specified. Assigning generic treatment names")

        trt_names <- paste0("trt", 1:n)

        colnames(samples) <- trt_names

      }

      trt_names <- colnames(samples)

    } else {

      colnames(samples) <- trt_names

    }

    # Calculate SUCRAs and global SIR
    rps <- get_rank_prob_mat(samples, largerbetter, trt_names)
    both <- calc_sir(ranking_prob_mat = rps, trt_names = trt_names)
    scores <- both$sucras[trt_names]
    sir <- both$sir

    id_order <- order(scores, decreasing = F)
    ordered_samples <- samples[,id_order]

    # LOO
    # recalculate rp matrix without i
    # recalculate sir

    loo_rps <- lapply(1:n, function(x) get_rank_prob_mat(ordered_samples[,-x], largerbetter, trt_names = colnames(ordered_samples[-x])))
    loo_sirs <- sapply(loo_rps, function(x) calc_sir(ranking_prob_mat = x)$sir)
    names(loo_sirs) <- colnames(ordered_samples)

  }



  sir_res <- sir-loo_sirs

  df <- data.frame(trt_name = names(sir_res),
                   sir_res = sir_res,
                   orig_score = scores[id_order],
                   score_type = score_type,
                   rank = rank(-scores[id_order]),
                   sir_res_ratio = sir_res/sir,
                   SIR= sir)

  df$res_pos <- df$sir_res > 0
  df$mag_order <- order(abs(df$sir_res))

  df$xplot <- paste0(interaction(df$rank, df$trt_name, sep = "\n("), ")")

  df$xplot <- factor(df$xplot, levels = df$xplot[order(df$rank)])

  g <- ggplot(df, aes(x = xplot, y = sir_res, fill = res_pos, alpha = abs(sir_res))) +
    geom_col(col = "black") +
    coord_cartesian(ylim = c(-max(abs(df$sir_res)), max(abs(df$sir_res)))) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(breaks = c(F, T), values = c("red4", "darkgreen"),
                      labels = c("Reduces certainty", "Increases certainty")) +
    labs(y = ifelse(F, # can change if I want to go back to the ratio thing
                    expression(frac(SIR-SIR^"-i",SIR)),
                    expression(SIR-SIR^"-i")),
         x = paste0("Global Rank Using ", score_type, "s (Treatment)"),
         fill = str_wrap(paste0("Contribution to SIR\n(Global SIR ", ifelse(sir < 0.001, "< 0.001", paste0("= ", round(sir, digits = 3))), ")"), width = 20)) +
    guides(alpha = "none") +
    theme_bw()

  return(list(plot = g,
              info = df))


}
