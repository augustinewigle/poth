#' Given a matrix of relative effects and standard errors, calculate the LOO SIR contribution
#' @param diffs Matrix of relative effects
#' @param ses Matrix of estimated standard errors for relative effects
#' @param trt_names optional; vector of treatment names matching order in diffs and sds
#' @param largerbetter Logical, T if larger outcomes are good. Defauls is T
#' @returns vector of LOO SIR contributions
#' @import ggplot2
#' @importFrom stringr str_wrap
#' @export
loo_sir <- function(diffs = NA,
                    ses = NA,
                    rank_mat = NA,
                    trt_names = NA,
                    largerbetter) {



  if(!is.null(diffs)) { # we are in Frequentist world

    n <- nrow(diffs)

    # name check
    if(length(trt_names) != n) {

      if(is.null(colnames(diffs))) {

        Warning("Please specify treatment names. Assigning generic treatment names")

        trt_names <- paste0("trt", 1:n)

        colnames(diffs) <- colnames(ses) <- trt_names

      }

      trt_names <- colnames(diffs)

    }

    # Calculate P-scores
    ps <- p_scores(diffs, ses, largerbetter, trt_names)

    # Regular SIR
    sir <-calc_sir(sucras = ps, treatment_names = trt_names)$sir

    # Calculate loo
    # treatment order
    id_order <- order(ps, decreasing = T)
    ordered_diffs <- diffs[id_order, id_order]
    ordered_ses <- ses[id_order, id_order]
    # loo_pscores[[1]] corresponds to the p-scores we calculate if we leave out the number one ranked treatment based on the base p-scores
    loo_pscores <- lapply(1:n, function(x) p_scores(ordered_diffs[-x,-x], ordered_ses[-x,-x], largerbetter, colnames(ordered_diffs)[-x]))

    names(loo_pscores) <- colnames(ordered_diffs)

    loo_sirs <- sapply(loo_pscores, function(x) calc_sir(sucras = x, treatment_names = names(x))$sir)
    names(loo_sirs) <- names(loo_pscores)

    sir_res <- sir-loo_sirs

    df <- data.frame(trt_name = names(sir_res),
                     sir_res = sir_res,
                     orig_pscore = ps[id_order],
                     rank = rank(-ps[id_order]))


    df$res_pos <- df$sir_res > 0
    df$mag_order <- order(abs(df$sir_res))

    df$xplot <- paste0(interaction(df$rank, df$trt_name, sep = "\n("), ")")

    df$xplot <- factor(df$xplot, levels = df$xplot[order(df$rank)])


    g <- ggplot(df, aes(x = xplot, y = sir_res, fill = res_pos, alpha = abs(sir_res))) +
      geom_col(col = "black") +
      geom_hline(yintercept = 0) +
      scale_fill_manual(breaks = c(F, T), values = c("red4", "darkgreen"),
                        labels = c("Reduces certainty", "Increases certainty")) +
      labs(y = expression(SIR-SIR^"-i"), x = "Global Rank (Treatment)",
           fill = str_wrap("Contribution to SIR", width = 20)) +
      guides(alpha = "none") +
      theme_bw()

    return(g)

  } else { # we are in Bayesian world

    print("implement this lol")

  }


}
