#' Calculate the measure for uncertainty for ranks
#' @param ranking_prob_mat square matrix where the rows represent ranks and the columns treatments.
#' @param sucras a vector of sucra values
#' @param treatment_names a vector of treatment names. Should match the order of treatments in ranking_prob_mat or sucras
#'
#' @returns The separation in ranking metric.
#' @export
calc_sir <- function(ranking_prob_mat = NA, sucras = NA, treatment_names) {

  if(length(ranking_prob_mat) > 1) {

    if(length(sucras)> 1) {

      warning("Calculating the SIR using the ranking probability matrix")

    }

    # error checking
    if(nrow(ranking_prob_mat) != ncol(ranking_prob_mat)) {

      stop("ranking_prob_mat must be a square matrix")

    } else if(any(apply(ranking_prob_mat, 1, sum) != 1)) {

      warning("The rows of ranking_prob_mat should sum to 1, check that it is a ranking matrix")

    } else if(any(apply(ranking_prob_mat, 2, sum) != 1)) {

      warning("The columns of ranking_prob_mat should sum to 1, check that it is a ranking matrix")

    }

    n <- ncol(ranking_prob_mat)

    rank_e <- apply(seq(1:n)*ranking_prob_mat, 2, sum)

    sucras <- (n - rank_e)/(n-1)


    rank_vars <- numeric(n)

    for(i in 1:n) {

      rank_vars[i] <- sum((seq(1:n) - rank_e[i])^2*ranking_prob_mat[,i])

    }

    names(sucras) <- treatment_names


    return(list(sucras = sort(sucras, decreasing = T),
                sir = 1-sum(rank_vars)*12/n/(n+1)/(n-1)))

  } else if(length(sucras) > 1) {

    # error checking
    if(abs(mean(sucras))-0.5 > 1e-7) {

      warning("The mean of the sucras should be 0.5. Check your values.")

    }

    n <- length(sucras)

    sir <- 3*((1-n)/(n+1) + 4*(n-1)/(n*(n+1))*sum(sucras^2))

    names(sucras) <- treatment_names

    return(list(sucras = sort(sucras, decreasing = T),
                sir = sir))

  } else {

    stop("One of ranking_prob_mat or sucras must be specified")

  }

}
