#' Determines the importance of each input
#'
#' Calculates the importance of each input model based on information criterion
#'  and returns which node is least important
#'
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param q Number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param inf_crit Information criterion
#' @param unif Random initial values max value
#' @param ... additional argument for nnet
#' @return The least important input node
#' @export
input_importance <- function(X, y, q, n_init, inf_crit = "BIC", unif = 3, ...) {
  p_full <- ncol(X)

  inf_crit_vec <- rep(NA, p_full)

  weights_min <- vector("list", length = p_full)

  for (p in 1:p_full) {
    X_new <- X[, -p, drop = FALSE]

    nn <- nn_fit_tracks(X_new, y, q, n_init, inf_crit, unif, ...)

    weights_min[[p]] <- nn$W_opt

    inf_crit_vec[p] <- nn$val
  }

  W_opt <- weights_min[[which.min(inf_crit_vec)]]

  return(list(
    "min" = which.min(inf_crit_vec),
    "val" = min(inf_crit_vec),
    "inf_crit_vec" = inf_crit_vec,
    "weights_min" = weights_min,
    "W_opt" = W_opt
  ))
}
