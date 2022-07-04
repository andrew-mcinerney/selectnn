#' Selects hidden layer structure
#'
#' Performs either bulk or step-wise hidden node selection.
#'
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param Q Candidate number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param type Selection type: `"bulk"` (the default) or `"step"`
#' @param inf_crit Information criterion: `"BIC"` (default), `"AIC"` or
#'  `"AICc"`
#' @param unif Random initial values max value
#' @param ... additional argument for nnet
#' @return Optimal number of hidden nodes
#' @export
hidden_node_sel <- function(X, y, Q, n_init, type = "bulk", inf_crit = "BIC", unif = 3, ...) {
  if (type == "bulk") {
    inf_crit_vec <- rep(NA, Q)

    names(inf_crit_vec) <- as.character(1:Q)

    weights_min <- vector("list", length = Q)

    for (q in 1:Q) {
      nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, unif, ...)

      weights_min[[q]] <- nn$W_opt

      inf_crit_vec[q] <- nn$val
    }

    W_opt <- weights_min[[which.min(inf_crit_vec)]]
  } else if (type == "step") {
    n_candidates <- 3

    inf_crit_vec <- rep(NA, n_candidates)

    names(inf_crit_vec) <- as.character((Q - 1):(Q + 1))

    n <- nrow(X) # sample size
    p <- ncol(X) # number of covariates

    weights_min <- vector(mode = "list", length = n_candidates)
    names(weights_min) <- as.character((Q - 1):(Q + 1))


    for (q in (Q - 1):(Q + 1)) {
      if (q == 0) {
        weights_min[[as.character(q)]] <- NULL

        inf_crit_vec[as.character(q)] <- NA
      } else {
        k <- (p + 1) * q + (q + 1)


        nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, unif, ...)

        weights_min[[as.character(q)]] <- nn$W_opt

        inf_crit_vec[as.character(q)] <- nn$val
      }
    }

    W_opt <- weights_min[[names(which.min(inf_crit_vec))]]
  }

  return(list(
    "min" = as.numeric(names(which.min(inf_crit_vec))),
    "val" = min(inf_crit_vec, na.rm = TRUE),
    "inf_crit_vec" = inf_crit_vec,
    "weights_min" = weights_min,
    "W_opt" = W_opt
  ))
}


#' Selects input layer structure
#'
#' Performs either bulk or step-wise input node selection.
#'
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param q Number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param type Selection type: `"bulk"` (the default) or `"step"`
#' @param inf_crit Information criterion: `"BIC"` (default), `"AIC"` or
#'  `"AICc"`
#' @param unif Random initial values max value
#' @param ... additional argument for nnet
#' @return Inputs dropped from model
#' @export
input_node_sel <- function(X, y, q, n_init, type = "bulk", inf_crit = "BIC", unif = 3, ...) {
  full_model <- nn_fit_tracks(X, y, q, n_init, inf_crit, unif, ...)

  full_inf_crit <- full_model$val

  p <- ncol(X)

  X_full <- X

  colnames(X_full) <- 1:p

  colnames(X) <- 1:p

  dropped <- c()

  W_opt <- full_model$W_opt

  min_inf_crit <- full_inf_crit

  continue_drop <- TRUE

  while (continue_drop == TRUE) {
    nn_in <- input_importance(
      X = X, y = y, q = q, n_init = n_init,
      inf_crit = inf_crit, unif = unif,
      ...
    )

    if (nn_in$val < min_inf_crit) {
      X <- X[, -nn_in$min, drop = FALSE]
      p <- ncol(as.matrix(X))


      W_opt <- nn_in$W_opt

      dropped <- colnames(X_full)[!colnames(X_full) %in% colnames(X)]
      min_inf_crit <- nn_in$val

      if (ncol(X) == 1) {
        continue_drop <- FALSE
      }

      if (type == "step") {
        continue_drop <- FALSE
      }
    } else {
      continue_drop <- FALSE
    }
  }
  return(list(
    "X" = X, "p" = p, "W_opt" = W_opt, "dropped" = dropped,
    "full_inf_crit" = full_inf_crit, "inf_crit" = min_inf_crit
  ))
}

#' Determines the importance of each input
#'
#' Calculates the importance of each input model based on information criterion
#'  and returns which node is least important
#'
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param q Number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param inf_crit Information criterion: `"BIC"` (default), `"AIC"` or
#'  `"AICc"`
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

