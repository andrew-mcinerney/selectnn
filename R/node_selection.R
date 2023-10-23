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
#' @param task `"regression"` (default) or `"classification"`
#' @param unif Random initial values max value
#' @param maxit maximum number of iterations for nnet (default = 100)
#' @param ... additional argument for nnet
#' @return Optimal number of hidden nodes
#' @export
hidden_node_sel <- function(X, y, Q, n_init, type = "bulk", inf_crit = "BIC",
                            task = "regression", unif = 3, maxit = 1000, ...) {
  if (type == "bulk") {
    inf_crit_vec <- rep(NA, Q)

    names(inf_crit_vec) <- as.character(1:Q)

    weights_min <- vector("list", length = Q)

    for (q in 1:Q) {
      nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, task, unif,
        maxit = maxit, ...
      )

      weights_min[[q]] <- nn$W_opt

      inf_crit_vec[q] <- nn$value
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


        nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, task, unif,
          maxit = maxit, ...
        )

        weights_min[[as.character(q)]] <- nn$W_opt

        inf_crit_vec[as.character(q)] <- nn$value
      }
    }

    W_opt <- weights_min[[names(which.min(inf_crit_vec))]]
  }

  return(list(
    "min" = as.numeric(names(which.min(inf_crit_vec))),
    "value" = min(inf_crit_vec, na.rm = TRUE),
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
#' @param task `"regression"` (default) or `"classification"`
#' @param unif Random initial values max value
#' @param X_full Full matrix of covariates if X has some dropped
#' @param maxit maximum number of iterations for nnet (default = 100)
#' @param ... additional argument for nnet
#' @return Inputs dropped from model
#' @export
input_node_sel <- function(X, y, q, n_init, type = "bulk", inf_crit = "BIC",
                           task = "regression", unif = 3, X_full = NULL,
                           maxit = 1000, ...) {
  init_model <- nn_fit_tracks(X, y, q, n_init, inf_crit, task, unif,
    maxit = maxit, ...
  )

  init_inf_crit <- init_model$value

  p <- ncol(X)

  X_init <- X

  delta_bic <- NULL

  dropped <- c()

  added <- c()

  W_opt <- init_model$W_opt

  min_inf_crit <- init_inf_crit

  if (type == "bulk") {
    continue_drop <- TRUE

    while (continue_drop == TRUE) {
      nn_in <- input_importance(
        X = X, y = y, q = q, n_init = n_init,
        inf_crit = inf_crit, task = task, unif = unif, maxit = maxit,
        ...
      )

      if (nn_in$value < min_inf_crit) {
        X <- X[, -nn_in$min, drop = FALSE]
        p <- ncol(as.matrix(X))


        W_opt <- nn_in$W_opt

        dropped <- colnames(X_init)[!colnames(X_init) %in% colnames(X)]
        min_inf_crit <- nn_in$value

        if (ncol(X) == 1) {
          continue_drop <- FALSE
        }
      } else {
        continue_drop <- FALSE
      }
    }
  } else if (type == "step") {
    if (is.null(X_full)) {
      stop("Error: X_full must not be null when performing stepwise")
    }


    nn_in <- input_importance(
      X = X, y = y, q = q, n_init = n_init,
      inf_crit = inf_crit, task = task, unif = unif,
      addition = TRUE, X_full = X_full, maxit = maxit,
      ...
    )

    if (nn_in$value < min_inf_crit) {
      if (names(nn_in$min) %in% colnames(X)) {
        X <- X[, -nn_in$min, drop = FALSE]
      } else {
        X <- cbind(X, X_full[, names(nn_in$min), drop = FALSE])
      }


      p <- ncol(as.matrix(X))


      W_opt <- nn_in$W_opt

      dropped <- colnames(X_init)[!colnames(X_init) %in% colnames(X)]
      added <- colnames(X)[!colnames(X) %in% colnames(X_init)]
      min_inf_crit <- nn_in$value
    }

    delta_bic <- nn_in$inf_crit_vec - init_inf_crit
  }
  return(list(
    "X" = X, "p" = p, "W_opt" = W_opt, "dropped" = dropped, "added" = added,
    "init_value" = init_inf_crit, "value" = min_inf_crit, "delta_bic" = delta_bic
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
#' @param task `"regression"` (default) or `"classification"`
#' @param unif Random initial values max value
#' @param addition Switch for addition step (default FALSE)
#' @param X_full Full matrix of covariates if X has some dropped
#' @param maxit maximum number of iterations for nnet (default = 100)
#' @param ... additional argument for nnet
#' @return The least important input node
#' @export
input_importance <- function(X, y, q, n_init, inf_crit = "BIC",
                             task = "regression", unif = 3,
                             addition = FALSE, X_full = NULL, maxit = 1000, ...) {
  if (addition == TRUE & is.null(X_full)) {
    stop("Error: X_full must not be null to allow addition step")
  }

  if (addition == TRUE & (is.null(colnames(X_full)) | is.null(colnames(X)))) {
    stop("Error: Column names must be supplied for X and X_full")
  }

  p_init <- ncol(X)

  inf_crit_vec <- rep(NA, p_init)

  weights_min <- vector("list", length = p_init)

  if (p_init > 1) {

    for (p in 1:p_init) {
      X_new <- X[, -p, drop = FALSE]

      nn <- nn_fit_tracks(X_new, y, q, n_init, inf_crit, task, unif,
                          maxit = maxit
      )

      weights_min[[p]] <- nn$W_opt

      inf_crit_vec[p] <- nn$value

      names(inf_crit_vec)[p] <- colnames(X)[p]
    }
  }

  if (addition == TRUE) {
    dropped <- colnames(X_full)[!colnames(X_full) %in% colnames(X)]

    i <- 1
    for (p in dropped) {
      X_new <- cbind(X, X_full[, p])

      nn <- nn_fit_tracks(X_new, y, q, n_init, inf_crit, task, unif,
                          maxit = maxit
      )

      weights_min[[p_init + i]] <- nn$W_opt

      inf_crit_vec[p_init + i] <- nn$value

      names(inf_crit_vec)[p_init + i] <- p

      i <- i + 1
    }
  }

  W_opt <- weights_min[[which.min(inf_crit_vec)]]

  return(list(
    "min" = which.min(inf_crit_vec),
    "value" = min(inf_crit_vec, na.rm = TRUE),
    "inf_crit_vec" = inf_crit_vec,
    "weights_min" = weights_min,
    "W_opt" = W_opt
  ))
}
