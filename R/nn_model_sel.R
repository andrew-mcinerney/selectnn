
#' Full model selection (input and hidden layer)
#'
#' Performs model selection procedure.
#'
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param Q Candidate number of hidden nodes
#' @param n_init Number of random initialisations (tracks)`
#' @param inf_crit Information criterion: `"BIC"` (default), `"AIC"` or
#'  `"AICc"`
#' @param unif Random initial values max value
#' @param maxit maximum number of iterations for nnet (default = 100)
#' @param ... additional argument for nnet
#' @return Optimal number of hidden nodes
#' @export
nn_model_sel <- function(X, y, Q, n_init, inf_crit = "BIC", unif = 3,
                         maxit = 1000, ...) {
  continue <- 1

  hidden_size <- c()

  n_rep_h <- 1

  n_rep_i <- 1

  nn_hidden <- vector(mode = "list", length = 1)

  nn_input <- vector(mode = "list", length = 1)

  nn_hidden[[1]] <- hidden_node_sel(
    X = X,
    y = y,
    Q = Q,
    n_init = n_init,
    type = "bulk",
    inf_crit = inf_crit,
    unif = unif,
    maxit = maxit,
    ...
  )

  hidden_size[1] <- nn_hidden[[1]]$min

  nn_input[[1]] <- input_node_sel(
    X = X,
    y = y,
    q = nn_hidden[[1]]$min,
    n_init = n_init,
    type = "bulk",
    inf_crit = inf_crit,
    unif = unif,
    maxit = maxit,
    ...
  )

  if (!is.null(nn_input[[1]]$dropped)) {
    X_new <- X[, -as.numeric(nn_input[[1]]$dropped), drop = FALSE]
    dropped <- colnames(X)[as.numeric(nn_input[[1]]$dropped)]
  } else {
    X_new <- X
    dropped <- c()
  }


  while (continue == 1) {
    nn_hidden[[n_rep_h + 1]] <- hidden_node_sel(
      X = X_new,
      y = y,
      Q = nn_hidden[[n_rep_h]]$min,
      n_init = n_init,
      type = "step",
      inf_crit = inf_crit,
      unif = unif,
      maxit = maxit,
      ...
    )

    n_rep_h <- n_rep_h + 1

    hidden_size[n_rep_h] <- nn_hidden[[n_rep_h]]$min

    if (hidden_size[n_rep_h] == hidden_size[n_rep_h - 1]) {
      continue <- 0
    } else {
      nn_input[[n_rep_i + 1]] <- input_node_sel(
        X = X_new,
        y = y,
        q = nn_hidden[[n_rep_h]]$min,
        n_init = n_init,
        type = "step",
        inf_crit = inf_crit,
        unif = unif,
        maxit = maxit,
        ...
      )

      n_rep_i <- n_rep_i + 1

      if (!is.null(nn_input[[n_rep_i]]$dropped) & ncol(X_new) > 1) {
        dropped <- c(dropped, colnames(X_new)[as.numeric(nn_input[[n_rep_i]]$dropped)])
        X_new <- X_new[, -as.numeric(nn_input[[n_rep_i]]$dropped)]
      } else {
        continue <- 0
      }

      if (is.null(ncol(X_new))) {
        continue <- 0
      }
    }
  }


  return(list(
    "nn_hidden" = nn_hidden,
    "nn_input" = nn_input,
    "n_rep_h" = n_rep_h,
    "n_rep_i" = n_rep_i,
    "X" = X_new,
    "X_full" = X,
    "dropped" = dropped,
    "hidden_size" = hidden_size
  ))
}
