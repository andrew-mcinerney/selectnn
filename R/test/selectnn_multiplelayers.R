selectnn1 <- function(X, y, Q, n_init, inf_crit = "BIC",
                             task = "regression", unif = 3, maxit = 1000,
                             ...) {
  if (any(is.na(X)) | any(is.na(y))) {
    stop("Error: Make sure data does not contain any NAs.")
  }

  if (!(inf_crit %in% c("BIC", "AIC", "AICc"))) {
    stop(sprintf(
      "Error: %s not recognised as information criterion.
      Please choose from AIC, AICc or BIC.", inf_crit
    ))
  }

  if (!(task %in% c("regression", "classification"))) {
    stop(sprintf(
      "Error: %s not recognised as task. Please choose from regression or classification.",
      task
    ))
  }

  cl <- match.call()


  if (is.null(colnames(X))) {
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = deparse(substitute(X)))
  }

  continue <- 1

  n_layers <- c()

  hidden_size <- c()

  n_rep_l <- 1

  n_rep_h <- 1

  n_rep_i <- 1

  nn_layers <- vector(mode = "list", length = 1)

  nn_hidden <- vector(mode = "list", length = 1)

  nn_input <- vector(mode = "list", length = 1)

  nn_layers[[1]] <- hidden_layer_sel(
    X = X,
    y = y,
    Q = Q,
    n_init = n_init,
    inf_crit = inf_crit,
    task = task,
    unif = unif,
    maxit = maxit,
    ...
  )

  nn_hidden[[1]] <- hidden_node_sel(
    X = X,
    y = y,
    Q = Q,
    n_init = n_init,
    type = "bulk",
    inf_crit = inf_crit,
    task = task,
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
    task = task,
    unif = unif,
    maxit = maxit,
    ...
  )

  if (!is.null(nn_input[[1]]$dropped)) {
    X_new <- X[, !colnames(X) %in% nn_input[[1]]$dropped, drop = FALSE]
    dropped <- nn_input[[1]]$dropped
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
      task = task,
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
        task = task,
        unif = unif,
        X_full = X,
        maxit = maxit,
        ...
      )

      n_rep_i <- n_rep_i + 1

      if (!is.null(nn_input[[n_rep_i]]$dropped) & ncol(X_new) > 1) {
        X_new <- nn_input[[n_rep_i]]$X
        dropped <- colnames(X)[!colnames(X) %in% colnames(X_new)]
      } else {
        continue <- 0
      }

      if (is.null(ncol(X_new))) {
        continue <- 0
      }
    }
  }

  delta_bic <- nn_input[[n_rep_i]]$delta_bic

  if (n_rep_i == 1) {
    nn_input_extra <- input_node_sel(
      X = X_new,
      y = y,
      q = nn_hidden[[n_rep_h]]$min,
      n_init = n_init,
      type = "step",
      inf_crit = inf_crit,
      task = task,
      unif = unif,
      X_full = X,
      maxit = maxit,
      ...
    )

    delta_bic <- nn_input_extra$delta_bic
  }

  p <- ncol(X_new)
  q <- hidden_size[n_rep_h]
  p_init <- ncol(X)
  q_init <- Q

  if (n_rep_h > n_rep_i) {
    W_opt <- nn_hidden[[n_rep_h]]$W_opt
    value <- nn_hidden[[n_rep_h]]$value
  } else {
    W_opt <- nn_input[[n_rep_i]]$W_opt
    value <- nn_input[[n_rep_i]]$value
  }


  net <- list(
    "W_opt" = W_opt,
    "value" = value,
    "nn_hidden" = nn_hidden,
    "nn_input" = nn_input,
    "n_rep_h" = n_rep_h,
    "n_rep_i" = n_rep_i,
    "X" = X_new,
    "X_full" = X,
    "dropped" = dropped,
    "hidden_size" = hidden_size,
    "p" = p,
    "p_init" = p_init,
    "q" = q,
    "q_init" = q_init,
    "delta_bic" = delta_bic,
    "call" = cl
  )

  class(net) <- "selectnn"


  return(net)
}
