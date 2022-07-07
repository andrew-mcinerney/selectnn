
#' Neural network model selection
#'
#' Performs both input and hidden layer selection for neural networks.
#'
#' @return A list with information of the optimal model.
#' \itemize{
#'   \item \code{nn_hidden} - list of hidden node selection results.
#'   \item \code{nn_input} - list of input node selection results.
#'   \item \code{n_rep_h} - number of hidden node selection steps.
#'   \item \code{n_rep_i} - number of input node selection steps.
#'   \item \code{X} - matrix of the important covariates found.
#'   \item \code{X_full} - matrix of all covariates.
#'   \item \code{dropped} - vector of unimportant covariates.
#'   \item \code{hidden_size} - vector of hidden layer size found at each step.
#'   }
#'
#' @export
nn_select <- function(...) UseMethod("nn_select")

#' @rdname nn_select
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param Q Candidate number of hidden nodes
#' @param n_init Number of random initialisations (tracks)`
#' @param inf_crit Information criterion: `"BIC"` (default), `"AIC"` or
#'  `"AICc"`
#' @param task `"regression"` (default) or `"classification"`
#' @param unif Random initial values max value
#' @param maxit maximum number of iterations for nnet (default = 100)
#' @param ... arguments passed to or from other methods
#' @export
nn_select.default <- function(X, y, Q, n_init, inf_crit = "BIC",
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

  p <- ncol(X_new)
  q <- hidden_size[n_rep_h]
  p_init <- ncol(X)
  q_init <- Q

  net <- list(
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
    "call" = cl
  )

  class(net) <- "nn_select"


  return(net)
}

#' @rdname nn_select
#' @param formula A formula of the form: response ~ x1 + x2 + ...
#' @param data Data frame from which variables specified in formula are to be
#'   taken
#' @export
nn_select.formula <- function(formula, data, ...) {

  cl <- match.call()

  mf <- stats::model.frame(formula, data = data)

  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(formula, data = data)[, -1]

  rownames(X) <- NULL

  nn <- nn_select.default(X, y, ...)

  nn$call <- cl

  return(nn)
}


#' @export
print.nn_select <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Model Selected: ", x$p, "-", x$q, "-", "1", " network",
      sep="")
  cat(" with", (x$p + 2) * x$q + 1,"weights\n")
  cat("Initial Model: ", x$p_init, "-", x$q_init, "-", "1", " network",
      sep="")
  cat(" with", (x$p_init + 2) * x$q_init + 1,"weights\n")
}
