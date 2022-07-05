#' Fits various tracks (different random starting values) and chooses best model
#'
#' Fits n_init tracks with different initial values and decides on best model
#' based on information criteria.
#'
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param q Number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param inf_crit Information criterion: `"BIC"` (default), `"AIC"` or
#'  `"AICc"`
#' @param task `"regression"` (default) or `"classification"`
#' @param unif Random initial values max value
#' @param maxit maximum number of iterations for nnet (default = 100)
#' @param ... additional argument for nnet
#' @return The best model from the different tracks
#' @export
nn_fit_tracks <- function(X, y, q, n_init, inf_crit = "BIC",
                          task = "regression", unif = 3, maxit = 1000, ...) {
  # Function with fits n_init tracks of model and finds best

  df <- data.frame(X, y)
  n <- nrow(X)
  p <- ncol(as.matrix(X)) # as.matrix() in case p = 1 (auto. becomes vector)

  k <- (p + 2) * q + 1

  weight_matrix_init <- matrix(stats::runif(n_init * k, min = -unif, max = unif), ncol = k)

  weight_matrix <- matrix(rep(NA, n_init * k), ncol = k)
  inf_crit_vec <- rep(NA, n_init)
  converge <- rep(NA, n_init)

  if (task == "regression") {
    linout <- TRUE
    entropy <- FALSE
  } else if (task == "classification") {
    linout <- FALSE
    entropy <- TRUE
  }

  for (iter in 1:n_init) {
    nn_model <- nnet::nnet(y ~ .,
                           data = df, size = q, trace = FALSE,
                           linout = linout, entropy = entropy,
                           Wts = weight_matrix_init[iter, ], maxit = maxit, ...
    )

    weight_matrix[iter, ] <- nn_model$wts

    if (task == "regression") {
      RSS <- nn_model$value
      sigma2 <- RSS / n

      log_likelihood <- (-n / 2) * log(2 * pi * sigma2) - RSS / (2 * sigma2)
    } else if (task == "classification") {
      log_likelihood <- -nn_model$value
    }

    inf_crit_vec[iter] <- ifelse(inf_crit == "AIC",
                                 (2 * (k + 1) - 2 * log_likelihood),
                                 ifelse(inf_crit == "BIC",
                                        (log(n) * (k + 1) - 2 * log_likelihood),
                                        ifelse(inf_crit == "AICc",
                                               (2 * (k + 1) * (n / (n - (k + 1) - 1)) - 2 * log_likelihood),
                                               NA
                                        )
                                 )
    )
    converge[iter] <- nn_model$convergence
  }
  W_opt <- weight_matrix[which.min(inf_crit_vec), ]

  return(list(
    "W_opt" = W_opt,
    "val" = min(inf_crit_vec),
    "inf_crit_vec" = inf_crit_vec,
    "converge" = converge,
    "weight_matrix" = weight_matrix
  ))
}


#' Neural network prediction
#'
#'
#' @param X Data
#' @param W Weight vector
#' @param q Number of hidden nodes
#' @param output Activation function for output unit: `"identity"` (default) or
#'  `"sigmoid"`
#' @return Prediction for given inputs
#' @export
nn_pred <- function(X, W, q, output = 'identity') {
  n <- nrow(X)
  p <- ncol(X)

  k <- (p + 2) * q + 1

  if (length(W) == k) {
    X <- cbind(rep(1, n), X)

    h_input <- X %*% t(matrix(W[1:((p + 1) * q)], nrow = q, byrow = TRUE))

    h_act <- cbind(rep(1, n), sigmoid(h_input))

    if (output == 'identity') {
      y_hat <- h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)
    } else if (output == 'sigmoid') {
      y_hat <- sigmoid(
        h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)
      )
    } else {
      stop(
        sprintf("Error: %s not recognised as available output function.",
                output))
    }

    return(y_hat)
  } else {
    stop(sprintf(
      "Error: Incorrect number of weights for NN structure. W should have
      %s weights (%s weights supplied).", k, length(W)))
  }
}
