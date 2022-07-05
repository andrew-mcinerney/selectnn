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
#' @param unif Random initial values max value
#' @param maxit maximum number of iterations for nnet (default = 100)
#' @param ... additional argument for nnet
#' @return The best model from the different tracks
#' @export
nn_fit_tracks <- function(X, y, q, n_init, inf_crit = "BIC", unif = 3,
                          maxit = 1000, ...) {
  # Function with fits n_init tracks of model and finds best

  df <- data.frame(X, y)
  n <- nrow(X)
  p <- ncol(as.matrix(X)) # as.matrix() in case p = 1 (auto. becomes vector)

  k <- (p + 2) * q + 1

  weight_matrix_init <- matrix(stats::runif(n_init * k, min = -unif, max = unif), ncol = k)

  weight_matrix <- matrix(rep(NA, n_init * k), ncol = k)
  inf_crit_vec <- rep(NA, n_init)
  converge <- rep(NA, n_init)

  for (iter in 1:n_init) {
    nn_model <- nnet::nnet(y ~ .,
      data = df, size = q, trace = F, linout = T,
      Wts = weight_matrix_init[iter, ], maxit = maxit, ...
    )
    weight_matrix[iter, ] <- nn_model$wts

    RSS <- nn_model$value
    sigma2 <- RSS / n

    log_likelihood <- (-n / 2) * log(2 * pi * sigma2) - RSS / (2 * sigma2)

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
#' @return Prediction for given inputs
#' @export
nn_pred <- function(X, W, q) {
  n <- nrow(X)
  p <- ncol(X)

  if (length(W) == ((p + 2) * q + 1)) {
    X <- cbind(rep(1, n), X)

    h_input <- X %*% t(matrix(W[1:((p + 1) * q)], nrow = q, byrow = TRUE))

    h_act <- cbind(rep(1, n), sigmoid(h_input))

    y_hat <- h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)

    return(y_hat)
  } else {
    return(print("Error: Incorrect number of weights for NN structure"))
  }
}