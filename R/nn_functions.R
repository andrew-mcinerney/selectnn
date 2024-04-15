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
#' @param nn_fn the neural network optimiser to use
#' @param ... additional argument for nnet
#' @return The best model from the different tracks
#' @export
nn_fit_tracks <- function(X, y, q, n_init, inf_crit = "BIC",
                              task = "regression", unif = 3, maxit = 1000,
                              nn_fn = "nnet", ...) {
  # Function with fits n_init tracks of model and finds best

  if ((length(q) > 1) & (nn_fn == "nnet")) {
    nn_fn <- "neuralnet"
    warning("nnet does not allow more than one hidden layer, neuralnet used instead")
  }

  df <- data.frame(X, y)
  n <- nrow(X)
  p <- ncol(as.matrix(X)) # as.matrix() in case p = 1 (auto. becomes vector)

  k <- sum(c(p + 1, q + 1) * c(q, 1))

  weight_matrix_init <- matrix(stats::runif(n_init * k, min = -unif, max = unif), ncol = k)

  weight_matrix <- matrix(rep(NA, n_init * k), ncol = k)
  inf_crit_vec <- rep(NA, n_init)

  if (task == "regression") {
    linear.output <- TRUE
  } else if (task == "classification") {
    linear.output <- FALSE
  } else {
    stop(sprintf(
      "Error: %s not recognised as task. Please choose regression or classification",
      task
    ))
  }

  for (iter in 1:n_init) {

    if (nn_fn == "nnet") {

      nn_model <- nnet::nnet(X, y, size = q, trace = FALSE,
                             linout = linear.output,
                             Wts = weight_matrix_init[iter, ],
                             maxit = maxit, ...
      )

      weight_matrix[iter, ] <- nn_model$wts

    } else if (nn_fn == "neuralnet") {

      nn_model <- neuralnet::neuralnet(y ~ .,
                                       data = df, hidden = q,
                                       linear.output = linear.output,
                                       startweights = weight_matrix_init[iter, ],
                                       stepmax = maxit, ...
      )

      weight_matrix[iter, ] <- unlist(sapply(nn_model$weights[[1]], as.vector))

    } else {
      stop("nn_fn not recognised")
    }

    if (task == "regression") {
      RSS <- sum((y - nn_pred(X, weight_matrix[iter, ], q))^2)
      sigma2 <- RSS / n

      log_likelihood <- (-n / 2) * log(2 * pi * sigma2) - RSS / (2 * sigma2)
    } else if (task == "classification") {
      # need to write likelihood value for binary case
      # log_likelihood <- -nn_model$value
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
  }
  W_opt <- weight_matrix[which.min(inf_crit_vec), ]

  return(list(
    "W_opt" = W_opt,
    "value" = min(inf_crit_vec),
    "inf_crit_vec" = inf_crit_vec,
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
nn_pred <- function(X, W, q, output = "identity") {
  n <- nrow(X)
  p <- ncol(X)

  k <- sum(c(p + 1, q + 1) * c(q, 1))

  layer_nodes <- c(0, cumsum(c(p + 1, q + 1) * c(q, 1)))

  if (length(W) == k) {

    X <- cbind(rep(1, n), X)

    temp <- X

    for(i in 1:length(q)) {

      h_input <- temp %*% t(matrix(W[(layer_nodes[i] + 1):layer_nodes[i + 1]],
                                   nrow = q[i], byrow = TRUE))

      h_act <- cbind(rep(1, n), sigmoid(h_input))

      temp <- h_act
    }


    if (output == "identity") {
      y_hat <- h_act %*%
        matrix(W[(layer_nodes[length(layer_nodes) - 1] + 1):
                   layer_nodes[length(layer_nodes)]],
               ncol = 1)
    } else if (output == "sigmoid") {
      y_hat <- sigmoid(
        h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)
      )
    } else {
      stop(
        sprintf(
          "Error: %s not recognised as available output function.",
          output
        )
      )
    }

    return(y_hat)
  } else {
    stop(sprintf(
      "Error: Incorrect number of weights for NN structure. W should have
      %s weights (%s weights supplied).", k, length(W)
    ))
  }
}
