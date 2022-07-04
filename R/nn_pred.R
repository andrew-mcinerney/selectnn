#' Network Output Function
#'
#'
#' @param X Data
#' @param W Weight vector
#' @param q Number of hidden units
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
