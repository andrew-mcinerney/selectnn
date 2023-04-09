#' Sigmoid activation function
#'
#'
#' @param x Input
#' @return Sigmoid function
#' @export
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' Number of weights in NN
#'
#'
#' @param p number of inputs
#' @param structure of hidden layers
#' @return number of parameters
#' @export
n_params <- function (p, q) {
  return(sum(c(p + 1, q + 1) * c(q, 1)))
}
