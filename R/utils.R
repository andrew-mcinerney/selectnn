#' Sigmoid activation function
#'
#'
#' @param x Input
#' @return Sigmoid function
#' @export
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}
