#' Sigmoid activation function
#'
#'
#' @param x Input
#' @return Sigmoid function
#' @export
sigmoid = function(x) 1/(1+exp(-x))


#' Selects which hidden unit to remove based off of BIC
#'
#' This function calculates the influence after of each hidden unit in the model
#' and determines which unit has the least influence on the output
#'
#' @param W Weight vector
#' @param q Number of hidden nodes in model
#' @param dataX Matrix of inputs
#' @param Y Output vector
#' @param inf_crit Information criterion
#' @return The node to remove
#' @export
remove_unit = function(W, dataX, Y, q, inf_crit = 'BIC'){
  n = nrow(dataX)
  p = ncol(dataX)
  k = ((p + 2)*q + 1)
  dataX = cbind(rep(1, n), dataX)
  influence = rep(NA, q)
  if(length(W) == k){
    for(j in 1:q){
      W_removed = W
      W_removed[((p + 1)*q + j + 1)] = 0

      h_input = dataX %*% t(matrix(W_removed[1:((p + 1)*q)], nrow = q, byrow = T))
      h_act = cbind(rep(1, n), sigmoid(h_input))

      y_hat = h_act %*% matrix(W_removed[c((length(W_removed) - q):length(W_removed))], ncol = 1)

      SSE = sum((y_hat - Y)^2)
      sigma2 = SSE/n
      log_likelihood = (-n/2)*log(2*pi*sigma2) - SSE/(2*sigma2)
      influence[j] = ifelse(inf_crit == 'BIC', -2*log_likelihood + log(n)*(k+1-p-2),
                            ifelse(inf_crit == 'AIC', -2*log_likelihood + 2*(k+1-p-2),
                                   ifelse(inf_crit == 'AICc', -2*log_likelihood + 2*(k+1)*(n/(n - (k+1-p-2) - 1)), NA)))
    }


    return(which.min(influence))

  }else{
    return(print('Error: Incorrect number of weights for NN structure'))
  }
}


#'Function to rearrange weights due to weight space symmetries
#' @export
weight_recover = function(W_true, W_pred, p, q){
  true  = cbind(matrix(W_true[1:(q*(p + 1))], byrow=T, ncol=(p + 1)), W_true[((p + 1)*q + 2):((p + 1)*q + q + 1)])
  pred  = cbind(matrix(W_pred[1:(q*(p + 1))], byrow=T, ncol=(p + 1)), W_pred[((p + 1)*q + 2):((p + 1)*q + q + 1)])

  bias_t = W_true[q*(p + 1) + 1]
  bias_p = W_pred[q*(p + 1) + 1]

  mat = matrix(NA, nrow = q, ncol = q*2)

  for(i in 1:q){
    mat[i,]=c(rowSums((matrix(rep(pred[i,], q),nrow = q, byrow = T) - true)^2),
              rowSums((matrix(rep(-pred[i,], q), nrow = q, byrow = T) - true)^2))
  }

  ind = apply(mat, 1 , which.min)

  for( i in 1:q){
    if(ind[i] %in% c((q + 1):(2*q))){
      ind[i] = ind[i] - q
      bias_p = bias_p + pred[i, ncol(pred)]
      pred[i,] = -pred[i,]
    }
  }
  predW = pred[order(ind),]
  predW_vec = c(t(predW[,1:(p + 1)]), bias_p, t(predW[, (p + 2)]))
  return(predW_vec)
}


#' @export
my_runif = function(n = 1, unif = 1, lower = 0.5){
  x = runif(n, lower, unif)
  y = sample(c(-1, 1), n, replace=T)
  return(x*y)
}

#Function is required for nn_fit optimization
#' Calculates normal log-likelihood of neural network
#'
#' @param W Weight vector
#' @param X Input data
#' @param Y Output data
#' @param q Number of hidden units
#' @return log-likelihood value
#' @export
log_likelihood = function(W, X, Y, q){
  n = nrow(X)
  p = ncol(X)

  if(length(W) == ((p+2)*q + 1)){
    X = cbind(rep(1, n), X)

    h_input = X %*% t(matrix(W[1:((p + 1)*q)], nrow = q, byrow = T))

    h_act = cbind(rep(1, n), sigmoid(h_input))

    y_hat = h_act %*% matrix(W[c((length(W) - q):length(W))], ncol = 1)
    SSE = sum((y_hat - Y)^2)
    sigma2 = SSE/n

    return(-((-n/2)*log(2*pi*sigma2) - SSE/(2*sigma2)))
  }else{
    return(print('Error: Incorrect number of weights for NN structure'))
  }
}
