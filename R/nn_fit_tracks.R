#' Fits various tracks (different random starting values) and choses best model
#'
#' Fits n_init tracks with different initial values and decides on best model
#' based on information criteria
#'
#' @param X Matrix of covarites
#' @param y Vector of response
#' @param q Number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param inf_crit Information criterion
#' @param unif Random initial values max value
#' @return The best model from the different tracks
#' @export
nn_fit_tracks <- function (X, y, q, n_init, inf_crit = 'BIC', unif = 3, ...){
  # Function with fits n_iter tracks of model and finds best

  df <- data.frame(X, y)
  n <- nrow(X)
  p <- ncol(as.matrix(X)) # as.matrix() incase p = 1 (auto. becomes vector)

  k <- (p + 2)*q + 1

  weight_matrix_init <- matrix(stats::runif(n_init*k, min = -unif, max = unif), ncol = k)

  weight_matrix <- matrix(rep(NA, n_init*k), ncol = k)
  inf_crit_vec <- rep(NA, n_init)
  converge <- rep(NA, n_init)

  for(iter in 1:n_init){
    nn_model <- nnet::nnet(y~., data = df, size = q, trace = F, linout = T,
                           Wts = weight_matrix_init[iter, ], ...)
    weight_matrix[iter,] = nn_model$wts

    RSS = nn_model$value
    sigma2 = RSS/n

    log_likelihood = (-n/2)*log(2*pi*sigma2) - RSS/(2*sigma2)

    inf_crit_vec[iter] = ifelse(inf_crit == 'AIC',
                                (2*(k+1) - 2*log_likelihood),
                                ifelse(inf_crit == 'BIC',
                                       (log(n)*(k+1) - 2*log_likelihood),
                                       ifelse(inf_crit == 'AICc',
                                              (2*(k+1)*(n/(n-(k+1)-1)) - 2*log_likelihood),
                                              NA)))
    converge[iter] = nn_model$convergence
  }
  W_opt <- weight_matrix[which.min(inf_crit_vec), ]

  return(list(
    'W_opt' = W_opt,
    'val' = min(inf_crit_vec),
    'inf_crit_vec' = inf_crit_vec,
    'converge' = converge,
    'weight_matrix' = weight_matrix))
}
