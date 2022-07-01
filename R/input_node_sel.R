#' Input Node Selection
#'
#' Performs either bulk or stepwise input node selection
#'
#' @param X Matrix of covarites
#' @param y Vector of response
#' @param q Number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param type bulk or step
#' @param inf_crit Information criterion
#' @param unif Random initial values max value
#' @return Inputs dropped from model
#' @export
input_node_sel <- function(X, y, q, n_init, type = 'bulk', inf_crit = 'BIC', unif = 3, ...){

  full_model <- nn_fit_tracks(X, y, q, n_init, inf_crit, unif, ...)

  full_inf_crit <- full_model$val

  p <- ncol(X)

  X_full <- X

  colnames(X_full) = 1:p

  colnames(X) = 1:p

  dropped <- c()

  W_opt <- full_model$W_opt

  min_inf_crit <- full_inf_crit

  continue_drop <- TRUE

  while(continue_drop == TRUE){

    nn_in <- input_importance(X = X, Y = y, q = q, n_iter = n_init,
                              inf_crit = inf_crit, unif = unif,
                              ...)

    if(nn_in$val < min_inf_crit){

      X = X[, -nn_in$min, drop = FALSE]
      p = ncol(as.matrix(X))


      W_opt = nn_in$W_opt

      dropped <- colnames(X_full)[!colnames(X_full) %in% colnames(X)]
      min_inf_crit <- nn_in$val

      if (ncol(X) == 1){
        continue_drop <- FALSE
      }

      if (type == 'step'){
        continue_drop <- FALSE
      }

    }else{
      continue_drop <- FALSE
    }
  }
  return(list('X' = X, 'p' = p, 'W_opt' = W_opt, 'dropped' = dropped,
              'full_inf_crit' = full_inf_crit, 'inf_crit' = min_inf_crit))
}
