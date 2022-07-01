#' Hidden Node Selection
#'
#' Performs either bulk or step-wise hidden node selection
#'
#' @param X Matrix of covariates
#' @param y Vector of response
#' @param Q Possible number of hidden nodes
#' @param n_init Number of random initialisations (tracks)
#' @param type bulk or step
#' @param inf_crit Information criterion
#' @param unif Random initial values max value
#' @param ... additional argument for nnet
#' @return Optimal number of hidden nodes
#' @export
hidden_node_sel = function(X, y, Q, n_init, type = 'bulk', inf_crit = 'BIC', unif = 3, ...){

  if (type == 'bulk') {

    inf_crit_vec <- rep(NA, Q)

    names(inf_crit_vec) <- as.character(1:Q)

    weights_min <- vector('list', length = Q)

    for(q in 1:Q){

      nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, unif, ...)

      weights_min[[q]] <- nn$W_opt

      inf_crit_vec[q] <- nn$val

    }

    W_opt <- weights_min[[which.min(inf_crit_vec)]]

  } else if (type == 'step') {

    n_candidates <- 3

    inf_crit_vec <- rep(NA, n_candidates)

    names(inf_crit_vec) <- as.character((Q - 1):(Q + 1))

    n = nrow(X) # sample size
    p = ncol(X) # number of covariates

    weights_min = vector(mode = "list", length = n_candidates)
    names(weights_min) <- as.character((Q - 1):(Q + 1))


    for (q in (Q - 1):(Q + 1)) {

      if (q == 0){
        weights_min[[as.character(q)]] <- NULL

        inf_crit_vec[as.character(q)] <- NA
      } else {

        k = (p + 1)*q + (q + 1)


        nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, unif, ...)

        weights_min[[as.character(q)]] <- nn$W_opt

        inf_crit_vec[as.character(q)] <- nn$val

      }

    }

    W_opt <- weights_min[[names(which.min(inf_crit_vec))]]

  }

  return(list('min' = as.numeric(names(which.min(inf_crit_vec))),
              'val' = min(inf_crit_vec, na.rm = TRUE),
              'inf_crit_vec' = inf_crit_vec,
              'weights_min' = weights_min,
              'W_opt' = W_opt))
}
