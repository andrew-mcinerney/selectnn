hidden_layer_sel <- function(X, y, Q, n_init, L = 1, type = "bulk", inf_crit = "BIC",
                            task = "regression", unif = 3, maxit = 1000, ...) {
  if ((length(Q) > 1) & (L == 1)) {
    L <- length(Q)
  } else if ((length(Q) == 1) & (L > 1)) {
    Q <- rep(Q, times = L)
  } else if (length(Q) != L){
    stop("Error: Q and L do not agree")
  }


  if (type == "bulk") {
    inf_crit_vec <- rep(NA, L)

    names(inf_crit_vec) <- as.character(1:L)

    weights_min <- vector("list", length = L)

    for (l in 1:L) {
      q <- Q[1:l]
      nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, task, unif,
                          maxit = maxit, ...
      )

      weights_min[[l]] <- nn$W_opt

      inf_crit_vec[l] <- nn$value
    }

    W_opt <- weights_min[[which.min(inf_crit_vec)]]
  } else if (type == "step") {
    n_candidates <- 3

    inf_crit_vec <- rep(NA, n_candidates)

    names(inf_crit_vec) <- as.character((Q - 1):(Q + 1))

    n <- nrow(X) # sample size
    p <- ncol(X) # number of covariates

    weights_min <- vector(mode = "list", length = n_candidates)
    names(weights_min) <- as.character((Q - 1):(Q + 1))


    for (q in (Q - 1):(Q + 1)) {
      if (q == 0) {
        weights_min[[as.character(q)]] <- NULL

        inf_crit_vec[as.character(q)] <- NA
      } else {
        k <- (p + 1) * q + (q + 1)


        nn <- nn_fit_tracks(X, y, q, n_init, inf_crit, task, unif,
                            maxit = maxit, ...
        )

        weights_min[[as.character(q)]] <- nn$W_opt

        inf_crit_vec[as.character(q)] <- nn$value
      }
    }

    W_opt <- weights_min[[names(which.min(inf_crit_vec))]]
  }

  return(list(
    "min" = as.numeric(names(which.min(inf_crit_vec))),
    "q" = Q[1:as.numeric(names(which.min(inf_crit_vec)))],
    "value" = min(inf_crit_vec, na.rm = TRUE),
    "inf_crit_vec" = inf_crit_vec,
    "weights_min" = weights_min,
    "W_opt" = W_opt
  ))
}
