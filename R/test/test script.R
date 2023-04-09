n <- 1000
p <- 4
q <- c(3, 1)
n_init <- 5
maxit <- 100000

k <- sum(c(p + 1, q + 1) * c(q, 1))

W <- nnic::my_runif(k, 3, 2)

X <- matrix(rnorm(n * p), nrow = n)

y <- nn_pred(X, W, q)

L <- 3

Q <- rep(7, L)

h <- hidden_layer_sel(X, y, Q, n_init = n_init, nn_fn = "neuralnet", maxit = maxit)


h$min

sum((nn_pred(X, h$weights_min[[1]], Q[1]) - y)^2)
sum((nn_pred(X, h$weights_min[[2]], Q[1:2]) - y)^2)
sum((nn_pred(X, h$weights_min[[3]], Q[1:3]) - y)^2)
h$inf_crit_vec



# test --------------------------------------------------------------------

nn_layers <- vector(mode = "list", length = 1)

nn_hidden <- vector(mode = "list", length = 1)

nn_layers[[1]] <- hidden_layer_sel(
  X = X,
  y = y,
  Q = Q,
  n_init = n_init,
  nn_fn = "neuralnet",
  maxit = maxit
  )

q_temp <- nn_layers[[1]]$q

for (l in length(nn_layers[[1]]$q):1) {
  nn_hidden[[1]] <- hidden_node_sel1(
    X = X,
    y = y,
    Q = q_temp,
    n_init = n_init,
    nn_fn = "neuralnet",
    maxit = maxit,
    l = l
  )

  q_temp[l] <- nn_hidden[[1]]$min
}

q_temp

nn_layers[[1]]
nn_hidden[[1]]
