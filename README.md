
<!-- README.md is generated from README.Rmd. Please edit that file -->

# selectnn

<!-- badges: start -->

[![R-CMD-check](https://github.com/andrew-mcinerney/selectnn/workflows/R-CMD-check/badge.svg)](https://github.com/andrew-mcinerney/selectnn/actions)
<!-- badges: end -->

The goal of selectnn is to determine the optimal architecture when
implementing neural networks as a statistical model. The function used
to fit neural networks is `nnet` (Venables and Ripley, 2002).

## Installation

You can install the development version of selectnn from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("andrew-mcinerney/selectnn")
```

## nn_model_sel()

This function is used to perform neural network model selection.

``` r
library(selectnn)
nn_model_sel(X, y, Q, n_init)
nn_model_sel(y ~ ., data = df, Q, n_init)
```
