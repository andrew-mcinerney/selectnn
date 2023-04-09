
<!-- README.md is generated from README.Rmd. Please edit that file -->

# selectnn <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/andrew-mcinerney/selectnn/workflows/R-CMD-check/badge.svg)](https://github.com/andrew-mcinerney/selectnn/actions)
<!-- badges: end -->

This package implements the model selection approach detailed in
McInerney and Burke (2022): “A Statistically-Based Approach to
Feedforward Neural Network Model Selection”. The preprint of this paper
is available on [arXiv](https://arxiv.org/abs/2207.04248). More
specifically, the algorithm alternates between selecting the hidden
layer complexity and the inputs with the objective of minimizing the
BIC. The neural network function used is `nnet`, which is available from
the R package of the same name (Ripley and Venables, 2022).

## Installation

You can install the development version of selectnn from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("andrew-mcinerney/selectnn")
```

## Model Selection

The primary function in this package is `selectnn()`.

``` r
library(selectnn)
selectnn(X, y, Q, n_init)
selectnn(y ~ ., data = df, Q, n_init)
```
