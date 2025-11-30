
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rerand

<!-- badges: start -->

<!-- badges: end -->

The **rerand** package provides minimal tools for implementing
rerandomization methods in randomized experiments. It includes functions
for generating treatment assignments under complete randomization (CRE)
and rerandomization (ReM), computing covariate-balance measures such as
the Mahalanobis distance, and evaluating whether an assignment satisfies
a prespecified balance criterion.

## Installation

You can install the development version of rerand like so:

``` r
# install.packages("remotes")
remotes::install_github("YOUR_USERNAME/rerand")
```

## Performance

The main rerandomization routines in **rerand** are implemented in C++
via Rcpp/RcppArmadillo, which can offer substantial speedups over
straightforward R implementations when the sample size or the number of
rerandomization attempts is large.

``` r
library(rerand)
library(microbenchmark)

set.seed(123)

n <- 5000
K <- 10
X <- matrix(rnorm(n * K), ncol = K)
Y <- rnorm(n)

# Example: compare a pure R implementation with the C++ implementation
bench <- microbenchmark(
  R   = rerand::ReM(X = X, Y = Y, n_1 = n / 2, p_accept = 0.1, engine = "R"),
  Cpp = rerand::ReM(X = X, Y = Y, n_1 = n / 2, p_accept = 0.1, engine = "cpp")
  times = 50
)

bench
```
