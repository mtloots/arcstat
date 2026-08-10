
# arcstat

[![R-CMD-check](https://github.com/mtloots/arcstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/arcstat/actions/workflows/R-CMD-check.yaml)

Inference from the **arc length** of statistical functions: a
goodness-of-fit test on the arc length of the probability plot,
distributions built from the arc length of their defining curve, the
characteristic-function version, the equivalence family, and a Bayesian
test.

Computation runs on a shared pure-C back-end that is also bound from
Python, and the two front ends are checked against each other value by
value.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mtloots/arcstat", subdir = "arcstat")
```

## Usage

``` r
library(arcstat)
set.seed(1)
al_test(runif(200))$p.value
#> [1] 0.3974529

## two closed forms that anchor the construction
c(normal = cf_arclength_family("normal"),
  exponential = cf_arclength_family("exponential", lambda = 1))
#>      normal exponential 
#>    2.000000    3.141593
```
