
# arcstat

[![CRAN
status](https://www.r-pkg.org/badges/version/arcstat)](https://CRAN.R-project.org/package=arcstat)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/arcstat)](https://CRAN.R-project.org/package=arcstat)
[![R-CMD-check](https://github.com/mtloots/arcstat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mtloots/arcstat/actions/workflows/R-CMD-check.yaml)

Inference from the **arc length** of statistical functions: a
goodness-of-fit test on the arc length of the probability plot,
distributions built from the arc length of their defining curve, the
characteristic-function version, the equivalence families, and a
Bayesian test; and the induction curves of the oxidation-stability
standards, fitted, bootstrapped and compared in batches of many runs.

Computation runs on a shared pure-C back-end, parallel by OpenMP
(`arc_threads()`), that is also bound from Python
(`pip install arcstat`), and the two front ends are checked against each
other value by value.

## Installation

Install the released version from CRAN:

``` r
install.packages("arcstat")
```

The development version, which may be ahead of CRAN:

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

Induction curves: a kappa curve fitted to each of the eleven runs of
`DNase`, and the two standard readings of each fit.

``` r
d <- split(DNase, DNase$Run)
f <- k4_fit(lapply(d, function(z) log(z$conc)), lapply(d, function(z) z$density),
            bounds = c(-0.98, 0.95, -5, 4), drift = FALSE)
round(head(f[, c("k", "h", "a", "b")], 4), 3)
#>         k      h      a      b
#> 10 -0.028 -0.672 -0.788 -0.250
#> 11  0.931  0.308 -0.263  0.851
#> 9   0.102 -1.231 -0.721  0.212
#> 1   0.063 -0.914 -0.623  0.157
```

## Citation

`citation("arcstat")`.
