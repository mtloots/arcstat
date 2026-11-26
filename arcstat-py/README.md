# arcstat

Arc-length statistics, equivalence families and induction curves — one shared, pure-C back-end,
parallel by OpenMP, bound identically from **R** (package `arcstat`) and **Python** (this package).

The C sources are compiled on first import; there are no third-party dependencies.

## Tools

- **Goodness of fit** — a test based on the arc length of the probability plot, with an analytic
  saddlepoint null, sensitive to local density structure (multimodality, clustering, heaping) that the
  Kolmogorov–Smirnov, Cramér–von Mises and Anderson–Darling tests miss:
  `al_statistic`, `al_pvalue`, `al_moments`.
- **Distributions** — build a distribution from the arc length of its defining curve: the quantile
  arc-length family and its L-moment fitting: `arcq_qd`, `arcq_arclength`, `sample_lmoments`.
- **Bayesian test** — a Dirichlet-process / Bayesian-bootstrap arc-length goodness-of-fit test:
  `bb_post_disc`, `bb_ref_disc`, `bb_evidence`.

- **Induction curves** — kappa, inverse Gaussian and beta-companion curves with drift, many runs in one
  call, with bootstraps and a permutation test of a shared shape: `k4_fit`, `k4eq_fit`, `k4_shared_fit`,
  `k4_shared_perm`, `k4_boot`, `k4eq_boot`, `wald_fit`, `bc_host_fit`; `set_threads(n)`.

```python
import arcstat
arcstat.al_pvalue(1.567, 50)                 # saddlepoint p-value
arcstat.arcq_arclength([0.3, -0.2], 1.5)     # arc length of a quantile curve
arcstat.bb_evidence(u)                        # Bayesian evidence against H0
```

The same C sources back the R package `arcstat`, and the two return identical values.

Author: M. Theodor Loots. Licence: GPL-3.
