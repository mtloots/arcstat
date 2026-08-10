# arcstat

Inference from the **arc length** of statistical functions. Most functionals summarise a curve by
its height somewhere; arc length summarises it by how far you travel along it, which makes it a
measure of curve complexity rather than of location or spread, and sensitive to structure the usual
functionals average away.

The package collects the arc-length programme:

* a **goodness-of-fit test** on the arc length of the probability plot, with an analytic saddlepoint
  null — powerful against local density structure (multimodality, clustering, heaping) where the
  empirical-distribution tests are weak, and weak against smooth location and scale departures where
  they are strong;
* **distributions built from arc length** — the arc-length generator and the quantile arc-length
  family, estimated by L-moments;
* the **characteristic-function** arc length, scale free, with closed forms for several families
  (the normal carries exactly 2, the exponential exactly pi);
* the **equivalence family**, where the shoulder equation collapses from transcendental to quadratic;
* the **kappa-4 induction curve** and the **circular** arc-length family;
* a **Bayesian** nonparametric arc-length test on the Dirichlet-process posterior.

One shared pure-C back-end serves two front ends:

	arcstat/      R package      (.C bindings; R CMD check: Status OK)
	arcstat-py/   Python package (ctypes; compiles the same C on first import)

The two are checked against each other by the `parity_*.sh` harnesses at the root of this repo.
