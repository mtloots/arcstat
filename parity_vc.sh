#!/bin/sh
# Parity harness for the variance-components and locus-geometry module (arcvc.c).
#
# Both routines back stated results in the induction-curve paper: the intraclass correlation
# decides how much evidence the replicate runs actually carry, and the distance to the equivalence
# locus is what the price of the constraint is compared against. R and Python compile the SAME C
# source, so every value below must agree byte for byte at full double precision, with no tolerance.
#
# The grouping vector is deliberately UNBALANCED (3, 4, 5 and 6 members) because that is the case
# the paper is in and it is the case the balanced shortcut gets wrong: dividing the between-group
# mean square by the mean group size rather than by k0 = (N - sum n_i^2 / N)/(G - 1) biases the
# ratio, and a balanced test would not detect it.
set -e
OUT=${TMPDIR:-/tmp}/parity_vc
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py" || exit 1

Rscript -e '
suppressMessages(library(arcstat))
## deterministic, unbalanced, and spanning a wide range of separation
y <- c(1.10,1.22,0.97,
       5.05,5.31,4.88,5.12,
       9.02,9.19,8.77,9.31,8.95,
       2.40,2.55,2.61,2.38,2.49,2.44)
g <- rep(c("a","b","c","d"), times=c(3,4,5,6))
## a near-null case: groups that do not separate at all, where the between-group variance estimate
## goes negative and is truncated. The two fronts must agree on the truncation, not merely on the
## easy case.
y2 <- c(1.0,3.0,2.0, 2.1,0.9,3.1, 1.9,3.2,0.8)
g2 <- rep(c("x","y","z"), each=3)
## built by the SAME arithmetic on both sides: seq() and a Python comprehension do not
## agree in the last bit, and a differing INPUT masquerades as a differing back end
lh <- 0.05 + 0.02*(0:17)
lk <- -0.037 + 0.05*(0:17)
hh <- c(lh[5], 0.30, 1.175, 0.00, 0.42)
kk <- c(lk[5], -0.90, -0.98, 0.27, 0.57)
vals <- c(k4_icc(y, g), k4_icc(y2, g2), k4_locus_dist(hh, kk, lh, lk))
cat(sprintf("%.17g\n", vals), sep="")
' > "$OUT/r.txt"

python3 - <<PY > "$OUT/py.txt"
import sys
sys.path.insert(0, "$PKG/arcstat-py/src")
from arcstat import _core as c
y = [1.10,1.22,0.97,
     5.05,5.31,4.88,5.12,
     9.02,9.19,8.77,9.31,8.95,
     2.40,2.55,2.61,2.38,2.49,2.44]
g = ["a"]*3 + ["b"]*4 + ["c"]*5 + ["d"]*6
y2 = [1.0,3.0,2.0, 2.1,0.9,3.1, 1.9,3.2,0.8]
g2 = ["x"]*3 + ["y"]*3 + ["z"]*3
lh = [0.05 + 0.02*i for i in range(18)]
lk = [-0.037 + 0.05*i for i in range(18)]
hh = [lh[4], 0.30, 1.175, 0.00, 0.42]
kk = [lk[4], -0.90, -0.98, 0.27, 0.57]
vals = c.icc_oneway(y, g) + c.icc_oneway(y2, g2) + c.locus_dist(hh, kk, lh, lk)
for v in vals:
    print("%.17g" % v)
PY

if cmp -s "$OUT/r.txt" "$OUT/py.txt"; then
  echo "PARITY OK: $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical (arcvc)"
else
  echo "PARITY FAILED (arcvc)"
  diff "$OUT/r.txt" "$OUT/py.txt" | head -20
  exit 1
fi
