#!/bin/sh
# Parity harness for the closed-form equivalence-family module (arceq2.c).
# R and Python compile the SAME C source; every value below must agree byte for byte.
set -e
OUT=${TMPDIR:-/tmp}/parity_eq2
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
# refuse to run against a stale copy of the shared back-end (see the guard for why)
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py" || exit 1
Rscript -e '
suppressMessages(library(arcstat))
u  <- (1:9)/10
al <- -0.60; be <- -0.3778
vals <- c(
  bc_q(u, al, be), bc_cdf(bc_q(u, al, be), al, be), bc_pdf(bc_q(u, al, be), al, be),
  eq_ub_quad(al, be), eq_ub_quad(-0.75, -0.5),
  eq_E(al, be), eq_E(-0.55, -0.10), eq_E(-0.62, -0.60),
  eq_bstar(-0.60), eq_bstar(-0.55), eq_bstar(-0.52))
cat(sprintf("%.17g\n", vals), sep="")
' > "$OUT/r.txt"

python3 - <<PY > "$OUT/py.txt"
import sys
sys.path.insert(0, "$PKG/arcstat-py/src")
from arcstat import _core as c
u  = [i/10 for i in range(1,10)]
al, be = -0.60, -0.3778
vals = []
q = c.bc_q(u, al, be)
vals += q
vals += c.bc_cdf(q, al, be)
vals += c.bc_pdf(q, al, be)
vals += c.eq_ub_quad(al, be)
vals += c.eq_ub_quad(-0.75, -0.5)
vals += [c.eq_E(al, be), c.eq_E(-0.55, -0.10), c.eq_E(-0.62, -0.60)]
vals += [c.eq_bstar(-0.60), c.eq_bstar(-0.55), c.eq_bstar(-0.52)]
print("\n".join("%.17g" % v for v in vals))
PY

n=$(wc -l < "$OUT/r.txt" | tr -d " ")
if diff -q "$OUT/r.txt" "$OUT/py.txt" > /dev/null; then
  echo "PARITY OK: $n values byte-identical"
else
  echo "PARITY FAIL"; diff "$OUT/r.txt" "$OUT/py.txt" | head
  exit 1
fi
