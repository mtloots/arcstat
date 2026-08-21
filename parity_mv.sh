#!/bin/sh
# Parity harness for the multiverse bootstrap module (arcmv.c).
# The two fronts call the same C entry with the same synthetic curve; the eight point estimates
# and every replicate reading must agree byte for byte at full double precision. The resampling
# never touches a host RNG -- per-replicate splitmix64 streams live in the C -- so thread count,
# front language and platform scheduling are all outside the value path by construction.
set -e
OUT=${TMPDIR:-/tmp}/parity_mv
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py/src/arcstat" || exit 1
Rscript -e '
suppressMessages(library(arcstat))
th <- c(1, 400, 8, 2.5, 0.2, 0.4)
n  <- 200
x  <- 20*(0:(n-1))/(n-1)
y  <- th[1] + th[2]*k4_cdf(x, th[3], th[4], th[5], th[6]) + 2*sin(9871*(1:n))
mv <- k4_mv_boot(x, y, cut = 0.4, B = 6, seed = 4207)
vals <- c(unname(mv$a_pt), as.vector(t(mv$A)))
cat(sprintf("%.17g", vals), sep="\n")
' > "$OUT/r.txt"
PYTHONPATH="$PKG/arcstat-py/src" python3 - <<PY > "$OUT/py.txt"
import math
from arcstat import k4_cdf, k4_mv_boot
th = [1, 400, 8, 2.5, 0.2, 0.4]
n = 200
x = [20*i/(n-1) for i in range(n)]
F = k4_cdf(x, th[2], th[3], th[4], th[5])
y = [th[0] + th[1]*F[i] + 2*math.sin(9871*(i+1)) for i in range(n)]
a_pt, A = k4_mv_boot(x, y, cut=0.4, B=6, seed=4207)
vals = list(a_pt) + [v for row in A for v in row]
print("\n".join("%.17g" % v for v in vals))
PY
if cmp -s "$OUT/r.txt" "$OUT/py.txt"; then
  echo "VERDICT: PARITY -- $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical"
else
  echo "VERDICT: MISMATCH"; diff "$OUT/r.txt" "$OUT/py.txt" | head; exit 1
fi
