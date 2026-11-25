#!/bin/sh
# Parity harness for the batched fits and bootstraps (arcbatch.c) and the thread setting
# (arcthreads.c). Every value must agree byte for byte between the fronts, and the whole Python
# vector again on 1, 2 and 8 threads.
set -e
OUT=${TMPDIR:-/tmp}/parity_batch
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py/src/arcstat" || exit 1
Rscript -e '
suppressMessages(library(arcstat))
xx <- function(n) 20*(0:(n-1))/(n-1)
x <- list(xx(120), xx(90))
y <- list(1 + 0.2*x[[1]] + 300*k4_cdf(x[[1]], 9, 2, -0.3, 0.2) + 2*sin(4321*(1:120)),
          2 - 0.1*x[[2]] + 250*k4_cdf(x[[2]], 8, 1.5, 0.1, -0.4) + 2*sin(977*(1:90)))
h <- -1.2 + 0.15*(0:11); tab <- list(cbind(h[1:8], k4_locus_k(h[1:8])), cbind(h[9:11], k4_locus_k(h[9:11])))
al <- -0.62 + 0.11*(0:7)/7; curve <- cbind(al, vapply(al, eq_bstar, 0))
f <- k4eq_fit(x, y); g <- k4eq_fit(x, y, tab)
b <- bc_host_fit(x, y, curve, keep = 3L)
box <- c(-0.98, 0.95, -60, 4)
fit <- rbind(c(1, 0.2, 300, 9, 2, -0.3, 0.2), c(2, -0.1, 250, 8, 1.5, 0.1, 0))
m <- k4_boot(x, y, fit, box, gev = c(FALSE, TRUE), B = 5L); i <- k4_boot(x, y, fit, box, B = 3L, block = 1L, seed = c(5L, 77L))
q <- k4eq_boot(x, y, f$fit, B = 5L); r <- k4eq_boot(x, y, g$fit, tab, B = 2L)
v <- function(a) as.vector(aperm(a, c(2, 1, 3)))
vals <- c(k4_locus_k(c(-60, -10, -1, -1e-9, 0, 0.3, 0.4999, 0.6)), t(f$fit), aperm(f$pieces, c(2, 3, 1)), t(g$fit), aperm(g$pieces, c(2, 3, 1)),
          t(b$free), t(b$con), v(m$boot), m$block, v(i$boot), i$block, v(q$boot), q$block, v(r$boot), r$block)
cat(ifelse(is.nan(vals), "NaN", sprintf("%.17g", vals)), sep="\n")
' > "$OUT/r.txt"
for T in "" 1 2 8; do
PYTHONPATH="$PKG/arcstat-py/src" ARCSTAT_THREADS=$T python3 - <<PY > "$OUT/py$T.txt"
import math
from arcstat import bc_host_fit, eq_bstar, k4_boot, k4_cdf, k4_locus_k, k4eq_boot, k4eq_fit
xx = lambda n: [20*i/(n-1) for i in range(n)]
x = [xx(120), xx(90)]
F1 = k4_cdf(x[0], 9, 2, -0.3, 0.2); F2 = k4_cdf(x[1], 8, 1.5, 0.1, -0.4)
y = [[1 + 0.2*x[0][i] + 300*F1[i] + 2*math.sin(4321*(i+1)) for i in range(120)],
     [2 - 0.1*x[1][i] + 250*F2[i] + 2*math.sin(977*(i+1)) for i in range(90)]]
h = [-1.2 + 0.15*i for i in range(12)]; k = k4_locus_k(h[:11])
tab = [list(zip(h[:8], k[:8])), list(zip(h[8:11], k[8:11]))]
al = [-0.62 + 0.11*i/7 for i in range(8)]; curve = [(a, eq_bstar(a)) for a in al]
f, fp = k4eq_fit(x, y); g, gp = k4eq_fit(x, y, tab)
fr, co = bc_host_fit(x, y, curve, keep=3)
box = (-0.98, 0.95, -60, 4)
fit = [(1, 0.2, 300, 9, 2, -0.3, 0.2), (2, -0.1, 250, 8, 1.5, 0.1, 0)]
m, mb = k4_boot(x, y, fit, box, gev=[False, True], B=5); i, ib = k4_boot(x, y, fit, box, B=3, block=1, seed=[5, 77])
q, qb = k4eq_boot(x, y, f, B=5); r, rb = k4eq_boot(x, y, g, tab, B=2)
def flat(o):
    return [v for e in o for v in flat(e)] if isinstance(o, (list, tuple)) else [o]
vals = flat([k4_locus_k([-60, -10, -1, -1e-9, 0, 0.3, 0.4999, 0.6]), f, fp, g, gp, fr, co, m, mb, i, ib, q, qb, r, rb])
print("\n".join("NaN" if math.isnan(v) else "%.17g" % v for v in vals))
PY
cmp -s "$OUT/r.txt" "$OUT/py$T.txt" || { echo "VERDICT: MISMATCH (threads '$T')"; diff "$OUT/r.txt" "$OUT/py$T.txt" | head; exit 1; }
done
echo "VERDICT: PARITY -- $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical, on the default, 1, 2 and 8 threads"
