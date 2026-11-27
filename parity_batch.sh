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
u <- k4_fit(x, y, gate = TRUE, gev = TRUE); w <- k4_fit(x, y, gate = "best", keep = 4L); z <- k4_fit(x, y, hfix = -1)
wd <- wald_fit(x, y); w0 <- wald_fit(x, y, drift = FALSE); f0 <- k4_fit(x, y, drift = FALSE); q0 <- k4eq_fit(x, y, drift = FALSE)
wr <- c(wald_cdf(c(0, 1, 5, 20), 9, 0.4), wald_cdf(c(0.1, 9), 9, 2e6), t(wald_readings(c(9, 1e8, 3), c(150, 15, 3e7))))
lw <- k4_lmom_window(x[1], y[1], rbind(c(1, 0.2, 300)), ngrid = 500L, nmodel = 500L)
eg <- k4eq_fit(x, y, gate = "best"); yf <- 300 - 250*k4_cdf(x[[1]], 9, 2, 0.1, -0.5) + sin(4321*(1:120)); en <- k4eq_fit(x[1], list(yf), gate = "best"); ki <- k4_fit(x, y, gate = "best", infimum = TRUE)
qg <- k4eq_boot(x, y, eg$fit, B = 4L, gate = "best")
xs <- c(x, x[1], x[2]); ys <- c(y, list(y[[1]] + cos(77*(1:120))), list(y[[2]] + cos(55*(1:90))))
sf <- k4_shared_fit(xs, ys, sets = list(1:4, c(1,3), c(2,4), 2:3), bounds = c(-3, 0.95, -3, 4)); sp <- k4_shared_perm(xs, ys, group = c(1,2,1,2), rss_free = 1000, slack = c(1e-9, 0.02))
kd <- k4_fit(x, y, bounds = c(-5, 0.95, -60, 4), gate = "best", dense = TRUE); b2 <- bc_host_fit(x, y, keep = 2L); mv <- k4_mv_boot(x[[2]], y[[2]], B = 2L, bounds = box, drift = TRUE)
nb0 <- k4eq_boot(x, y, q0$fit, B = 4L); nb1 <- k4eq_boot(x, y, q0$fit, B = 3L, drift = TRUE); nf0 <- k4_boot(x, y, f0, box, B = 3L, gate = "best")
v <- function(a) as.vector(aperm(a, c(2, 1, 3)))
vals <- c(k4_locus_k(c(-60, -10, -1, -1e-9, 0, 0.3, 0.4999, 0.6)), t(f$fit), aperm(f$pieces, c(2, 3, 1)), t(g$fit), aperm(g$pieces, c(2, 3, 1)),
          t(b$free), t(b$con), v(m$boot), m$block, v(i$boot), i$block, v(q$boot), q$block, v(r$boot), r$block, t(u), t(w), t(z), t(b2$free), t(b2$con), mv$a_pt, t(mv$A), t(wd), t(w0), t(f0), t(q0$fit), wr, lw, t(eg$fit), aperm(eg$pieces, c(2, 3, 1)), t(eg$physical), t(en$fit), t(en$physical), t(ki), v(qg$boot), k4_mv_share(mv$A[, 1:4], 0.1), k4_mv_share(i$boot[, c("a","k","h"), 1]), t(sf$fit), t(do.call(rbind, sf$runs)), t(sp$stat), sp$p, t(kd), v(nb0$boot), nb0$block, v(nb1$boot), v(nf0$boot), nf0$block)
cat(ifelse(is.nan(vals), "NaN", sprintf("%.17g", vals)), sep="\n")
' > "$OUT/r.txt"
for T in "" 1 2 8; do
PYTHONPATH="$PKG/arcstat-py/src" ARCSTAT_THREADS=$T python3 - <<PY > "$OUT/py$T.txt"
import math
from arcstat import bc_host_fit, eq_bstar, k4_boot, k4_cdf, k4_fit, k4_lmom_window, k4_locus_k, k4_mv_boot, k4_mv_share, k4_shared_fit, k4_shared_perm, k4eq_boot, k4eq_fit, wald_cdf, wald_fit, wald_readings
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
u = k4_fit(x, y, gate=True, gev=True); w = k4_fit(x, y, gate="best", keep=4); z = k4_fit(x, y, hfix=-1)
wd = wald_fit(x, y); w0 = wald_fit(x, y, drift=False); f0 = k4_fit(x, y, drift=False); q0 = k4eq_fit(x, y, drift=False)[0]
wr = wald_cdf([0, 1, 5, 20], 9, 0.4) + wald_cdf([0.1, 9], 9, 2e6) + wald_readings([9, 1e8, 3], [150, 15, 3e7])
lw = k4_lmom_window(x[:1], y[:1], [(1, 0.2, 300)], ngrid=500, nmodel=500)
eg = k4eq_fit(x, y, gate="best"); F3 = k4_cdf(x[0], 9, 2, 0.1, -0.5); yf = [300 - 250*F3[j] + math.sin(4321*(j+1)) for j in range(120)]
en = k4eq_fit(x[:1], [yf], gate="best"); ki = k4_fit(x, y, gate="best", infimum=True)
qg = k4eq_boot(x, y, eg[0], B=4, gate="best")[0]
xs = x + [x[0], x[1]]; ys = y + [[y[0][j] + math.cos(77*(j+1)) for j in range(120)], [y[1][j] + math.cos(55*(j+1)) for j in range(90)]]
sf = k4_shared_fit(xs, ys, sets=[[0,1,2,3],[0,2],[1,3],[1,2]], bounds=(-3, 0.95, -3, 4)); sp = k4_shared_perm(xs, ys, group=[1,2,1,2], rss_free=[1000], slack=(1e-9, 0.02))
kd = k4_fit(x, y, bounds=(-5, 0.95, -60, 4), gate="best", dense=True); b2 = bc_host_fit(x, y, keep=2); mv = k4_mv_boot(x[1], y[1], B=2, bounds=box, drift=True)
nb0 = k4eq_boot(x, y, q0, B=4); nb1 = k4eq_boot(x, y, q0, B=3, drift=True); nf0 = k4_boot(x, y, f0, box, B=3, gate="best")
def flat(o):
    return [v for e in o for v in flat(e)] if isinstance(o, (list, tuple)) else [o]
vals = flat([k4_locus_k([-60, -10, -1, -1e-9, 0, 0.3, 0.4999, 0.6]), f, fp, g, gp, fr, co, m, mb, i, ib, q, qb, r, rb, u, w, z, b2, mv, wd, w0, f0, q0, wr, lw, eg, en[0], en[2], ki, qg, k4_mv_share([r[:4] for r in mv[1]], 0.1), k4_mv_share([[r[0], r[2], r[3]] for r in i[0]]), sf[0], [r for s in sf[1] for r in s], sp[0], sp[1], kd, nb0, nb1[0], nf0])
print("\n".join("NaN" if math.isnan(v) else "%.17g" % v for v in vals))
PY
cmp -s "$OUT/r.txt" "$OUT/py$T.txt" || { echo "VERDICT: MISMATCH (threads '$T')"; diff "$OUT/r.txt" "$OUT/py$T.txt" | head; exit 1; }
done
echo "VERDICT: PARITY -- $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical, on the default, 1, 2 and 8 threads"
