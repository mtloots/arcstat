#!/bin/sh
# Parity harness for the equivalence-fitting module (arceqfit.c): the beta-companion fit, the
# estimator study, and the two grid sweeps, byte-identical between the R and Python fronts.
# Randomness lives in per-replicate splitmix64 streams inside the C, so thread count and host
# RNGs are outside the value path by construction.
set -e
OUT=${TMPDIR:-/tmp}/parity_eqfit
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py/src/arcstat" || exit 1
Rscript -e '
suppressMessages(library(arcstat))
n <- 120
x <- 20*(0:(n-1))/(n-1)
curve <- cbind(-0.64 + 0.12*(0:12)/12,
               -0.90 + 0.85*(0:12)/12)
th <- c(1, 300, 0, 20, -0.60, -0.40)
u  <- (1:n)/(n+1)
y  <- th[1] + th[2]*u + 3*sin(4321*(1:n))
st <- rbind(c(1, log(300), 0, log(20), 0.3, 0.5),
            c(1, log(250), 1, log(15), -0.2, 0.2))
f  <- eqfit_bc(x, y, st, curve, maxit = 400L)
es <- eqfit_estsim(-0.60+1, -0.35+1, R = 4L, nsizes = c(60L, 120L), seed = 11L, maxit = 200L)
bs <- k4_b_sweep(seq(0.90, 1.10, by = 0.05), 0.5)
Es <- eq_E_sweep(-0.60, seq(-0.9, -0.1, by = 0.1))
bl <- c(eqfit_blocklen(sin(0.05*(1:300))), eqfit_blocklen(sin(4321*(1:300))))
## the manifold-test simulation entry points, at a toy configuration: the value of the check
## is byte-agreement between the fronts on the SAME code path that the manifold table uses,
## not statistical meaning
th0 <- c(1, 300, 0, 20, -0.60, -0.40)
thref <- rbind(c(0,20,-0.58,-0.35), c(0,20,-0.60,-0.40), c(0,20,-0.62,-0.55))
msv <- eqfit_msim(60L, 1.0, 2L, th0, c(0,0,0), curve, c(-0.64,-0.52), thref,
                  seed=17L, maxit=200L)
ntv <- eqfit_nullT(60L, 1.0, 2L, th0, c(0,0,0), curve, c(-0.64,-0.52), seed=18L, maxit=200L)
vals <- c(f$th, f$sse, f$a, f$b, as.vector(es), bs, Es, bl, unlist(msv), ntv)
cat(sprintf("%.17g", vals), sep="\n")
' > "$OUT/r.txt"
PYTHONPATH="$PKG/arcstat-py/src" python3 - <<PY > "$OUT/py.txt"
import math
from arcstat import eqfit_bc, eqfit_estsim, k4_b_sweep, eq_E_sweep
n = 120
x = [20*i/(n-1) for i in range(n)]
curve = [(-0.64 + 0.12*i/12, -0.90 + 0.85*i/12) for i in range(13)]
th = [1, 300, 0, 20, -0.60, -0.40]
y = [th[0] + th[1]*((i+1)/(n+1)) + 3*math.sin(4321*(i+1)) for i in range(n)]
st = [[1, math.log(300), 0, math.log(20), 0.3, 0.5],
      [1, math.log(250), 1, math.log(15), -0.2, 0.2]]
thf, sse, a, b = eqfit_bc(x, y, st, curve, maxit=400)
es = eqfit_estsim(0.40, 0.65, 4, [60, 120], seed=11, maxit=200)
bs = k4_b_sweep([0.90 + 0.05*i for i in range(5)], 0.5)
Es = eq_E_sweep(-0.60, [-0.9 + 0.1*i for i in range(9)])
from arcstat import eqfit_blocklen
bl = [eqfit_blocklen([math.sin(0.05*(i+1)) for i in range(300)]),
      eqfit_blocklen([math.sin(4321*(i+1)) for i in range(300)])]
from arcstat import eqfit_msim, eqfit_nullT
th0 = [1, 300, 0, 20, -0.60, -0.40]
thref = [[0,20,-0.58,-0.35], [0,20,-0.60,-0.40], [0,20,-0.62,-0.55]]
msv = eqfit_msim(60, 1.0, 2, th0, [0,0,0], curve, [-0.64,-0.52], thref, seed=17, maxit=200)
ntv = eqfit_nullT(60, 1.0, 2, th0, [0,0,0], curve, [-0.64,-0.52], seed=18, maxit=200)
mflat = [v for k in ["Tlev","Taud","Tp2","Tp15","Tp4","Tref2","Tref15","Tref4"] for v in msv[k]]
vals = list(thf) + [sse, a, b] + es + bs + Es + bl + mflat + ntv
print("\n".join(("NaN" if math.isnan(v) else "%.17g" % v) for v in vals))
PY
if cmp -s "$OUT/r.txt" "$OUT/py.txt"; then
  echo "VERDICT: PARITY -- $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical"
else
  echo "VERDICT: MISMATCH"; diff "$OUT/r.txt" "$OUT/py.txt" | head; exit 1
fi
