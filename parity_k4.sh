#!/bin/sh
# Parity harness for the four-parameter kappa module (arck4.c).
# R and Python compile the SAME C source; every value below must agree byte for byte.
set -e
OUT=${TMPDIR:-/tmp}/parity_k4
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
# refuse to run against a stale copy of the shared back-end (see the guard for why)
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py" || exit 1
Rscript -e '
suppressMessages(library(arcstat))
u  <- (1:9)/10
xg <- seq(2, 12, by=0.5)
th <- c(1, 400, 8, 2.5, 0.2, 0.4)
n  <- 200
x  <- 20*(0:(n-1))/(n-1)
y  <- th[1] + th[2]*k4_cdf(x, th[3], th[4], th[5], th[6]) + 2*sin(9871*(1:n))
br <- 20*(0:12)/12; br[13] <- br[13] + 2e-8
st <- c(1, log(380), 7.5, log(2.2), 0.15, log(0.5))
vals <- c(
  k4_q(u, 0, 1, 0.2, 0.4), k4_cdf(xg, 8, 2.5, 0.2, 0.4), k4_pdf(xg, 8, 2.5, 0.2, 0.4),
  k4_tau34(0.2, 0.4), k4_tau34(-0.15, 0.8),
  k4_runmed(y[1:40], 9L),
  k4_band_model(th, br), k4_band_sample(x, y, br),
  k4_readings(th),
  k4_fit_lmom(0.13737, 0.102328),
  k4_fit_aleq(y, rbind(c(.05,.20),c(.20,.40),c(.40,.60),c(.60,.80),c(.80,.90))),
  k4_fit_nls(x, y, st)$theta, k4_fit_nls(x, y, st)$obj,
  k4_fit_nalr(x, y, st)$theta, k4_fit_nalr(x, y, st)$obj,
  eq_readings(-0.60, -0.3778), eq_readings(-0.70, -0.55, theta=c(0.2)),
  eq_readings_vsl(0.25, 0.5966), eq_readings_vsl(0.35, 0.3459))
cat(sprintf("%.17g\n", vals), sep="")
' > "$OUT/r.txt"

python3 - <<PY > "$OUT/py.txt"
import sys, math
sys.path.insert(0, "$PKG/arcstat-py/src")
from arcstat import _core as c
u  = [i/10 for i in range(1,10)]
xg = [2 + 0.5*i for i in range(21)]
th = [1, 400, 8, 2.5, 0.2, 0.4]
n  = 200
x  = [20*i/(n-1) for i in range(n)]
base = c.k4_cdf(x, 8, 2.5, 0.2, 0.4)
y  = [1 + 400*base[i] + 2*math.sin(9871*(i+1)) for i in range(n)]
br = [20*j/12 for j in range(13)]; br[12] += 2e-8
st = [1, math.log(380), 7.5, math.log(2.2), 0.15, math.log(0.5)]
bands = [.05,.20,.20,.40,.40,.60,.60,.80,.80,.90]
vals = []
vals += c.k4_q(u,0,1,0.2,0.4) + c.k4_cdf(xg,8,2.5,0.2,0.4) + c.k4_pdf(xg,8,2.5,0.2,0.4)
vals += c.k4_tau34(0.2,0.4) + c.k4_tau34(-0.15,0.8)
vals += c.k4_runmed(y[:40], 9)
vals += c.k4_band_model(th, br) + c.k4_band_sample(x, y, br)
vals += c.k4_readings(th)
vals += c.k4_fit_lmom(0.13737, 0.102328)
vals += c.k4_fit_aleq(sorted(y), bands, 5, [0,0,math.log(0.3)])
r = c.k4_fit_nls(x, y, st); vals += r[:6] + [r[6]]
r = c.k4_fit_nalr(x, y, st); vals += r[:6] + [r[6]]
vals += c.eq_readings(-0.60, -0.3778) + c.eq_readings(-0.70, -0.55, theta=[0.2])
vals += c.eq_readings_vsl(0.25, 0.5966) + c.eq_readings_vsl(0.35, 0.3459)
for v in vals: print("%.17g" % v)
PY

if diff -q "$OUT/r.txt" "$OUT/py.txt" >/dev/null; then
  echo "PARITY OK: $(wc -l < "$OUT/r.txt" | tr -d ' ') values byte-identical"
else
  echo "PARITY FAIL:"; diff "$OUT/r.txt" "$OUT/py.txt" | head -10
fi
