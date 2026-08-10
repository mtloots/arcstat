#!/bin/sh
# Parity harness for the circular (wrapped) arc-length family.
#
# The R package arcstat and the Python package arcstat-py compile the SAME C source, arccirc.c, so
# every value below must agree to the last printed digit. The comparison is byte-identical: the two
# sides print with the same format and the files are diffed, not compared numerically. Simulation is
# driven by uniforms supplied from the caller, so the C routine is deterministic and its output is
# comparable across the two front ends without seeding two different generators.
#
# Usage: sh parity_arccirc.sh
set -e
OUT=${TMPDIR:-/tmp}/parity_arccirc
mkdir -p "$OUT"
PKG="$(cd "$(dirname "$0")" && pwd)"
# refuse to run against a stale copy of the shared back-end (see the guard for why)
sh "$PKG/parity_sync_guard.sh" "$PKG/arcstat/src" "$PKG/arcstat-py" || exit 1
# ---- deterministic inputs, shared by both sides --------------------------------------------------
# u grid, angles, shape values, and the uniforms used for simulation
Rscript -e '
u   <- (0:10)/10
th  <- (0:8) * (2*pi/8)
c3s <- c(-1.0392304845, -0.6, -0.2, 0, 0.3, 0.75, 1.0392304845)
unf <- (seq_len(20) - 0.5)/20
cat(sprintf("%.17g\n", c(u, th, c3s, unf)), sep="")
' > "$OUT/inputs.txt"

# ---- R side --------------------------------------------------------------------------------------
Rscript -e '
suppressMessages(library(arcstat))
u   <- (0:10)/10
th  <- (0:8) * (2*pi/8)
c3s <- c(-1.0392304845, -0.6, -0.2, 0, 0.3, 0.75, 1.0392304845)
unf <- (seq_len(20) - 0.5)/20
p   <- function(x) cat(sprintf("%.17g\n", x), sep="")

p(arcc_c3max())
for (c3 in c3s) {
  o <- arccirc(c3 = c3, mu = 0.7)
  p(qd_arccirc(u, o))
  p(Q_arccirc(u, o))
  p(darccirc(th, o))
  for (k in 1:4) { z <- trigmom_arccirc(k, o, nodes = 2048L); p(c(Re(z), Im(z))) }
  p(.C("arcc_rand3", unif = as.double(unf), n = 20L, c3 = as.double(c3),
       mu = as.double(0.7), out = double(20), PACKAGE = "arcstat")$out)
}
# method-of-moments fit on a deterministic sample
for (c3 in c(-0.7, 0.4)) {
  uu <- (seq_len(4000) - 0.5)/4000
  ang <- .C("arcc_rand3", unif = as.double(uu), n = 4000L, c3 = as.double(c3),
            mu = as.double(1.1), out = double(4000), PACKAGE = "arcstat")$out
  f <- fit_arccirc(ang, nodes = 2048L)
  p(c(f$mu, f$c3, f$rho, as.numeric(f$interior)))
}
# Fejer-Riesz form
pre <- c(1.0, 0.4, -0.3); pim <- c(0.0, 0.25, 0.1)
o <- arccirc_fr(complex(real = pre, imaginary = pim), mu = 0.35)
p(qd_arccirc_fr(u, o))
for (k in 1:3) { z <- trigmom_arccirc_fr(k, o, nodes = 2048L); p(c(Re(z), Im(z))) }

# Fejer-Riesz block: quantile function, density, simulation, circular L-moments, fit
pre2 <- c(1.0, 0.4, -0.3); pim2 <- c(0.0, 0.25, 0.1)
o2 <- arccirc_fr(complex(real = pre2, imaginary = pim2), mu = 0.35)
p(Q_arccirc_fr(u, o2, nodes = 1024L))
p(darccirc_fr(th, o2, nodes = 1024L))
p(.C("arcc_rand_fr", unif = as.double(unf), n = 20L, pre = as.double(pre2), pim = as.double(pim2),
     np = 3L, mu = as.double(0.35), nodes = 1024L, out = double(20), PACKAGE = "arcstat")$out)
xs <- (seq_len(500) - 0.5)/500
el <- lmom_circ(xs, 3L);        p(Re(el)); p(Im(el))
rr <- rho_from_lmom(el);        p(Re(rr)); p(Im(rr))
p(qmin_rho(rr, nodes = 1024L))
ad <- admiss_rho(rr, margin = 0, nodes = 1024L); p(Re(ad$rho)); p(Im(ad$rho)); p(ad$shrink)
fz <- factorise_rho(rr);        p(Re(fz)); p(Im(fz))
uu2 <- (seq_len(3000) - 0.5)/3000
ang2 <- .C("arcc_rand_fr", unif = as.double(uu2), n = 3000L, pre = as.double(pre2),
           pim = as.double(pim2), np = 3L, mu = as.double(0.9), nodes = 1024L,
           out = double(3000), PACKAGE = "arcstat")$out
ff <- fit_arccirc_fr(ang2, nm = 3L, nodes = 1024L)
p(Re(ff$p)); p(Im(ff$p)); p(c(ff$mu, ff$shrink, as.numeric(ff$admissible)))
' > "$OUT/r.txt"

# ---- Python side -----------------------------------------------------------------------------------
PYTHONPATH="$PKG/arcstat-py/src" python3 -c '
import arcstat._core as c
import math
u   = [i/10 for i in range(11)]
th  = [i*(2*math.pi/8) for i in range(9)]
c3s = [-1.0392304845, -0.6, -0.2, 0, 0.3, 0.75, 1.0392304845]
unf = [(i+0.5)/20 for i in range(20)]
def p(xs):
    for x in xs: print("%.17g" % x)

p([c.arcc_c3max()])
for c3 in c3s:
    p(c.arcc_qd3(u, c3))
    p(c.arcc_Q3(u, c3))
    p(c.arcc_dens3(th, c3, 0.7))
    for k in (1,2,3,4):
        z = c.arcc_trigmom3(k, c3, 0.7, 2048)
        p([z.real, z.imag])
    p(c.arcc_rand3(unf, c3, 0.7))
for c3 in (-0.7, 0.4):
    uu = [(i+0.5)/4000 for i in range(4000)]
    ang = c.arcc_rand3(uu, c3, 1.1)
    f = c.arcc_fit3(ang, 2048)
    p([f["mu"], f["c3"], f["rho"], 1.0 if f["interior"] else 0.0])
pre = [1.0, 0.4, -0.3]; pim = [0.0, 0.25, 0.1]
p(c.arcc_qd_fr(u, pre, pim))
for k in (1,2,3):
    z = c.arcc_trigmom_fr(k, pre, pim, 0.35, 2048)
    p([z.real, z.imag])

pre2 = [1.0, 0.4, -0.3]; pim2 = [0.0, 0.25, 0.1]
p(c.arcc_Q_fr(u, pre2, pim2, 1024))
p(c.arcc_dens_fr(th, pre2, pim2, 0.35, 1024))
p(c.arcc_rand_fr(unf, pre2, pim2, 0.35, 1024))
xs = [(i+0.5)/500 for i in range(500)]
el = c.arcc_lmom(xs, 3);        p([z.real for z in el]); p([z.imag for z in el])
rr = c.arcc_rho_from_lmom(el);  p([z.real for z in rr]); p([z.imag for z in rr])
p([c.arcc_qmin_rho(rr, 1024)])
adr, adt = c.arcc_admiss(rr, 0.0, 1024)
p([z.real for z in adr]); p([z.imag for z in adr]); p([adt])
fz = c.arcc_factorise(rr);      p([z.real for z in fz]); p([z.imag for z in fz])
uu2 = [(i+0.5)/3000 for i in range(3000)]
ang2 = c.arcc_rand_fr(uu2, pre2, pim2, 0.9, 1024)
ff = c.arcc_fit_fr(ang2, 3, 0.0, 1024)
p([z.real for z in ff["p"]]); p([z.imag for z in ff["p"]])
p([ff["mu"], ff["shrink"], 1.0 if ff["admissible"] else 0.0])
' > "$OUT/py.txt"

# ---- compare ---------------------------------------------------------------------------------------
NR=$(wc -l < "$OUT/r.txt" | tr -d ' ')
NP=$(wc -l < "$OUT/py.txt" | tr -d ' ')
echo "values compared: R $NR, Python $NP"
if [ "$NR" != "$NP" ]; then
  echo "PARITY FAIL: different number of values"
  exit 1
fi
if diff -q "$OUT/r.txt" "$OUT/py.txt" > /dev/null; then
  echo "PARITY OK: $NR values byte-identical between the R and Python front ends"
else
  echo "PARITY FAIL: first differing values"
  diff "$OUT/r.txt" "$OUT/py.txt" | head -20
  exit 1
fi
