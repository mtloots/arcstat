#!/bin/sh
# Byte-parity harness for ellstat: every entry point, both fronts, %.15e.
#
# Grids are built by the SAME arithmetic on both sides (i/(n-1), never seq()),
# because seq() and a computed grid differ in the last bit and that difference
# propagates. Seeds are pushed as doubles and counts as ints on both fronts, so
# the two languages hand the C identical bits.
set -e
cd "$(dirname "$0")"
sh parity_sync_guard.sh ellstat ellstat-py/src/ellstat

R_OUT=$(mktemp); PY_OUT=$(mktemp)
trap 'rm -f "$R_OUT" "$PY_OUT"' EXIT

Rscript -e '
suppressMessages(library(ellstat))
p <- function(lab, v) for (i in seq_along(v)) cat(sprintf("%s[%d] %.15e\n", lab, i, v[i]))
ks <- c(0.1, 0.5, 0.9, 0.99)
ke <- ell_KE(ks); p("K", ke$K); p("E", ke$E); p("nome", ell_nome(ks))
u  <- vapply(0:20, function(i) 4*i/20, 0)
j  <- ell_jacobi(u, 0.7); p("sn", j$sn); p("cn", j$cn); p("dn", j$dn)
t  <- vapply(0:31, function(i) 2*pi*i/32, 0)
b  <- ell_basis(t, 3L, 0.9); p("basC", as.vector(t(b$C))); p("basS", as.vector(t(b$S)))
p("smom", as.vector(t(ell_moments(t, 3L, 0.9))))
p("mvar", as.vector(t(ell_moment_var(t, 3L, 0.9))))
p("tran", as.vector(t(ell_transfer(0.9, 3L, 5L))))
p("dn2",  ell_dn2_moments(0.8, 5L))
p("dvm",  dellvm(t, mu = 1, kappa = 2, k = 0.9, ng = 256L))
p("rvm",  rellvm(12L, mu = 1, kappa = 2, k = 0.9, seed = 4242))
p("pmom", as.vector(t(ell_evm_moments(1, 2, 0.9, 2L, 256L))))
p("fit",  fit_ellvm(rellvm(200L, mu = 1, kappa = 3, k = 0.9, seed = 77), 0.9, 2L, 256L))
uu <- ell_exact_uniform(rellvm(60L, 0, 2, 0.9, seed = 11), 0.9, 199L, seed = 5)
p("unif", c(uu$statistic, uu$p.value))
ci <- ell_exact_ci(rellvm(80L, 1, 3, 0.9, seed = 13), 0.9,
                   kgrid = vapply(0:7, function(i) 1 + 5*i/7, 0),
                   B = 99L, level = 0.05, seed = 17, ng = 256L)
p("ci", c(ci$lower, ci$upper, ci$statistic, ci$curve$p.value))
p("eff", ell_efficiency(0.9, 3, 120L, 200L, 256L, seed = 23))
p("pjd", dprojell(t, tau = 0.3, d = 0.6, beta = 2, ng = 120L))
p("pjr", rprojell(12L, tau = 0.3, d = 0.6, beta = 2, seed = 909))
ec <- ell_eff_curve(0.3, 0.6, 2, kgrid = vapply(0:4, function(i) i/5, 0),
                    nt = 256L, ng = 120L)
p("pje", c(ec$efficiency, ec$fisher, ec$godambe))
p("pjk", ell_kopt(0.3, 0.6, 2, nt = 256L, ng = 120L))
p("pfn", dprojfam(t, tau = 0.6, d = 1.2, family = "normal"))
p("pfc", dprojfam(t, tau = 0.6, d = 1.2, family = "cauchy"))
p("pfr", rprojfam(10L, tau = 0.6, d = 1.2, family = "cauchy", seed = 606))
pe <- ell_projfam_eff(0.6, 1.2, "normal", kgrid = vapply(0:3, function(i) i/4, 0), nt = 256L)
p("pfe", c(pe$efficiency, pe$fisher, pe$mass))
ct <- ell_catalogue("wrappedcauchy", 0.7, kgrid = vapply(0:3, function(i) i/4, 0), nt = 256L)
p("cat", c(ct$efficiency, ct$fisher, ct$mass, ct$kopt, ct$eff_kopt))
fj <- dellvm(vapply(0:255, function(i) 2*pi*(i+0.5)/256, 0), 1, 3, 0.9, 256L)
jj <- ell_joint(fj, 2L, 0.9)
p("jco", c(jj$mean, as.vector(t(jj$Sigma))))
js <- ell_joint_sample(rellvm(50L, 1, 3, 0.9, seed = 4141), 2L, 0.9)
p("jsa", c(js$mean, as.vector(t(js$Sigma))))
ju <- ell_joint_uniform(rellvm(40L, 1, 2, 0.9, seed = 31), 2L, 0.9, 99L, 256L, seed = 5)
p("jun", c(ju$statistic, ju$p.value))
' > "$R_OUT" 2>/dev/null

PYTHONPATH=ellstat-py/src python3 -c '
import math, ellstat as e
def p(lab, v):
    v = v if hasattr(v, "__len__") else [v]
    for i, x in enumerate(v): print("%s[%d] %.15e" % (lab, i+1, x))
ks = [0.1, 0.5, 0.9, 0.99]
K, E = e.ell_KE(ks); p("K", K); p("E", E); p("nome", e.ell_nome(ks))
u = [4*i/20 for i in range(21)]
sn, cn, dn = e.ell_jacobi(u, 0.7); p("sn", sn); p("cn", cn); p("dn", dn)
t = [2*math.pi*i/32 for i in range(32)]
C, S = e.ell_basis(t, 3, 0.9)
p("basC", [x for row in C for x in row]); p("basS", [x for row in S for x in row])
p("smom", [x for pr in e.ell_moments(t, 3, 0.9) for x in pr])
p("mvar", [x for tr in e.ell_moment_var(t, 3, 0.9) for x in tr])
p("tran", [x for row in e.ell_transfer(0.9, 3, 5) for x in row])
p("dn2",  e.ell_dn2_moments(0.8, 5))
p("dvm",  e.dellvm(t, 1, 2, 0.9, 256))
p("rvm",  e.rellvm(12, 1, 2, 0.9, 4242))
p("pmom", [x for pr in e.ell_evm_moments(1, 2, 0.9, 2, 256) for x in pr])
p("fit",  list(e.fit_ellvm(e.rellvm(200, 1, 3, 0.9, 77), 0.9, 2, 256)))
uu = e.ell_exact_uniform(e.rellvm(60, 0, 2, 0.9, 11), 0.9, 199, 5)
p("unif", [uu["statistic"], uu["p_value"]])
ci = e.ell_exact_ci(e.rellvm(80, 1, 3, 0.9, 13), 0.9,
                    [1 + 5*i/7 for i in range(8)], 99, 0.05, 17, 256)
p("ci", [ci["lower"], ci["upper"], ci["statistic"]] + [q for _, q in ci["curve"]])
p("eff", list(e.ell_efficiency(0.9, 3, 120, 200, 256, 23)))
p("pjd", e.dprojell(t, 0.3, 0.6, 2, 120))
p("pjr", e.rprojell(12, 0.3, 0.6, 2, 909))
ec = e.ell_eff_curve(0.3, 0.6, 2, [i/5 for i in range(5)], 256, 120)
p("pje", ec["efficiency"] + [ec["fisher"], ec["godambe"]])
p("pjk", list(e.ell_kopt(0.3, 0.6, 2, 256, 120)))
p("pfn", e.dprojfam(t, 0.6, 1.2, "normal"))
p("pfc", e.dprojfam(t, 0.6, 1.2, "cauchy"))
p("pfr", e.rprojfam(10, 0.6, 1.2, "cauchy", 606))
pe = e.ell_projfam_eff(0.6, 1.2, "normal", [i/4 for i in range(4)], 256)
p("pfe", pe["efficiency"] + [pe["fisher"], pe["mass"]])
ct = e.ell_catalogue("wrappedcauchy", 0.7, 1.0, 1, [i/4 for i in range(4)], 256)
p("cat", ct["efficiency"] + [ct["fisher"], ct["mass"], ct["kopt"], ct["eff_kopt"]])
fj = e.dellvm([2*math.pi*(i+0.5)/256 for i in range(256)], 1, 3, 0.9, 256)
jj = e.ell_joint(fj, 2, 0.9)
p("jco", jj["mean"] + [x for row in jj["Sigma"] for x in row])
js = e.ell_joint_sample(e.rellvm(50, 1, 3, 0.9, 4141), 2, 0.9)
p("jsa", js["mean"] + [x for row in js["Sigma"] for x in row])
ju = e.ell_joint_uniform(e.rellvm(40, 1, 2, 0.9, 31), 2, 0.9, 99, 256, 5)
p("jun", [ju["statistic"], ju["p_value"]])
' > "$PY_OUT" 2>/dev/null

NR=$(wc -l < "$R_OUT"); NP=$(wc -l < "$PY_OUT")
if [ "$NR" -ne "$NP" ]; then
  echo "VERDICT: FAIL -- $NR values from R against $NP from Python"; exit 1
fi
if diff -q "$R_OUT" "$PY_OUT" > /dev/null; then
  echo "VERDICT: PASS -- $NR values byte-identical across the R and Python fronts"
else
  echo "VERDICT: FAIL -- differing values:"; diff "$R_OUT" "$PY_OUT" | head -20; exit 1
fi
