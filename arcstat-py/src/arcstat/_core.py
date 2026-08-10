"""arcstat: Python (ctypes) binding to the shared arc-length C back-end. Pure standard library.
The three bundled C sources (goodness of fit, distributions, Bayesian test) are compiled together on
first import; the same sources back the R package 'arcstat'."""
import ctypes as _ct, os as _os
_here = _os.path.dirname(_os.path.abspath(__file__))
_ext = "dylib" if _os.uname().sysname == "Darwin" else "so"
_libpath = _os.path.join(_here, "libarcstat." + _ext)
if not _os.path.exists(_libpath):                       # compile the three cores together on first import
    import subprocess as _sp
    _sp.check_call(["cc", "-O2", "-fPIC", "-shared", "-I", _here,
                    _os.path.join(_here, "arclen.c"), _os.path.join(_here, "arcdistc.c"),
                    _os.path.join(_here, "bayesarc.c"), _os.path.join(_here, "cfarc.c"), _os.path.join(_here, "arccirc.c"),
                    _os.path.join(_here, "arck4.c"), _os.path.join(_here, "arceq.c"),
                    _os.path.join(_here, "arceq2.c"),
                    "-lm", "-o", _libpath])
_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)
_lib.al_statistic.argtypes = [_ip, _dp, _dp]
_lib.al_pvalue.argtypes = [_dp, _ip, _dp]
_lib.al_moments.argtypes = [_ip, _dp]
_lib.arcq_qd.argtypes = [_dp, _ip, _dp, _ip, _dp, _dp]
_lib.arcq_arclength.argtypes = [_dp, _ip, _dp, _ip, _dp]
_lib.al_band_model.argtypes = [_dp, _dp, _dp, _ip, _dp]
_lib.al_band_sample.argtypes = [_dp, _ip, _dp, _dp, _dp]
_lib.al_scale.argtypes = [_dp, _ip, _dp, _dp, _dp]
_lib.sample_lmoments_c.argtypes = [_dp, _ip, _ip, _dp]
_lib.arcq_fit_cf_c.argtypes = [_dp, _ip, _dp]
_lib.bb_post.argtypes = [_ip, _dp, _ip, _ip, _dp]
_lib.bb_ref.argtypes = [_ip, _ip, _ip, _ip, _dp]
_lib.cf_arclength_emp.argtypes = [_dp, _ip, _dp, _ip, _dp]
_lib.arcc_c3max.argtypes = [_dp]
_lib.arcc_qd3.argtypes = [_dp, _ip, _dp, _dp]
_lib.arcc_Q3.argtypes = [_dp, _ip, _dp, _dp]
_lib.arcc_dens3.argtypes = [_dp, _ip, _dp, _dp, _dp]
_lib.arcc_trigmom3.argtypes = [_ip, _dp, _dp, _ip, _dp]
_lib.arcc_rand3.argtypes = [_dp, _ip, _dp, _dp, _dp]
_lib.arcc_fit3.argtypes = [_dp, _ip, _ip, _dp]
_lib.arcc_qd_fr.argtypes = [_dp, _ip, _dp, _dp, _ip, _dp]
_lib.arcc_trigmom_fr.argtypes = [_ip, _dp, _dp, _ip, _dp, _ip, _dp]
_lib.arcc_Q_fr.argtypes = [_dp, _ip, _dp, _dp, _ip, _ip, _dp]
_lib.arcc_dens_fr.argtypes = [_dp, _ip, _dp, _dp, _ip, _dp, _ip, _dp]
_lib.arcc_rand_fr.argtypes = [_dp, _ip, _dp, _dp, _ip, _dp, _ip, _dp]
_lib.arcc_lmom.argtypes = [_dp, _ip, _ip, _dp]
_lib.arcc_rho_from_lmom.argtypes = [_dp, _dp, _ip, _dp]
_lib.arcc_qmin_rho.argtypes = [_dp, _dp, _ip, _ip, _dp]
_lib.arcc_admiss.argtypes = [_dp, _dp, _ip, _dp, _ip, _dp]
_lib.arcc_factorise.argtypes = [_dp, _dp, _ip, _dp]
_lib.arcc_fit_fr.argtypes = [_dp, _ip, _ip, _dp, _ip, _dp]
_lib.arcc_temper_vm.argtypes = [_dp, _ip, _dp, _dp, _dp, _ip, _dp]
_lib.arcc_temper_vm_trigmom.argtypes = [_ip, _dp, _dp, _dp, _ip, _dp]

def _d(seq):
    a = (_ct.c_double * max(len(seq), 1))()
    for i, v in enumerate(seq): a[i] = v
    return a
def _i(v): a = (_ct.c_int * 1)(); a[0] = v; return a
def _d1(v): a = (_ct.c_double * 1)(); a[0] = v; return a

# ---- goodness of fit ----
def al_statistic(u):
    """Arc-length goodness-of-fit statistic (C back-end)"""
    out = _d1(0.0); _lib.al_statistic(_i(len(u)), _d(u), out); return out[0]
def al_pvalue(s, n):
    """Saddlepoint right-tail probability of the arc-length statistic (C back-end)"""
    out = _d1(0.0); _lib.al_pvalue(_d1(s), _i(n), out); return out[0]
def al_moments(n):
    """Exact mean, variance and support of the arc-length statistic (C back-end)"""
    out = (_ct.c_double * 4)(); _lib.al_moments(_i(n), out)
    return {"mean": out[0], "var": out[1], "support": (out[2], out[3])}

# ---- distributions ----
def arcq_qd(u, coef=(), sigma=1.0):
    """Quantile density of the arcq family. Binding to the C routine arcq_qd."""
    n = len(u); out = (_ct.c_double * n)()
    _lib.arcq_qd(_d(u), _i(n), _d(coef), _i(len(coef)), _d1(sigma), out); return list(out)
def arcq_arclength(coef=(), sigma=1.0, nodes=24):
    """Arc length of the quantile function of the arcq family, by Gauss-Legendre quadrature.
    A shape functional: shifting the distribution cannot change it. Binding to the C routine
    arcq_arclength; the R front calls this arclength()."""
    out = _d1(0.0); _lib.arcq_arclength(_d(coef), _i(len(coef)), _d1(sigma), _i(nodes), out); return out[0]
def sample_lmoments(x, nmom=4):
    """Sample L-moments of data"""
    out = (_ct.c_double * nmom)(); _lib.sample_lmoments_c(_d(x), _i(len(x)), _i(nmom), out); return list(out)


def fit_arcq_cf(x):
    """Closed-form order-two L-moment fit of the arcq family.

    The family's L-moments are an exact linear function of its shape coefficients, so matching
    inverts explicitly and no numerical search is involved.  The inversion is derived on the
    admissible set only, so the last return value flags whether the estimate stayed inside it.
    Returns (c1, c2, mu, sigma, admissible).
    """
    out = (_ct.c_double * 5)()
    _lib.arcq_fit_cf_c(_d(x), _i(len(x)), out)
    return out[0], out[1], out[2], out[3], out[4] > 0.5

# ---- Bayesian test ----
def bb_post_disc(u, M, seed=1):
    """Bayesian-bootstrap posterior arc-length discrepancies"""
    u = sorted(u); out = (_ct.c_double * (2 * M))()
    _lib.bb_post(_i(len(u)), _d(u), _i(M), _i(seed), out)
    return [(out[2 * i], out[2 * i + 1]) for i in range(M)]
def bb_ref_disc(n, D=300, m=16, seed=7):
    """Null reference distribution of the arc-length discrepancy"""
    out = (_ct.c_double * (2 * D * m))(); _lib.bb_ref(_i(n), _i(D), _i(m), _i(seed), out)
    return [(out[2 * i], out[2 * i + 1]) for i in range(D * m)]
def bb_evidence(u, which="arc", level=0.95, M=1500, ref=None, seed=1):
    """Bayesian arc-length goodness-of-fit evidence"""
    j = 0 if which == "arc" else 1
    if ref is None: ref = bb_ref_disc(len(u))
    col = sorted(r[j] for r in ref); thr = col[int(level * (len(col) - 1))]
    return sum(1 for r in bb_post_disc(u, M, seed) if r[j] > thr) / M

# ---- characteristic-function arc length ----
import math as _math
def cf_arclength(x, T=6.0, ngrid=1200):
    """Windowed arc length of the empirical characteristic function of a sample (standardised
    internally, so scale free). A symmetric law tends to two; asymmetry adds length."""
    xf = [float(v) for v in x if v == v and abs(v) != _math.inf]
    out = _d1(0.0)
    _lib.cf_arclength_emp(_d(xf), _i(len(xf)), _d1(float(T)), _i(int(ngrid)), out)
    return out[0]

def cf_arclength_family(family="exponential", k=1.0, delta=0.0, alpha=1.5, beta=0.0, lam=1.0):
    """Closed-form total arc length of the characteristic-function curve for a named family:
    'exponential', 'gamma' (shape k), 'normal' (drift ratio delta), 'cauchy' (delta),
    'skewstable' (alpha, beta), 'poisson' (rate lam, per period). Scale free; a symmetric monotone
    law carries two. The normal and Poisson forms use modified Bessel functions."""
    from math import pi, sqrt, exp, tan, gamma
    def _besselI0(z):
        s, term, k2 = 1.0, 1.0, 0
        while True:
            k2 += 1; term *= (z / 2) ** 2 / (k2 * k2)
            s += term
            if term < 1e-16 * s: break
        return s
    def _besselK01(z):
        # K0, K1 by numerical integration of the standard integral representations
        n = 2000; hi = 40.0 / max(z, 0.1); h = hi / n
        k0 = k1 = 0.0
        for i in range(n + 1):
            t = i * h; w = 0.5 if (i == 0 or i == n) else 1.0
            k0 += w * exp(-z * _math.cosh(t)) * h
            k1 += w * exp(-z * _math.cosh(t)) * _math.cosh(t) * h
        return k0, k1
    if family == "exponential": return pi
    if family == "gamma": return 2 * sqrt(pi) * gamma(k / 2 + 1) / gamma((k + 1) / 2)
    if family == "normal":
        x = delta * delta / 4
        if x == 0: return 2.0
        k0, k1 = _besselK01(x); return 2 * x * exp(x) * (k0 + k1)
    if family == "cauchy": return 2 * sqrt(1 + delta * delta)
    if family == "skewstable": return 2 * sqrt(1 + (beta * tan(pi * alpha / 2)) ** 2)
    if family == "poisson": return 2 * pi * lam * exp(-lam) * _besselI0(lam)
    raise ValueError("unknown family")


def al_band_model(sigma, a=0.05, b=0.95, nodes=400):
    """Band arc length of the normal reference curve."""
    out = _d1(0.0); _lib.al_band_model(_d1(sigma), _d1(a), _d1(b), _i(nodes), out); return out[0]

def al_band_sample(x, a=0.05, b=0.95):
    """Sample band arc length of the empirical distribution curve."""
    out = _d1(0.0); _lib.al_band_sample(_d(x), _i(len(x)), _d1(a), _d1(b), out); return out[0]

def al_scale(x, a=0.05, b=0.95):
    """Scale by arc-length band matching, in the scale-equivariant standardised form."""
    out = _d1(0.0); _lib.al_scale(_d(x), _i(len(x)), _d1(a), _d1(b), out); return out[0]


# ---- circular (wrapped) arc-length family ----
# The linear family has Q(0)=0 and Q(1)=1 identically, so theta = 2 pi Q(U) covers the circle exactly
# once: the wrapping series has one term. Order three is the first smooth non-uniform member.
def arcc_c3max():
    """Largest admissible shape parameter of the smooth circular family"""
    out = _d1(0.0); _lib.arcc_c3max(out); return out[0]

def arcc_qd3(u, c3):
    """The smooth circular arc-length family. Binding to the C routine arcc_qd3; the R front
    calls this qd_arccirc()."""
    n = len(u); out = (_ct.c_double * max(n,1))()
    _lib.arcc_qd3(_d(u), _i(n), _d1(c3), out); return list(out)[:n]

def arcc_Q3(u, c3):
    """The smooth circular arc-length family. Binding to the C routine arcc_Q3; the R front
    calls this Q_arccirc()."""
    n = len(u); out = (_ct.c_double * max(n,1))()
    _lib.arcc_Q3(_d(u), _i(n), _d1(c3), out); return list(out)[:n]

def arcc_dens3(theta, c3, mu=0.0):
    """The smooth circular arc-length family. Binding to the C routine arcc_dens3; the R front
    calls this darccirc()."""
    n = len(theta); out = (_ct.c_double * max(n,1))()
    _lib.arcc_dens3(_d(theta), _i(n), _d1(c3), _d1(mu), out); return list(out)[:n]

def arcc_trigmom3(p, c3, mu=0.0, nodes=4096):
    """The smooth circular arc-length family. Binding to the C routine arcc_trigmom3; the R
    front calls this trigmom_arccirc()."""
    out = (_ct.c_double * 2)()
    _lib.arcc_trigmom3(_i(p), _d1(c3), _d1(mu), _i(nodes), out)
    return complex(out[0], out[1])

def arcc_rho3(c3, mu=0.0, nodes=4096):
    """Mean resultant length of the circular arc-length family, the modulus of its first
    trigonometric moment."""
    return abs(arcc_trigmom3(1, c3, mu, nodes))

def arcc_rand3(unif, c3, mu=0.0):
    """The smooth circular arc-length family. Binding to the C routine arcc_rand3; the R front
    calls this rarccirc()."""
    n = len(unif); out = (_ct.c_double * max(n,1))()
    _lib.arcc_rand3(_d(unif), _i(n), _d1(c3), _d1(mu), out); return list(out)[:n]

def arcc_fit3(theta, nodes=4096):
    """Method-of-moments fit of the smooth circular arc-length family. Binding to the C routine
    arcc_fit3; the R front calls this fit_arccirc()."""
    out = (_ct.c_double * 4)()
    _lib.arcc_fit3(_d(theta), _i(len(theta)), _i(nodes), out)
    return {"mu": out[0], "c3": out[1], "rho": out[2], "interior": out[3] == 1.0}

def arcc_qd_fr(u, pre, pim):
    """Fejer-Riesz form of the circular arc-length family. Binding to the C routine arcc_qd_fr;
    the R front calls this qd_arccirc_fr()."""
    n = len(u); out = (_ct.c_double * max(n,1))()
    _lib.arcc_qd_fr(_d(u), _i(n), _d(pre), _d(pim), _i(len(pre)), out); return list(out)[:n]

def arcc_trigmom_fr(p, pre, pim, mu=0.0, nodes=4096):
    """Fejer-Riesz form of the circular arc-length family. Binding to the C routine
    arcc_trigmom_fr; the R front calls this trigmom_arccirc_fr()."""
    out = (_ct.c_double * 2)()
    _lib.arcc_trigmom_fr(_i(p), _d(pre), _d(pim), _i(len(pre)), _d1(mu), _i(nodes), out)
    return complex(out[0], out[1])


def arcc_Q_fr(u, pre, pim, nodes=4096):
    """Fejer-Riesz form of the circular arc-length family. Binding to the C routine arcc_Q_fr;
    the R front calls this Q_arccirc_fr()."""
    n = len(u); out = (_ct.c_double * max(n,1))()
    _lib.arcc_Q_fr(_d(u), _i(n), _d(pre), _d(pim), _i(len(pre)), _i(nodes), out); return list(out)[:n]

def arcc_dens_fr(theta, pre, pim, mu=0.0, nodes=4096):
    """Fejer-Riesz form of the circular arc-length family. Binding to the C routine
    arcc_dens_fr; the R front calls this darccirc_fr()."""
    n = len(theta); out = (_ct.c_double * max(n,1))()
    _lib.arcc_dens_fr(_d(theta), _i(n), _d(pre), _d(pim), _i(len(pre)), _d1(mu), _i(nodes), out)
    return list(out)[:n]

def arcc_rand_fr(unif, pre, pim, mu=0.0, nodes=4096):
    """Fejer-Riesz form of the circular arc-length family. Binding to the C routine
    arcc_rand_fr; the R front calls this rarccirc_fr()."""
    n = len(unif); out = (_ct.c_double * max(n,1))()
    _lib.arcc_rand_fr(_d(unif), _i(n), _d(pre), _d(pim), _i(len(pre)), _d1(mu), _i(nodes), out)
    return list(out)[:n]

def arcc_lmom(x, nm):
    """Circular L-moments. Binding to the C routine arcc_lmom; the R front calls this
    lmom_circ()."""
    out = (_ct.c_double * (2*nm))()
    _lib.arcc_lmom(_d(x), _i(len(x)), _i(nm), out)
    return [complex(out[2*m], out[2*m+1]) for m in range(nm)]

def arcc_rho_from_lmom(ell):
    """Circular L-moments. Binding to the C routine arcc_rho_from_lmom; the R front calls this
    rho_from_lmom()."""
    nm = len(ell); out = (_ct.c_double * (2*nm))()
    _lib.arcc_rho_from_lmom(_d([e.real for e in ell]), _d([e.imag for e in ell]), _i(nm), out)
    return [complex(out[2*m], out[2*m+1]) for m in range(nm)]

def arcc_qmin_rho(rho, nodes=4096):
    """Circular L-moments. Binding to the C routine arcc_qmin_rho; the R front calls this
    qmin_rho()."""
    out = _d1(0.0)
    _lib.arcc_qmin_rho(_d([r.real for r in rho]), _d([r.imag for r in rho]), _i(len(rho)),
                       _i(nodes), out)
    return out[0]

def arcc_admiss(rho, margin=0.0, nodes=4096):
    """Circular L-moments. Binding to the C routine arcc_admiss; the R front calls this
    admiss_rho()."""
    nm = len(rho); out = (_ct.c_double * (2*nm+1))()
    _lib.arcc_admiss(_d([r.real for r in rho]), _d([r.imag for r in rho]), _i(nm),
                     _d1(margin), _i(nodes), out)
    return [complex(out[2*m], out[2*m+1]) for m in range(nm)], out[2*nm]

def arcc_factorise(rho):
    """Circular L-moments. Binding to the C routine arcc_factorise; the R front calls this
    factorise_rho()."""
    nm = len(rho); out = (_ct.c_double * (2*(nm+1)))()
    _lib.arcc_factorise(_d([r.real for r in rho]), _d([r.imag for r in rho]), _i(nm), out)
    return [complex(out[2*j], out[2*j+1]) for j in range(nm+1)]

def arcc_fit_fr(theta, nm, margin=0.0, nodes=4096):
    """Closed-form fit of the Fejer-Riesz circular arc-length family. Binding to the C routine
    arcc_fit_fr; the R front calls this fit_arccirc_fr()."""
    out = (_ct.c_double * (2*(nm+1)+3))()
    _lib.arcc_fit_fr(_d(theta), _i(len(theta)), _i(nm), _d1(margin), _i(nodes), out)
    p = [complex(out[2*j], out[2*j+1]) for j in range(nm+1)]
    b = 2*(nm+1)
    return {"p": p, "mu": out[b], "shrink": out[b+1], "admissible": out[b+2] == 1.0}


def arcc_temper_vm(theta, kappa, mu=0.0, s=1.0, nodes=4096):
    """Arc-length tempering of a von Mises base. Binding to the C routine arcc_temper_vm; the R
    front calls this dtemper_vm()."""
    n=len(theta); out=(_ct.c_double*max(n,1))()
    _lib.arcc_temper_vm(_d(theta), _i(n), _d1(kappa), _d1(mu), _d1(s), _i(nodes), out)
    return list(out)[:n]

def arcc_temper_vm_trigmom(p, kappa, mu=0.0, s=1.0, nodes=4096):
    """Arc-length tempering of a von Mises base. Binding to the C routine
    arcc_temper_vm_trigmom; the R front calls this trigmom_temper_vm()."""
    out=(_ct.c_double*2)()
    _lib.arcc_temper_vm_trigmom(_i(p), _d1(kappa), _d1(mu), _d1(s), _i(nodes), out)
    return complex(out[0], out[1])


# ---- four-parameter kappa module (induction-curve paper) ----------------------------------------
_lib.arck4_q.argtypes          = [_dp, _ip, _dp, _dp, _dp, _dp, _dp]
_lib.arck4_cdf.argtypes        = [_dp, _ip, _dp, _dp, _dp, _dp, _dp]
_lib.arck4_pdf.argtypes        = [_dp, _ip, _dp, _dp, _dp, _dp, _dp]
_lib.arck4_tau34.argtypes      = [_dp, _dp, _ip, _dp]
_lib.arck4_runmed.argtypes     = [_dp, _ip, _ip, _dp]
_lib.arck4_band_model.argtypes = [_dp, _dp, _ip, _ip, _dp]
_lib.arck4_band_sample.argtypes= [_dp, _dp, _ip, _dp, _ip, _dp]
_lib.arck4_readings.argtypes   = [_dp, _ip, _dp]
_lib.arck4_fit_lmom.argtypes   = [_dp, _dp, _ip, _dp]
_lib.arck4_fit_aleq.argtypes   = [_dp, _ip, _dp, _ip, _dp, _dp]
_lib.arck4_fit_nls.argtypes    = [_dp, _dp, _ip, _dp, _dp]
_lib.arck4_fit_nalr.argtypes   = [_dp, _dp, _ip, _ip, _dp, _ip, _dp, _dp]

def _vec(v): return (_ct.c_double * len(v))(*[float(t) for t in v])

def k4_q(u, mu, sg, k, h):
    """Four-parameter kappa quantile, distribution and density functions"""
    u = list(u); out = (_ct.c_double * len(u))()
    _lib.arck4_q(_vec(u), _i(len(u)), _d1(mu), _d1(sg), _d1(k), _d1(h), out); return list(out)

def k4_cdf(x, mu, sg, k, h):
    """Four-parameter kappa quantile, distribution and density functions"""
    x = list(x); out = (_ct.c_double * len(x))()
    _lib.arck4_cdf(_vec(x), _i(len(x)), _d1(mu), _d1(sg), _d1(k), _d1(h), out); return list(out)

def k4_pdf(x, mu, sg, k, h):
    """Four-parameter kappa quantile, distribution and density functions"""
    x = list(x); out = (_ct.c_double * len(x))()
    _lib.arck4_pdf(_vec(x), _i(len(x)), _d1(mu), _d1(sg), _d1(k), _d1(h), out); return list(out)

def k4_tau34(k, h, nodes=200):
    """Theoretical L-moment ratios of the standard kappa distribution"""
    out = (_ct.c_double * 4)(); _lib.arck4_tau34(_d1(k), _d1(h), _i(nodes), out); return list(out)

def k4_runmed(y, w=9):
    """Running median with shrinking symmetric windows at the edges"""
    y = list(y); out = (_ct.c_double * len(y))()
    _lib.arck4_runmed(_vec(y), _i(len(y)), _i(w), out); return list(out)

def k4_band_model(theta, breaks, nodes=60):
    """Band arc lengths of a scaled kappa curve and of a data polyline"""
    J = len(breaks) - 1; out = (_ct.c_double * J)()
    _lib.arck4_band_model(_vec(theta), _vec(breaks), _i(J), _i(nodes), out); return list(out)

def k4_band_sample(x, y, breaks):
    """Band arc lengths of a scaled kappa curve and of a data polyline"""
    J = len(breaks) - 1; out = (_ct.c_double * J)()
    _lib.arck4_band_sample(_vec(x), _vec(y), _i(len(x)), _vec(breaks), _i(J), out); return list(out)

def k4_readings(theta, grid=4000):
    """The two standard induction-period readings of a fitted kappa curve"""
    out = (_ct.c_double * 3)(); _lib.arck4_readings(_vec(theta), _i(grid), out); return list(out)

def k4_fit_lmom(t3, t4, nodes=200):
    """Deterministic kappa fits: L-moment shape inversion, quantile-domain arc-length shape fit,
    curve-domain NLS and NALR (banded arc lengths of the running-median presmoothed polyline)"""
    out = (_ct.c_double * 3)(); _lib.arck4_fit_lmom(_d1(t3), _d1(t4), _i(nodes), out); return list(out)

def k4_fit_aleq(y_sorted, bands_flat, J, start):
    """Deterministic kappa fits: L-moment shape inversion, quantile-domain arc-length shape fit,
    curve-domain NLS and NALR (banded arc lengths of the running-median presmoothed polyline)"""
    out = (_ct.c_double * 4)()
    _lib.arck4_fit_aleq(_vec(y_sorted), _i(len(y_sorted)), _vec(bands_flat), _i(J), _vec(start), out)
    return list(out)

def k4_fit_nls(x, y, start):
    """Deterministic kappa fits: L-moment shape inversion, quantile-domain arc-length shape fit,
    curve-domain NLS and NALR (banded arc lengths of the running-median presmoothed polyline)"""
    out = (_ct.c_double * 7)()
    _lib.arck4_fit_nls(_vec(x), _vec(y), _i(len(x)), _vec(start), out); return list(out)

def k4_fit_nalr(x, y, start, J=12, lam=1.0, w=9):
    """Deterministic kappa fits: L-moment shape inversion, quantile-domain arc-length shape fit,
    curve-domain NLS and NALR (banded arc lengths of the running-median presmoothed polyline)"""
    out = (_ct.c_double * 7)()
    _lib.arck4_fit_nalr(_vec(x), _vec(y), _i(len(x)), _i(J), _d1(lam), _i(w), _vec(start), out)
    return list(out)

_lib.arceq_readings.argtypes = [_dp, _dp, _ip, _dp, _ip, _dp]
def eq_readings(alpha, beta, theta=(), ngrid=20001):
    """Quantile-domain induction readings and equivalence discrepancy for tilted beta-kernel
    quantile densities q(u) = u^alpha (1-u)^beta exp(sum theta_j P_j(u))"""
    th = list(theta); out = (_ct.c_double * 5)()
    _lib.arceq_readings(_d1(alpha), _d1(beta), _i(len(th)),
                        _vec(th) if th else (_ct.c_double * 1)(), _i(ngrid), out)
    return list(out)

_lib.arceq_readings_vsl.argtypes = [_dp, _dp, _ip, _dp]
def eq_readings_vsl(lam, delta, ngrid=20001):
    """Quantile-domain induction readings and equivalence discrepancy for tilted beta-kernel
    quantile densities q(u) = u^alpha (1-u)^beta exp(sum theta_j P_j(u))"""
    out = (_ct.c_double * 5)()
    _lib.arceq_readings_vsl(_d1(lam), _d1(delta), _i(ngrid), out)
    return list(out)


def bc_q(u, alpha, beta):
    """Beta-companion quantile function (standardised support)."""
    n = len(u); out = (_ct.c_double * n)()
    _lib.arceq2_bc_q(_i(n), _d(u), _d1(alpha), _d1(beta), out)
    return list(out)

def bc_cdf(x, alpha, beta):
    """Beta-companion distribution function."""
    n = len(x); out = (_ct.c_double * n)()
    _lib.arceq2_bc_cdf(_i(n), _d(x), _d1(alpha), _d1(beta), out)
    return list(out)

def bc_pdf(x, alpha, beta):
    """Beta-companion density."""
    n = len(x); out = (_ct.c_double * n)()
    _lib.arceq2_bc_pdf(_i(n), _d(x), _d1(alpha), _d1(beta), out)
    return list(out)

def eq_ub_quad(alpha, beta):
    """Quadratic shoulder and mode of the kernel family: (u_b, u_c)."""
    out = (_ct.c_double * 2)()
    _lib.arceq2_ub_quad(_d1(alpha), _d1(beta), out)
    return [out[0], out[1]]

def eq_E(alpha, beta):
    """Closed-form equivalence discrepancy E(alpha, beta)."""
    out = _d1(0.0)
    _lib.arceq2_E(_d1(alpha), _d1(beta), out)
    return out[0]

def eq_bstar(alpha):
    """Equivalence curve beta*(alpha)."""
    out = _d1(0.0)
    _lib.arceq2_bstar(_d1(alpha), out)
    return out[0]
