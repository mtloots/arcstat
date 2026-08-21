"""arcstat: Python (ctypes) binding to the shared arc-length C back-end. Pure standard library.
The three bundled C sources (goodness of fit, distributions, Bayesian test) are compiled together on
first import; the same sources back the R package 'arcstat'."""
import ctypes as _ct
import warnings as _warnings, os as _os
_here = _os.path.dirname(_os.path.abspath(__file__))
_ext = "dylib" if _os.uname().sysname == "Darwin" else "so"
_libpath = _os.path.join(_here, "libarcstat." + _ext)
if not _os.path.exists(_libpath):                       # compile the three cores together on first import
    import subprocess as _sp
    _sp.check_call(["cc", "-O2", "-fPIC", "-shared", "-I", _here,
                    _os.path.join(_here, "arclen.c"), _os.path.join(_here, "arcdistc.c"),
                    _os.path.join(_here, "bayesarc.c"), _os.path.join(_here, "cfarc.c"), _os.path.join(_here, "arccirc.c"),
                    _os.path.join(_here, "arck4.c"), _os.path.join(_here, "arceq.c"),
                    _os.path.join(_here, "arceq2.c"), _os.path.join(_here, "arck4fit.c"), _os.path.join(_here, "arcvc.c"),
                    _os.path.join(_here, "arcmv.c"),
                    _os.path.join(_here, "arceqfit.c"),
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
_lib.arcc_gof.argtypes = [_dp, _ip, _dp, _dp, _ip, _ip, _dp]
def arcc_gof(theta, c3, mu, B=999, seed=4207):
    """Goodness of fit of a fitted circular arc-length member: arc-length spacings statistic on the
    probability-integral transform, whose reference law is distribution-free. Same bytes as R."""
    out = (_ct.c_double * 2)()
    _lib.arcc_gof(_vec(theta), _i(len(theta)), _d1(c3), _d1(mu), _i(B), _i(seed), out)
    return {"stat": out[0], "p": out[1]}

_lib.arcc_exact_ci.argtypes = [_dp, _ip, _ip, _ip, _dp, _ip, _ip, _ip, _dp, _dp]
def arcc_exact_ci(theta, cgrid, B=999, stat=0, group=0, level=0.10, seed=4207):
    """Exact confidence interval for |c3| of the wrapped order-three family, by inverting a Monte
    Carlo test on a rotation-invariant statistic. Level exact for any n and any B; the rotation is
    removed by invariance. Same bytes as the R front."""
    ng = len(cgrid)
    out = (_ct.c_double * (3 + ng))()
    _lib.arcc_exact_ci(_vec(theta), _i(len(theta)), _i(B), _i(seed), _vec(cgrid), _i(ng),
                       _i(stat), _i(group), _d1(level), out)
    o = list(out)
    return {"lower": o[0], "upper": o[1], "stat": o[2], "pcurve": o[3:3+ng]}

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
def _i2(seq):
    a = (_ct.c_int * len(seq))()
    a[:] = [int(v) for v in seq]
    return a

def _i(v): a = (_ct.c_int * 1)(); a[0] = v; return a
def _iv(v): a = (_ct.c_int * len(v))(); a[:] = [int(t) for t in v]; return a
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

def arcc_factorise(rho, margin=0.0, nodes=4096):
    """Fejer-Riesz spectral factorisation. The R front calls this factorise_rho().

    A factorisation exists only where the implied quantile density is non-negative, which
    arcc_qmin_rho measures; for a single moment with rho_0 = 1 that reduces to |rho| <= 1/2.
    Asked for an inadmissible spectrum this used to run the root finder anyway and return the
    coefficients of a factorisation that does not exist, silently and at any modulus. It now does
    what the fitting path in the C back end has always done: shrink towards the circular uniform
    until the spectrum is admissible, and REPORT that it did.

    Returns {"p": coefficients, "shrink": factor applied (1.0 if none), "admissible": whether the
    spectrum AS SUPPLIED was admissible}, mirroring the R front's list.
    """
    admissible = arcc_qmin_rho(rho, nodes) >= 0.0
    shrink = 1.0
    if not admissible:
        rho, shrink = arcc_admiss(rho, margin, nodes)
        _warnings.warn(
            "spectrum is not admissible: no Fejer-Riesz factorisation exists for it, so it was "
            "shrunk towards the circular uniform by a factor of %.6g. The returned coefficients "
            "describe the SHRUNKEN spectrum, not the one supplied." % shrink,
            RuntimeWarning, stacklevel=2)
    nm = len(rho); out = (_ct.c_double * (2*(nm+1)))()
    _lib.arcc_factorise(_d([r.real for r in rho]), _d([r.imag for r in rho]), _i(nm), out)
    return {"p": [complex(out[2*j], out[2*j+1]) for j in range(nm+1)],
            "shrink": shrink, "admissible": admissible}

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
_lib.arck4_band_model.argtypes = [_dp, _dp, _ip, _ip, _dp, _dp]
_lib.arck4_band_sample.argtypes= [_dp, _dp, _ip, _dp, _ip, _dp, _dp]
_lib.arck4_readings.argtypes   = [_dp, _ip, _dp]
_lib.arck4_fit_lmom.argtypes   = [_dp, _dp, _ip, _dp]
_lib.arck4_fit_aleq.argtypes   = [_dp, _ip, _dp, _ip, _dp, _dp]
_lib.arck4_fit_nls.argtypes    = [_dp, _dp, _ip, _dp, _dp]
_lib.arck4_fit_nalr.argtypes   = [_dp, _dp, _ip, _ip, _dp, _ip, _dp, _dp, _dp]

def _pcode(p):
    """The tropical norm is p = inf; the C signals it with any non-positive p."""
    import math
    return -1.0 if (isinstance(p,(int,float)) and math.isinf(p) and p > 0) else float(p)

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

def k4_band_model(theta, breaks, nodes=60, p=2.0):
    """Band arc lengths of a scaled kappa curve, under the p-norm (p=inf is the tropical one)"""
    J = len(breaks) - 1; out = (_ct.c_double * J)()
    _lib.arck4_band_model(_vec(theta), _vec(breaks), _i(J), _i(nodes), _d1(_pcode(p)), out)
    return list(out)

def k4_band_sample(x, y, breaks, p=2.0):
    """Band arc lengths of a data polyline, under the p-norm (p=inf is the tropical one)"""
    J = len(breaks) - 1; out = (_ct.c_double * J)()
    _lib.arck4_band_sample(_vec(x), _vec(y), _i(len(x)), _vec(breaks), _i(J), _d1(_pcode(p)), out)
    return list(out)

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

def k4_fit_nalr(x, y, start, J=12, lam=1.0, w=9, p=2.0):
    """Deterministic kappa fits: L-moment shape inversion, quantile-domain arc-length shape fit,
    curve-domain NLS and NALR (banded arc lengths of the running-median presmoothed polyline)"""
    out = (_ct.c_double * 7)()
    _lib.arck4_fit_nalr(_vec(x), _vec(y), _i(len(x)), _i(J), _d1(lam), _i(w),
                        _d1(_pcode(p)), _vec(start), out)
    return list(out)

_lib.arceqfit_bc.argtypes = [_dp, _dp, _ip, _dp, _ip, _ip, _ip, _dp, _dp, _ip, _dp, _dp, _dp]
def eqfit_bc(x, y, starts, curve, maxit=4000):
    """Beta-companion sigmoid fit (free: 6-column starts; equivalence-constrained: 5-column),
    deterministic Nelder-Mead in the shared C, two passes per start, best kept.
    Returns (theta[6], sse, a, b). Same semantics and bytes as the R front."""
    flat = [float(v) for row in starts for v in row]
    np_ = len(starts[0]); ns = len(starts)
    cal = [float(r[0]) for r in curve]; cbe = [float(r[1]) for r in curve]
    out = (_ct.c_double * 9)()
    _lib.arceqfit_bc(_vec(x), _vec(y), _i(len(x)), _vec(flat), _i(ns), _i(np_), _i(maxit),
                     _vec(cal), _vec(cbe), _i(len(cal)), _d1(min(cal)), _d1(max(cal)), out)
    return list(out[0:6]), out[6], out[7], out[8]

_lib.arceqfit_estsim.argtypes = [_dp, _dp, _ip, _ip, _ip, _ip, _ip, _dp]
def eqfit_estsim(a0, b0, R, nsizes, seed=4207, maxit=500):
    """Three-estimator study of the equivalence member: MLE, moments, L-moments, OpenMP over
    replicates with per-replicate splitmix64 streams. Returns the flat R x 3 x 2 x len(nsizes)
    array in R's column-major order."""
    nn = len(nsizes)
    out = (_ct.c_double * (R * 3 * 2 * nn))()
    _lib.arceqfit_estsim(_d1(a0), _d1(b0), _i(R), _iv(nsizes), _i(nn), _i(seed), _i(maxit), out)
    return list(out)

_lib.arceqfit_msim.argtypes = [_ip, _dp, _ip, _ip, _dp, _dp, _ip, _dp, _dp, _ip, _dp, _dp, _dp, _dp]
def eqfit_msim(nM, sde, Rm, th0, st5, curve, awin, thref, seed=4207, maxit=2000):
    """The manifold-test simulation of the equivalence paper: four arms (truth, three
    displacements), two-stage profiled minimisation over the manifold interpolant inside the
    identifiable nuisance set, truth-arm audit and displaced-arm reference certificates.
    thref is a 3x4 sequence of (mu, sigma, alpha, beta) reference members. Same bytes as R."""
    cal = [float(r[0]) for r in curve]; cbe = [float(r[1]) for r in curve]
    tr = [float(v) for row in thref for v in row]
    out = (_ct.c_double * (8 * Rm))()
    _lib.arceqfit_msim(_i(nM), _d1(sde), _i(Rm), _i(seed), _vec(th0), _vec(st5), _i(maxit),
                       _vec(cal), _vec(cbe), _i(len(cal)), _d1(awin[0]), _d1(awin[1]),
                       _vec(tr), out)
    o = list(out)
    return {k: o[j*Rm:(j+1)*Rm] for j, k in enumerate(
        ["Tlev","Taud","Tp2","Tp15","Tp4","Tref2","Tref15","Tref4"])}

_lib.arceqfit_nullT.argtypes = [_ip, _dp, _ip, _ip, _dp, _dp, _ip, _dp, _dp, _ip, _dp, _dp, _dp]
def eqfit_nullT(nM, sde, R2, th0, st5, curve, awin, seed=4207, maxit=2000):
    """Null draws of the profiled manifold statistic at one member; same bytes as R."""
    cal = [float(r[0]) for r in curve]; cbe = [float(r[1]) for r in curve]
    out = (_ct.c_double * R2)()
    _lib.arceqfit_nullT(_i(nM), _d1(sde), _i(R2), _i(seed), _vec(th0), _vec(st5), _i(maxit),
                        _vec(cal), _vec(cbe), _i(len(cal)), _d1(awin[0]), _d1(awin[1]), out)
    return list(out)

_lib.arceqfit_blocklen.argtypes = [_dp, _ip, _ip]
def eqfit_blocklen(res):
    """Moving-block length rule of the residual bootstraps: first lag at which the sample
    autocorrelation enters the 2/sqrt(n) noise band, with an n^(1/3) fallback, floor 2."""
    out = (_ct.c_int * 1)()
    _lib.arceqfit_blocklen(_vec(res), _i(len(res)), out)
    return out[0]

_lib.arceqfit_bsweep.argtypes = [_dp, _ip, _dp, _dp]
def k4_b_sweep(k, h):
    """Derivative reading of the standardised kappa member over a k grid at fixed h."""
    out = (_ct.c_double * len(k))()
    _lib.arceqfit_bsweep(_vec(k), _i(len(k)), _d1(h), out)
    return list(out)

_lib.arceqfit_Esweep.argtypes = [_dp, _dp, _ip, _dp]
def eq_E_sweep(alpha, betas):
    """Closed-form equivalence discrepancy E over a beta grid at fixed alpha."""
    out = (_ct.c_double * len(betas))()
    _lib.arceqfit_Esweep(_d1(alpha), _vec(betas), _i(len(betas)), out)
    return list(out)

_lib.arck4_mv_boot.argtypes = [_dp, _dp, _ip, _dp, _ip, _ip, _dp, _ip, _dp, _dp]
def k4_mv_boot(x, y, cut=float("-inf"), B=50, seed=4207, bounds=(-0.98, 0.95, 0.0, 4.0),
               maxit=1500):
    """The fitting-method multiverse of one induction run: eight admissible pipelines under one
    moving-block residual bootstrap, all in the shared C back end. Same pipelines, replicate
    policy and per-replicate splitmix64 streams as the R front; byte-identical output whatever
    the thread count. Returns (a_pt, A): the eight point estimates in method order (LS, med9-LS,
    med21-LS, arc, GEV-LS, transient-excised, grid-odd, grid-even), and B rows of replicate
    readings (None when B == 0)."""
    n = len(x)
    a_pt = (_ct.c_double * 8)()
    A = (_ct.c_double * (max(1, B) * 8))()
    _lib.arck4_mv_boot(_vec(x), _vec(y), _i(n),
                       _d1(cut if cut == cut and cut > float("-inf") else -1e300),
                       _i(B), _i(seed), _vec(bounds), _i(maxit), a_pt, A)
    rows = [list(A[r * 8:(r + 1) * 8]) for r in range(B)] if B > 0 else None
    return list(a_pt), rows

_lib.arck4_fit_varpro.argtypes = [_dp, _dp, _ip, _dp, _ip, _ip, _dp, _dp, _dp]
def k4_fit_varpro(x, y, starts, maxit=1500, bounds=(-0.98, 0.95, 0.0, 4.0), hfix=None):
    """Variable-projection fit of the drifted kappa response y = g0 + m x + g1 F(x; mu, sg, k, h).
    The three coefficients are linear and solve exactly at any shape, so only the four shape
    parameters are searched, from a grid of starts; the C back end runs the starts in parallel and
    picks the best serially, so the answer does not depend on the thread count. starts is a sequence
    of 4-tuples (mu, log sigma, k, log h); bounds is the admissible shape box
    (k_lo, k_hi, h_lo, h_hi), defaulting to the one the edible-oil analysis settled on.
    hfix holds h at a given value and fits the other three, which is how the submodels are
    fitted in their own right: h = 0 is the generalised extreme value distribution and, with
    k = 0, the Gumbel. Returns (g0, m, g1, mu, sigma, k, h, rss)."""
    flat = [float(v) for row in starts for v in row]
    ns = len(flat) // 4
    if ns * 4 != len(flat):
        raise ValueError("starts must have four columns: (mu, log sigma, k, log h)")
    out = (_ct.c_double * 8)()
    _lib.arck4_fit_varpro(_vec(x), _vec(y), _i(len(x)), _vec(flat), _i(ns), _i(maxit), _vec(bounds),
                          _vec([-1.0 if hfix is None else float(hfix)]), out)
    return list(out)

_lib.arck4_band_trop_mean.argtypes = [_dp, _ip, _dp, _dp, _dp, _ip, _dp]
def k4_band_trop_mean(x, theta, sigma, breaks):
    """Exact mean of the tropical band arc length under Gaussian error. The max-plus element
    max(dx,|dy|) has an expectation elementary in Phi and phi, where the Euclidean element's is a
    confluent hypergeometric; the observed band arc length can therefore be compared with its own
    mean. Only the mean: consecutive increments share an observation, so a band's variance is not
    the sum of its elements' variances and summing as if independent understates it by about a
    fifth."""
    J = len(breaks) - 1
    if len(theta) != 6 or J < 1:
        raise ValueError("theta must have six entries and breaks at least two")
    out = (_ct.c_double * J)()
    _lib.arck4_band_trop_mean(_vec(x), _i(len(x)), _vec(theta), _d1(float(sigma)),
                              _vec(breaks), _i(J), out)
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

def icc_oneway(y, g):
    """One-way variance components: (icc, sd_within, sd_between).

    Groups may be unequal in size; the between-group mean square is divided by the unbalanced
    constant k0 = (N - sum n_i^2 / N) / (G - 1), not by the mean group size. A negative
    between-group variance estimate is truncated at zero.
    """
    lab = {}
    gi = []
    for v in g:
        if v not in lab:
            lab[v] = len(lab)
        gi.append(lab[v])
    n = len(y)
    out = (_ct.c_double * 3)()
    _lib.arcvc_icc_oneway(_d(y), _i2(gi), _i(n), _i(len(lab)), out)
    return [out[0], out[1], out[2]]

def locus_dist(h, k, locus_h, locus_k):
    """Shortest distance in the (h,k) plane from each shape to the locus polyline.

    Measured to the polyline's segments rather than its vertices.
    """
    n = len(h)
    out = (_ct.c_double * n)()
    _lib.arcvc_locus_dist(_d(h), _d(k), _i(n), _d(locus_h), _d(locus_k), _i(len(locus_h)), out)
    return list(out)
