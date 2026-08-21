"""ellstat: Python (ctypes) binding to the shared arc-length C back-end. Pure standard library.

The same C sources back the R package 'ellstat'; parity harnesses in the repository require every
exported quantity to agree byte for byte between the two front ends.

LOADING, and why it has two paths. An INSTALLED package carries a library that pip built at
install time from setup.py -- no compiler is needed here and nothing is written to disk. A
SOURCE CHECKOUT has the .c files but no built library, so it is compiled on demand into a cache
directory to keep the edit-a-C-file-and-re-import development loop working. The previous version
compiled on every fresh import AND wrote into its own site-packages directory, which fails on
read-only installs and needs a compiler at import rather than install.
"""
import ctypes as _ct
import glob as _glob
import os as _os
import sys as _sys
import warnings as _warnings

_here = _os.path.dirname(_os.path.abspath(__file__))

def _find_built():
    """A library placed beside this file by pip at install time."""
    pats = ("_libellstat*.so", "_libellstat*.dylib", "_libellstat*.pyd", "_libellstat*.dll")
    for pat in pats:
        hits = sorted(_glob.glob(_os.path.join(_here, pat)))
        if hits:
            return hits[0]
    return None

def _cache_dir():
    """Somewhere writable, preferring the package directory only if it actually is."""
    if _os.access(_here, _os.W_OK):
        return _here
    base = _os.environ.get("XDG_CACHE_HOME") or _os.path.join(_os.path.expanduser("~"), ".cache")
    d = _os.path.join(base, "ellstat")
    _os.makedirs(d, exist_ok=True)
    return d

def _compile_from_source():
    """Development fallback: build the bundled sources. Never reached in an installed package."""
    import subprocess as _sp
    srcs = sorted(f for f in _glob.glob(_os.path.join(_here, "*.c"))
                  if _os.path.basename(f) != "init.c")
    if not srcs:
        return None
    ext = "dll" if _sys.platform == "win32" else ("dylib" if _sys.platform == "darwin" else "so")
    out = _os.path.join(_cache_dir(), "libellstat." + ext)
    cc = _os.environ.get("CC") or "cc"
    try:
        _sp.check_call([cc, "-O2", "-fPIC", "-shared", "-I", _here] + srcs + ["-lm", "-o", out])
    except (OSError, _sp.CalledProcessError) as exc:
        raise ImportError(
            "ellstat could not build its C back-end from the bundled sources (%s). "
            "In an installed copy this should never happen -- pip builds the library at install "
            "time. In a source checkout, a C compiler is required; set CC if it is not 'cc'." % exc
        )
    return out

_libpath = _find_built()
if _libpath is None:
    _cached = _os.path.join(_cache_dir(), "libellstat." +
                            ("dll" if _sys.platform == "win32" else
                             "dylib" if _sys.platform == "darwin" else "so"))
    _libpath = _cached if _os.path.exists(_cached) else _compile_from_source()
if _libpath is None:
    raise ImportError("ellstat: no compiled back-end and no C sources to build one from.")

_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)

for _nm, _at in {
    "C_ell_ke":         [_dp, _ip, _dp, _dp],
    "C_ell_nome":       [_dp, _ip, _dp],
    "C_ell_jac":        [_dp, _ip, _dp, _dp, _dp, _dp],
    "C_ell_basis":      [_dp, _ip, _ip, _dp, _dp, _dp],
    "C_ell_smom":       [_dp, _ip, _ip, _dp, _dp],
    "C_ell_transfer":   [_dp, _ip, _ip, _dp],
    "C_ell_dn2mom":     [_dp, _ip, _dp],
    "C_ell_evm_d":      [_dp, _ip, _dp, _dp, _dp, _ip, _dp],
    "C_ell_evm_r":      [_ip, _dp, _dp, _dp, _dp, _dp],
    "C_ell_moment_var": [_dp, _ip, _ip, _dp, _dp],
    "C_ell_evm_pop":    [_dp, _dp, _dp, _ip, _ip, _dp],
    "C_ell_evm_fit":    [_dp, _ip, _dp, _ip, _ip, _dp],
    "C_ell_exact_unif": [_dp, _ip, _dp, _ip, _dp, _dp],
    "C_ell_exact_ci":   [_dp, _ip, _dp, _dp, _ip, _ip, _dp, _dp, _ip, _ip, _dp],
    "C_ell_effstudy":   [_dp, _dp, _ip, _ip, _ip, _dp, _dp],
    "C_ell_pj_dens":    [_dp, _ip, _dp, _dp, _dp, _ip, _dp],
    "C_ell_pj_rand":    [_ip, _dp, _dp, _dp, _dp, _dp],
    "C_ell_pj_effcurve":[_dp, _dp, _dp, _dp, _ip, _ip, _ip, _dp, _dp],
    "C_ell_pj_kopt":    [_dp, _dp, _dp, _ip, _ip, _dp, _dp],
    "C_ell_projfam_dens":[_dp, _ip, _dp, _dp, _ip, _dp],
    "C_ell_projfam_rand":[_ip, _dp, _dp, _ip, _dp, _dp],
    "C_ell_projfam_eff": [_dp, _dp, _ip, _dp, _ip, _ip, _dp, _dp],
    "C_ell_projfam_kopt":[_dp, _dp, _ip, _ip, _dp, _dp],
    "C_ell_cat_eff":     [_ip, _dp, _dp, _ip, _dp, _ip, _ip, _dp, _dp],
    "C_ell_cat_kopt":    [_ip, _dp, _dp, _ip, _ip, _dp, _dp],
    "C_ell_joint_cov":   [_dp, _ip, _ip, _dp, _dp, _dp],
    "C_ell_joint_samp":  [_dp, _ip, _ip, _dp, _dp, _dp],
    "C_ell_joint_unif":  [_dp, _ip, _ip, _dp, _ip, _ip, _dp, _dp],
}.items():
    getattr(_lib, _nm).argtypes = _at

def _d(x):
    if hasattr(x, "__len__"):
        return (_ct.c_double * len(x))(*[float(v) for v in x])
    return _ct.c_double(float(x))
def _i(x): return _ct.c_int(int(x))
def _o(n): return (_ct.c_double * n)()

def ell_KE(k):
    """Complete elliptic integrals K and E of modulus k."""
    ks = _d([k] if not hasattr(k, "__len__") else k)
    n = len(ks); K, E = _o(n), _o(n)
    _lib.C_ell_ke(ks, _ct.byref(_i(n)), K, E)
    return (list(K), list(E)) if hasattr(k, "__len__") else (K[0], E[0])

def ell_nome(k):
    """Nome q = exp(-pi K'/K)."""
    ks = _d([k] if not hasattr(k, "__len__") else k)
    n = len(ks); out = _o(n)
    _lib.C_ell_nome(ks, _ct.byref(_i(n)), out)
    return list(out) if hasattr(k, "__len__") else out[0]

def ell_jacobi(u, k):
    """Jacobi sn, cn, dn at argument u and modulus k."""
    us = _d(u if hasattr(u, "__len__") else [u]); n = len(us)
    sn, cn, dn = _o(n), _o(n), _o(n)
    _lib.C_ell_jac(us, _ct.byref(_i(n)), _ct.byref(_d(k)), sn, cn, dn)
    return list(sn), list(cn), list(dn)

def ell_basis(t, P, k):
    """The elliptic moment basis c_p, s_p; returns two P-by-len(t) lists."""
    ts = _d(t); nt = len(ts); P = int(P)
    C, S = _o(P * nt), _o(P * nt)
    _lib.C_ell_basis(ts, _ct.byref(_i(nt)), _ct.byref(_i(P)), _ct.byref(_d(k)), C, S)
    return ([list(C[p*nt:(p+1)*nt]) for p in range(P)],
            [list(S[p*nt:(p+1)*nt]) for p in range(P)])

def ell_moments(t, P=2, k=0.9):
    """Sample elliptic moments; returns a list of (c, s) pairs."""
    ts = _d(t); P = int(P); out = _o(2 * P)
    _lib.C_ell_smom(ts, _ct.byref(_i(len(ts))), _ct.byref(_i(P)), _ct.byref(_d(k)), out)
    return [(out[2*p], out[2*p+1]) for p in range(P)]

def ell_moment_var(t, P=3, k=0.9):
    """Exact sample theory: (|hat E_p|, realised trace variance, (1-|E|^2)/n)."""
    ts = _d(t); P = int(P); out = _o(3 * P)
    _lib.C_ell_moment_var(ts, _ct.byref(_i(len(ts))), _ct.byref(_i(P)),
                          _ct.byref(_d(k)), out)
    return [(out[3*p], out[3*p+1], out[3*p+2]) for p in range(P)]

def ell_transfer(k, P=3, M=5):
    """Transfer coefficients from trigonometric to elliptic moments."""
    P, M = int(P), int(M); out = _o(P * M)
    _lib.C_ell_transfer(_ct.byref(_d(k)), _ct.byref(_i(P)), _ct.byref(_i(M)), out)
    return [list(out[p*M:(p+1)*M]) for p in range(P)]

def ell_dn2_moments(k, P=4):
    """Fourier coefficients of the arc-length family on an ellipse."""
    P = int(P); out = _o(P + 1)
    _lib.C_ell_dn2mom(_ct.byref(_d(k)), _ct.byref(_i(P)), out)
    return list(out)

def dellvm(t, mu=0.0, kappa=1.0, k=0.9, ng=1024):
    """Elliptic von Mises density."""
    ts = _d(t); out = _o(len(ts))
    _lib.C_ell_evm_d(ts, _ct.byref(_i(len(ts))), _ct.byref(_d(mu)),
                     _ct.byref(_d(kappa)), _ct.byref(_d(k)), _ct.byref(_i(ng)), out)
    return list(out)

def rellvm(n, mu=0.0, kappa=1.0, k=0.9, seed=1):
    """Simulate from the elliptic von Mises."""
    n = int(n); out = _o(n)
    _lib.C_ell_evm_r(_ct.byref(_i(n)), _ct.byref(_d(mu)), _ct.byref(_d(kappa)),
                     _ct.byref(_d(k)), _ct.byref(_d(seed)), out)
    return list(out)

def ell_evm_moments(mu, kappa, k, P=2, ng=1024):
    """Population elliptic moments of the elliptic von Mises."""
    P = int(P); out = _o(2 * P)
    _lib.C_ell_evm_pop(_ct.byref(_d(mu)), _ct.byref(_d(kappa)), _ct.byref(_d(k)),
                       _ct.byref(_i(P)), _ct.byref(_i(ng)), out)
    return [(out[2*p], out[2*p+1]) for p in range(P)]

def fit_ellvm(t, k=0.9, P=2, ng=512):
    """Minimum-distance fit of the elliptic von Mises; returns (mu, kappa)."""
    ts = _d(t); out = _o(2)
    _lib.C_ell_evm_fit(ts, _ct.byref(_i(len(ts))), _ct.byref(_d(k)),
                       _ct.byref(_i(P)), _ct.byref(_i(ng)), out)
    return out[0], out[1]

def ell_exact_uniform(t, k=0.9, B=999, seed=1):
    """Exact test of circular uniformity by rotation invariance."""
    ts = _d(t); out = _o(2)
    _lib.C_ell_exact_unif(ts, _ct.byref(_i(len(ts))), _ct.byref(_d(k)),
                          _ct.byref(_i(B)), _ct.byref(_d(seed)), out)
    return {"statistic": out[0], "p_value": out[1]}

def ell_exact_ci(t, k=0.9, kgrid=None, B=499, level=0.05, seed=1, ng=512, nmu=8):
    """Exact confidence interval for the concentration by test inversion.

    *** ARGUMENT LIST FIXED 21 Aug 2026. *** C_ell_exact_ci takes ELEVEN parameters; this front
    passed ten, omitting nmu. The C therefore read the output pointer as `int *nmu` and wrote its
    results through whatever followed -- a wild write, so this routine had never worked from
    Python and segfaulted the interpreter. The R front passes nmu = 8L and was correct. The
    allocation, 3 + ngr, was right all along; only the argument list was short.
    """
    if kgrid is None:
        kgrid = [0.05 + (10.0 - 0.05) * i / 59.0 for i in range(60)]
    ts = _d(t); kg = _d(kgrid); ngr = len(kgrid); out = _o(3 + ngr)
    _lib.C_ell_exact_ci(ts, _ct.byref(_i(len(ts))), _ct.byref(_d(k)), kg,
                        _ct.byref(_i(ngr)), _ct.byref(_i(B)), _ct.byref(_d(level)),
                        _ct.byref(_d(seed)), _ct.byref(_i(ng)), _ct.byref(_i(nmu)), out)
    return {"lower": out[0], "upper": out[1], "statistic": out[2],
            "curve": list(zip(list(kgrid), list(out)[3:]))}

def ell_efficiency(k, kappa=3.0, n=250, R=2000, ng=1024, seed=1):
    """Matched-modulus efficiency.

    Returns (eta, var_trig, var_elliptic, gain, solved_trig, solved_elliptic), matching the R
    front's named vector.

    *** BUFFER SIZE FIXED 21 Aug 2026. *** C_ell_effstudy writes out[0] through out[4], five
    doubles, and this allocated THREE -- a sixteen-byte heap overflow on every call. The R front
    allocates double(5) and was correct. It did not crash at the call site: it corrupted the heap
    and brought the interpreter down later, which is why it survived until a test suite made
    several calls in one process. Any change to what the C writes must be matched here.
    """
    out = _o(5)
    _lib.C_ell_effstudy(_ct.byref(_d(k)), _ct.byref(_d(kappa)), _ct.byref(_i(n)),
                        _ct.byref(_i(R)), _ct.byref(_i(ng)), _ct.byref(_d(seed)), out)
    return out[0], out[1], out[2], out[1] / out[2], out[3], out[4]

def dprojell(t, tau=0.3, d=0.6, beta=2.0, ng=400):
    """Density of the projected family (offset ellipse, generalised gamma scale)."""
    ts = _d(t); out = _o(len(ts))
    _lib.C_ell_pj_dens(ts, _ct.byref(_i(len(ts))), _ct.byref(_d(tau)),
                       _ct.byref(_d(d)), _ct.byref(_d(beta)), _ct.byref(_i(ng)), out)
    return list(out)

def rprojell(n, tau=0.3, d=0.6, beta=2.0, seed=1):
    """Simulate from the projected family."""
    n = int(n); out = _o(n)
    _lib.C_ell_pj_rand(_ct.byref(_i(n)), _ct.byref(_d(tau)), _ct.byref(_d(d)),
                       _ct.byref(_d(beta)), _ct.byref(_d(seed)), out)
    return list(out)

def ell_eff_curve(tau=0.3, d=0.6, beta=2.0, kgrid=None, nt=4096, ng=400, eps=1e-4):
    """Asymptotic efficiency corr^2(c_1(.,k), score) of the moment estimator."""
    if kgrid is None:
        kgrid = [0.0, 0.9, 0.99, 0.999, 0.99999]
    kg = _d(kgrid); ngr = len(kgrid); out = _o(ngr + 2)
    _lib.C_ell_pj_effcurve(_ct.byref(_d(tau)), _ct.byref(_d(d)), _ct.byref(_d(beta)),
                           kg, _ct.byref(_i(ngr)), _ct.byref(_i(nt)),
                           _ct.byref(_i(ng)), _ct.byref(_d(eps)), out)
    return {"efficiency": list(out)[:ngr], "fisher": out[ngr],
            "godambe": out[ngr + 1], "kgrid": list(kgrid)}

def ell_kopt(tau=0.3, d=0.6, beta=2.0, nt=4096, ng=400, eps=1e-4):
    """The optimal modulus and its efficiency. The optimum may be interior."""
    out = _o(2)
    _lib.C_ell_pj_kopt(_ct.byref(_d(tau)), _ct.byref(_d(d)), _ct.byref(_d(beta)),
                       _ct.byref(_i(nt)), _ct.byref(_i(ng)), _ct.byref(_d(eps)), out)
    return out[0], out[1]

_FAMS = ("vonmises", "wrappednormal", "wrappedcauchy", "cardioid",
         "projnormal", "projcauchy", "ellvonmises", "arcellipse")

def dprojfam(t, tau=1.0, d=1.0, family="normal"):
    """Density of the projected normal or projected Cauchy."""
    w = 0 if family == "normal" else 1
    ts = _d(t); out = _o(len(ts))
    _lib.C_ell_projfam_dens(ts, _ct.byref(_i(len(ts))), _ct.byref(_d(tau)),
                            _ct.byref(_d(d)), _ct.byref(_i(w)), out)
    return list(out)

def rprojfam(n, tau=1.0, d=1.0, family="normal", seed=1):
    """Simulate from the projected normal or projected Cauchy."""
    w = 0 if family == "normal" else 1
    n = int(n); out = _o(n)
    _lib.C_ell_projfam_rand(_ct.byref(_i(n)), _ct.byref(_d(tau)), _ct.byref(_d(d)),
                            _ct.byref(_i(w)), _ct.byref(_d(seed)), out)
    return list(out)

def ell_projfam_eff(tau=1.0, d=1.0, family="normal", kgrid=None, nt=4096, eps=1e-4):
    """Efficiency curve for a cited projected family."""
    if kgrid is None:
        kgrid = [0.0, 0.9, 0.99, 0.999, 0.99999]
    w = 0 if family == "normal" else 1
    kg = _d(kgrid); ngr = len(kgrid); out = _o(ngr + 2)
    _lib.C_ell_projfam_eff(_ct.byref(_d(tau)), _ct.byref(_d(d)), _ct.byref(_i(w)),
                           kg, _ct.byref(_i(ngr)), _ct.byref(_i(nt)),
                           _ct.byref(_d(eps)), out)
    return {"efficiency": list(out)[:ngr], "fisher": out[ngr], "mass": out[ngr + 1]}

def ell_projfam_kopt(tau=1.0, d=1.0, family="normal", nt=4096, eps=1e-4):
    """The modulus that maximises efficiency for a cited projected family, and that efficiency.

    The R front has exported this since the package was written; the Python binding was simply
    never added, even though the C entry point is compiled into this library and its argtypes
    were already declared above. Returns (k, efficiency), matching R's named vector.
    """
    w = 0 if family == "normal" else 1
    out = _o(2)
    _lib.C_ell_projfam_kopt(_ct.byref(_d(tau)), _ct.byref(_d(d)), _ct.byref(_i(w)),
                            _ct.byref(_i(nt)), _ct.byref(_d(eps)), out)
    return out[0], out[1]

def ell_catalogue(family, p1, p2=1.0, order=1, kgrid=None, nt=4096, eps=1e-4):
    """The catalogue of circular distributions, after Hosking's table."""
    if kgrid is None:
        kgrid = [0.0, 0.9, 0.99, 0.999, 0.99999]
    f = _FAMS.index(family)
    kg = _d(kgrid); ngr = len(kgrid); out = _o(ngr + 2); ko = _o(2)
    _lib.C_ell_cat_eff(_ct.byref(_i(f)), _ct.byref(_d(p1)), _ct.byref(_d(p2)),
                       _ct.byref(_i(order)), kg, _ct.byref(_i(ngr)),
                       _ct.byref(_i(nt)), _ct.byref(_d(eps)), out)
    _lib.C_ell_cat_kopt(_ct.byref(_i(f)), _ct.byref(_d(p1)), _ct.byref(_d(p2)),
                        _ct.byref(_i(order)), _ct.byref(_i(nt)),
                        _ct.byref(_d(eps)), ko)
    return {"efficiency": list(out)[:ngr], "fisher": out[ngr], "mass": out[ngr + 1],
            "kopt": ko[0], "eff_kopt": ko[1]}

def ell_joint(f, P=2, k=0.9):
    """Exact joint mean and covariance of the sample elliptic moments."""
    fs = _d(f); D = 2 * int(P)
    mean = _o(D); Sig = _o(D * D)
    _lib.C_ell_joint_cov(fs, _ct.byref(_i(len(fs))), _ct.byref(_i(int(P))),
                         _ct.byref(_d(k)), mean, Sig)
    return {"mean": list(mean),
            "Sigma": [list(Sig)[r * D:(r + 1) * D] for r in range(D)]}

def ell_joint_sample(t, P=2, k=0.9):
    """Realised mean and covariance of the basis on a sample."""
    ts = _d(t); D = 2 * int(P)
    mean = _o(D); Sig = _o(D * D)
    _lib.C_ell_joint_samp(ts, _ct.byref(_i(len(ts))), _ct.byref(_i(int(P))),
                          _ct.byref(_d(k)), mean, Sig)
    return {"mean": list(mean),
            "Sigma": [list(Sig)[r * D:(r + 1) * D] for r in range(D)]}

def ell_joint_uniform(t, P=2, k=0.9, B=999, ng=4096, seed=1):
    """Simultaneous exact test of uniformity from the joint moment vector."""
    ts = _d(t); out = _o(2)
    _lib.C_ell_joint_unif(ts, _ct.byref(_i(len(ts))), _ct.byref(_i(int(P))),
                          _ct.byref(_d(k)), _ct.byref(_i(B)), _ct.byref(_i(ng)),
                          _ct.byref(_d(seed)), out)
    return {"statistic": out[0], "p_value": out[1]}
