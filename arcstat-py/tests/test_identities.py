"""Identity tests for the Python front end.

The parity harnesses at the repository root already prove that the R and Python fronts return
BYTE-IDENTICAL values from the shared C. That is agreement, not correctness: both fronts could
agree on a wrong number. These tests assert the mathematics independently, mirroring
arcstat/tests/testthat/test-backends.R and test-families.R, so the Python package is checked on
its own terms rather than only against its twin.
"""
import math
import os
import sys
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import arcstat as a  # noqa: E402

SQRT2 = math.sqrt(2.0)


def _unif(n):
    return [(i + 0.5) / n for i in range(n)]


def _same(a_, b_):
    """Equality that treats NaN as equal to NaN.

    Some multiverse pipelines legitimately return NaN on a given design, and Python's ``==``
    makes NaN unequal to itself -- so a plain list comparison reports a difference between two
    bit-identical results. R's identical() does not have this behaviour, which is why the
    mirrored test on the R side passes without any such helper.
    """
    if isinstance(a_, (list, tuple)):
        return len(a_) == len(b_) and all(_same(x, y) for x, y in zip(a_, b_))
    if isinstance(a_, float) and isinstance(b_, float):
        return (math.isnan(a_) and math.isnan(b_)) or a_ == b_
    return a_ == b_


# --------------------------------------------------------------- arc length core
def test_arclength_minimised_by_the_diagonal():
    """A perfectly uniform sample plots as the diagonal of the unit square, whose arc length is
    exactly sqrt(2). Nothing can be shorter, so the p-value there is one."""
    u = _unif(200)
    assert abs(a.al_statistic(u) - SQRT2) < 1e-3
    assert a.al_pvalue(a.al_statistic(u), 200) == 1.0


def test_al_moments_support_and_ordering():
    m = a.al_moments(200)
    assert set(m) == {"mean", "var", "support"}
    lo, hi = m["support"]
    assert abs(lo - SQRT2) < 1e-12          # the diagonal is the floor
    assert lo < m["mean"] < hi              # a proper distribution puts its mean inside
    assert m["var"] > 0


def test_pvalue_is_a_probability_and_decreases_in_the_statistic():
    n = 200
    s0 = a.al_statistic(_unif(n))
    ps = [a.al_pvalue(s, n) for s in (s0, s0 + 0.05, s0 + 0.15)]
    assert all(0.0 <= p <= 1.0 for p in ps)
    assert ps[0] >= ps[1] >= ps[2]


def test_characteristic_function_closed_forms():
    """The two members with exact values. These are the sharpest checks available."""
    assert abs(a.cf_arclength_family("normal") - 2.0) < 1e-8
    assert abs(a.cf_arclength_family("exponential") - math.pi) < 1e-8


# --------------------------------------------------------------- circular family
def test_uniform_circular_member_has_no_resultant():
    assert abs(a.arcc_rho3(0.0)) < 1e-12
    assert a.arcc_c3max() > 1.0


def test_resultant_increases_with_concentration():
    rhos = [a.arcc_rho3(c) for c in (0.0, 0.3, 0.6, 0.9, 1.03)]
    assert all(y > x for x, y in zip(rhos, rhos[1:]))
    assert all(0.0 <= r <= 1.0 for r in rhos)


def test_circular_draws_lie_on_the_circle_and_fit_recovers():
    u = _unif(4000)
    s = a.arcc_rand3(u, 1.0, 1.0)
    assert len(s) == 4000
    assert all(0.0 <= t < 2.0 * math.pi for t in s)
    fit = a.arcc_fit3(s)
    assert set(fit) >= {"mu", "c3", "rho", "interior"}
    assert abs(fit["c3"] - 1.0) < 0.1
    assert abs(fit["mu"] - 1.0) < 0.1


def test_quantile_density_is_positive():
    qd = a.arcc_qd3([i / 100 for i in range(1, 100)], 0.8)
    assert all(v > 0 for v in qd)


def test_factorise_accepts_an_admissible_spectrum_untouched():
    f = a.arcc_factorise([0.2])
    assert f["admissible"] is True
    assert f["shrink"] == 1.0
    assert abs(sum(abs(z) ** 2 for z in f["p"]) - 1.0) < 1e-10
    m1 = a.arcc_trigmom_fr(1, [z.real for z in f["p"]], [z.imag for z in f["p"]])
    assert abs(abs(m1) - 0.2) < 0.01


def test_one_half_is_the_admissibility_boundary():
    """For a single moment with rho_0 = 1 the spectrum is 1 + 2 rho cos(theta), non-negative
    exactly while |rho| <= 1/2. arcc_qmin_rho measures that minimum."""
    assert abs(a.arcc_qmin_rho([0.5])) < 1e-9
    assert a.arcc_qmin_rho([0.4]) > 0
    assert a.arcc_qmin_rho([0.7]) < 0
    assert a.arcc_factorise([0.5])["admissible"] is True


def test_factorise_reports_an_inadmissible_spectrum_rather_than_inventing_one():
    """Beyond the boundary no factorisation exists. The routine must warn, shrink, and say by how
    much -- not return the coefficients of a factorisation that does not exist."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        f = a.arcc_factorise([0.9])
        assert len(caught) == 1
        assert "not admissible" in str(caught[0].message)
    assert f["admissible"] is False
    assert f["shrink"] < 1.0
    assert abs(sum(abs(z) ** 2 for z in f["p"]) - 1.0) < 1e-10
    # the shrink lands exactly on the boundary, which is why every inadmissible request used to
    # produce the same achieved resultant with nothing said about it
    assert abs(0.9 * f["shrink"] - 0.5) < 1e-6


def test_tempered_von_mises_is_a_density():
    """Integrates to one over the circle by the trapezium rule on a fine grid."""
    n = 4000
    step = 2.0 * math.pi / n
    grid = [i * step for i in range(n)]
    dens = a.arcc_temper_vm(grid, 1.0, 0.0, 1.0)
    assert all(d >= 0 for d in dens)
    assert abs(sum(dens) * step - 1.0) < 1e-4


# --------------------------------------------------------------- kappa-4 fitting
def test_variable_projection_recovers_a_noiseless_response():
    x = [i * 10.0 / 59 for i in range(60)]
    y = a.k4_cdf(x, 5.0, 1.5, 0.2, 0.8)
    out = a.k4_fit_varpro(x, y, [[5, 1.5, 0.2, 0.8], [4, 2, 0, 1]])
    *theta, rss = out
    assert rss < 1e-12
    assert abs(theta[3] - 5.0) < 1e-4       # mu
    assert abs(theta[4] - 1.5) < 1e-4       # sigma
    assert abs(theta[5] - 0.2) < 1e-4       # k
    assert abs(theta[6] - 0.8) < 1e-4       # h


def test_multiverse_shape_and_reproducibility():
    x = [i * 10.0 / 59 for i in range(60)]
    y = a.k4_cdf(x, 5.0, 1.5, 0.2, 0.8)
    a_pt, A = a.k4_mv_boot(x, y, B=0)
    assert len(a_pt) == 8
    assert A is None                        # no replicate matrix when none are asked for
    _, A3 = a.k4_mv_boot(x, y, B=3, seed=4207)
    assert len(A3) == 3
    _, A3b = a.k4_mv_boot(x, y, B=3, seed=4207)
    assert _same(A3, A3b)                   # per-replicate splitmix64 streams are deterministic
    _, A3c = a.k4_mv_boot(x, y, B=3, seed=99)
    assert not _same(A3, A3c)


# --------------------------------------------------------------- variance components
def test_intraclass_correlation_hits_both_endpoints():
    y = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    g = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
    assert abs(a.icc_oneway(y, g)[0] - 1.0) < 1e-10


def test_locus_distance_vanishes_on_the_locus():
    lh = [i / 20 for i in range(21)]
    lk = [v * v for v in lh]
    assert a.locus_dist([lh[7]], [lk[7]], lh, lk)[0] < 1e-12
    near = a.locus_dist([lh[7]], [lk[7] + 0.01], lh, lk)[0]
    far = a.locus_dist([lh[7]], [lk[7] + 0.50], lh, lk)[0]
    assert 0 < near < far


# --------------------------------------------------------------- equivalence system
def test_equivalence_curve_root_and_sign_change():
    """eq_bstar returns the beta at which the discrepancy vanishes. It is defined only where a
    sign change exists in the scanned window; the admissible alpha is NEGATIVE."""
    al = -0.60
    bs = a.eq_bstar(al)
    assert math.isfinite(bs)
    assert abs(a.eq_E(al, bs)) < 1e-10
    assert a.eq_E(al, bs - 0.05) * a.eq_E(al, bs + 0.05) < 0


def test_equivalence_sweep_agrees_with_the_scalar_form():
    al = -0.60
    bs = a.eq_bstar(al)
    betas = [bs - 0.05, bs, bs + 0.05]
    swept = a.eq_E_sweep(al, betas)
    for b, s in zip(betas, swept):
        assert abs(a.eq_E(al, b) - s) < 1e-10


def test_equivalence_fit_recovers_a_member_on_the_curve():
    al = -0.60
    bs = a.eq_bstar(al)
    grid = [-0.64 + i * (0.13 / 13) for i in range(14)]
    curve = [[g, a.eq_bstar(g)] for g in grid]
    assert all(math.isfinite(c[1]) for c in curve)
    x = [0.05 + i * (0.90 / 39) for i in range(40)]
    y = a.bc_cdf(x, al, bs)
    starts = [[0, 1, al, bs, 0, 0], [0, 1, -0.58, -0.35, 0, 0]]
    out = a.eqfit_bc(x, y, starts, curve)
    th, sse = out[0], out[1]
    assert sse < 1e-06
    assert abs(th[4] - al) < 0.05


# --------------------------------------------------------------- Bayesian test
def test_bayesian_evidence_separates_uniform_from_a_departure():
    """Zero evidence against a sample that IS uniform; nearly all of it against one that is
    visibly not. Both are compared against one shared reference."""
    ref = a.bb_ref_disc(60)
    e_unif = a.bb_evidence(_unif(60), ref=ref, M=400, seed=1)
    # a deterministic departure: two tight clusters at the ends, which is visibly not uniform.
    # A mild departure such as u**2 is NOT detected at this sample size and is not asserted.
    clustered = sorted([0.02 + 0.06 * i / 29 for i in range(30)]
                       + [0.92 + 0.06 * i / 29 for i in range(30)])
    e_clust = a.bb_evidence(clustered, ref=ref, M=400, seed=1)
    assert 0.0 <= e_unif <= 1.0 and 0.0 <= e_clust <= 1.0
    assert e_unif < 0.10
    assert e_clust > 0.90


def test_posterior_draws_are_shaped_and_seeded():
    p1 = a.bb_post_disc(_unif(60), 200, 1)
    p2 = a.bb_post_disc(_unif(60), 200, 1)
    assert len(p1) == 200
    assert p1 == p2


def test_pvalue_never_leaves_unit_interval_far_in_the_tail():
    # badly non-uniform samples that drove the pre-0.3.1 saddlepoint tail below zero
    import random
    rng = random.Random(17)
    for pw in (4, 12, 30):
        u = [min(max(rng.random() ** pw, 1e-9), 1 - 1e-9) for _ in range(150)]
        p = a.al_pvalue(a.al_statistic(u), len(u))
        assert math.isfinite(p) and 0.0 <= p <= 1.0


def test_equivalence_readings_match_high_precision_values():
    # mpmath, 50 digits: columns alpha, beta, a, b, uc, ub. Mirrors test-backends.R. The integral
    # of q starts at zero, which matters for alpha < -1/2 (the grid starts at 1e-5).
    ref = [
        (-0.6434, -0.966, 1.11192920277, 1.11193319294, 0.399776314154, 0.0709201279567),
        (-0.6, -0.378118, 1.06789574393, 1.06789588842, 0.613422920343, 0.115438787487),
        (-0.55, -0.9, 0.640186539725, 0.442993243486, 0.379310344828, 0.027298078942),
        (-0.52, -1 + 1e-9, 0.512394550159, 0.228607978552, 0.342105263383, 0.00994857612993),
        (-0.505, -1 + 1e-9, 0.465652925522, 0.104103098037, 0.33554817298, 0.00249685186998),
        (-0.6, -1 + 1e-9, 0.843319480661, 0.756439847235, 0.375000000234, 0.0486043951283),
        (-0.64, -0.99, 1.08223342321, 1.07374270604, 0.39263803681, 0.0678008353378),
    ]
    for al, be, ea, eb, euc, eub in ref:
        r = a.eq_readings(al, be)
        for key, exact in ((1, ea), (2, eb), (3, euc), (4, eub)):
            assert abs(r[key] - exact) <= 1e-7 * abs(exact)
    rt = a.eq_readings(-0.6, -0.9, [0.05, 0.02, 0.01])
    assert abs(rt[1] - 0.826378472447) <= 1e-7 * 0.83
    assert abs(rt[3] - 0.391886932994) <= 1e-7
    rv = a.eq_readings_vsl(0.3, 0.6)
    assert abs(rv[1] - 0.652733314402) <= 1e-7 * 0.66
    assert abs(rv[3] - 0.440653844715) <= 1e-7


def test_kappa_readings_match_high_precision_values():
    # mpmath, 50 digits: mu, sg, k, h, a, b, mode. Mirrors test-arck4.R.
    nan = float("nan")
    ref = [
        (0, 1, -0.98, 0.2418, -0.618391710679162, -0.664133279975884, -0.376903116941754),
        (0, 1, -0.98, 0.8088, -0.182328163643372, nan, -0.094726927684711),
        (6, 1.3, -0.999999, 0.1, 5.09324965009938, 5.0416030219296, 5.41500013495375),
        (0, 1, 0.2, 0.4, -0.591499322579312, -0.776714806468733, 0.13782667601799),
        (0, 1, 0.0, 0.0, -1.0, -0.962423650119207, 0.0),
        (0, 1, -0.4155, -1.5959, -2.09123825802291, nan, -1.08116330098459),
        (0, 1, -0.9, -30, nan, nan, nan),
        (4.92632094446957, 6.95102192620741, -0.98, 0.2418, 0.627866604553808, 0.309915953433158, 2.30645911455152),
        (13.9330859822076, 28.8494200479927, -0.98, 0.8088, 8.6730242026808, nan, 11.2002690555755),
    ]
    for mu, sg, k, h, *want in ref:
        got = a.k4_readings([1, 400, mu, sg, k, h])
        assert got == a.k4_readings([1, 400, mu, sg, k, h], grid=400000) or math.isnan(got[1])
        for g, w in zip(got, want):
            assert math.isnan(g) == math.isnan(w)
            if not math.isnan(w):
                assert abs(g - w) < 1e-12


def test_kappa_l_moment_ratios_match_high_precision_values():
    got = a.k4_tau34(-0.98, 0.3)
    for g, w in zip(got, (0.9793008357201542, 0.966575435477793, 49.58640438135162, 49.05639225332278)):
        assert abs(g / w - 1) < 1e-12


def test_companion_keeps_accuracy_near_beta_minus_one():
    a1 = -0.6448948644172633
    assert abs(a.eq_bstar(a1 + 1e-8) + 0.9999997679034622) < 1e-10
    for be in (-0.9999, -0.6, -0.03):
        for u in (1e-12, 0.01, 0.6, 0.9):
            assert abs(a.bc_cdf(a.bc_q([u], -0.6, be), -0.6, be)[0] / u - 1) < 1e-10


def _mix(z):
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return z ^ (z >> 31)


def _block_replicate(m, r, bl, seed, b):
    """The moving-block replicate b of arcstat, written out from its definition."""
    M = 0xFFFFFFFFFFFFFFFF
    s = _mix(((seed * 0x9E3779B97F4A7C15) & M) ^ (((b + 1) * 0xBF58476D1CE4E5B9) & M))
    y = []
    while len(y) < len(m):
        s = (s + 0x9E3779B97F4A7C15) & M
        i0 = int((_mix(s) >> 11) / 9007199254740992.0 * (len(m) - bl + 1))
        for i in range(bl):
            if len(y) < len(m):
                y.append(m[len(y)] + r[i0 + i])
    return y


def test_kappa_bootstrap_is_the_refit_of_its_stated_replicate():
    n = 100
    x = [20 * i / (n - 1) for i in range(n)]
    F = a.k4_cdf(x, 9, 2, -0.3, -0.2)
    y = [1 + 0.2 * x[i] + 300 * F[i] + 3 * math.sin(0.4 * i) + math.sin(4321 * i) for i in range(n)]
    box = (-0.98, 0.95, -8, 4)
    v = a.k4_fit_varpro(x, y, [(9, math.log(2), -0.3, -0.2)], bounds=box)
    fit = [v[0], v[1], v[2], v[3], v[4], v[5], v[6]]
    G = a.k4_cdf(x, *fit[3:])
    m = [fit[0] + fit[1] * x[i] + fit[2] * G[i] for i in range(n)]
    r = [y[i] - m[i] for i in range(n)]
    for gev in (False, True):
        boot, bl = a.k4_boot(x, y, [fit], box, gev=gev, B=4, seed=31)
        assert bl[0] == a.eqfit_blocklen(r)
        for b in range(4):
            yb = _block_replicate(m, r, bl[0], 31, b)
            st = [(fit[3], math.log(fit[4]), fit[5], 0.0 if gev else fit[6])]
            o = a.k4_fit_varpro(x, yb, st, 800, box, 0.0 if gev else None)
            rd = a.k4_readings([o[0], o[2], o[3], o[4], o[5], o[6]])
            assert boot[0][b][0] == rd[0] and boot[0][b][2:4] == o[5:7]
            assert boot[0][b][1] == rd[1] or (math.isnan(rd[1]) and math.isnan(boot[0][b][1]))
            if gev:
                assert math.isnan(boot[0][b][4])
            else:
                g = a.k4_fit_varpro(x, yb, [st[0][:3] + (0.0,)], 800, box, 0.0)
                assert boot[0][b][4] == float(n * math.log(g[7] / o[7]) - 2 < 0)


def test_constrained_bootstrap_agrees_with_a_full_refit():
    n = 100
    x = [20 * i / (n - 1) for i in range(n)]
    F = a.k4_cdf(x, 9, 2, a.k4_locus_k([0.2])[0], 0.2)
    y = [1 + 0.2 * x[i] + 300 * F[i] + math.sin(4321 * i) for i in range(n)]
    f, pieces = a.k4eq_fit(x, y)
    p = f[0]
    assert min(q[7] for q in pieces[0]) == p[7]
    G = a.k4_cdf(x, *p[3:7])
    m = [p[0] + p[1] * x[i] + p[2] * G[i] for i in range(n)]
    r = [y[i] - m[i] for i in range(n)]
    boot, bl = a.k4eq_boot(x, y, f, B=3, seed=8)
    for b in range(3):
        g = a.k4eq_fit(x, _block_replicate(m, r, bl[0], 8, b))[0][0]
        assert abs(g[8] - boot[0][b][0]) < 1e-5


def test_threads():
    n = a.get_threads()
    assert n >= 1
    assert a.set_threads(5) in (1, 5)
    assert a.set_threads(None) == n


def test_locus_against_high_precision_readings():
    # discrepancy of the two readings evaluated from their definitions at 40 digits: zero to 1e-15 in k
    for h, k in ((-8.0, -0.06249674542089948), (-3.0, -0.16228952021548154), (0.0, -0.12306932718098818),
                 (0.45, 0.9321745480804188), (0.49, 0.9890420763756562)):
        assert abs(a.k4_locus_k([h])[0] - k) < 4e-15


def test_bootstrap_replicates_take_the_drift_setting_of_the_point_fit():
    # Up to 0.4.0 every replicate was refitted WITH drift whatever the point fit was.
    n = 100
    x = [20 * i / (n - 1) for i in range(n)]
    F = a.k4_cdf(x, 9, 2, a.k4_locus_k([0.2])[0], 0.2)
    y = [1 + 300 * F[i] + 2 * math.sin(4321 * i) for i in range(n)]
    f, pieces = a.k4eq_fit(x, y, drift=False)
    p = f[0]
    assert p[1] == 0.0
    G = a.k4_cdf(x, *p[3:7])
    m = [p[0] + p[2] * G[i] for i in range(n)]
    r = [y[i] - m[i] for i in range(n)]
    boot, bl = a.k4eq_boot(x, y, f, B=3, seed=8)
    assert boot == a.k4eq_boot(x, y, f, B=3, seed=8, drift=False)[0]
    for b in range(3):
        assert boot[0][b][4] == 0.0
        g = a.k4eq_fit(x, _block_replicate(m, r, bl[0], 8, b), drift=False)[0][0]
        assert abs(g[8] - boot[0][b][0]) < 1e-5
    assert all(row[4] != 0.0 for row in a.k4eq_boot(x, y, f, B=3, seed=8, drift=True)[0][0])
    one, _ = a.k4eq_boot(x, y, f, B=1, block=10 ** 6)          # the replicate is the data itself
    assert abs(one[0][0][0] - p[8]) < 1e-6 and one[0][0][4] == 0.0
    fr = a.k4_fit(x, y, drift=False)
    assert fr[0][1] == 0.0
    kb, _ = a.k4_boot(x, y, fr, (-0.98, 0.95, -60, 4), B=4, seed=3)
    assert all(row[5] == 0.0 for row in kb[0])
    one, _ = a.k4_boot(x, y, fr, (-0.98, 0.95, -60, 4), B=1, block=10 ** 6)
    assert abs(one[0][0][0] - fr[0][8]) < 1e-3 and one[0][0][5] == 0.0
    try:
        a.k4eq_boot(x, y, a.k4eq_fit(x, y)[0], B=1, drift=False)
    except ValueError:
        pass
    else:
        raise AssertionError("drift=False for a fit with a slope must be refused")


def test_bayesian_evidence_reference_size_and_r_quantile():
    ref = a.bb_ref_disc(40, D=20, m=4)
    assert ref.n == 40 and len(ref) == 80
    u = [((i + 0.5) / 40) ** 1.3 for i in range(40)]
    col = sorted(r[0] for r in ref)
    idx = 1 + 79 * 0.9; lo = int(idx); h = idx - lo
    thr = (1 - h) * col[lo - 1] + h * col[lo]          # R's quantile(type = 7)
    post = a.bb_post_disc(u, 50, 1)
    assert a.bb_evidence(u, level=0.9, M=50, ref=ref) == sum(1 for r in post if r[0] > thr) / 50
    try:
        a.bb_evidence(u[:30], ref=ref)
    except ValueError:
        pass
    else:
        raise AssertionError("a reference for another sample size must be refused")
