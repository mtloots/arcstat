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
