"""Identity tests for the ellstat Python front end.

parity_ellstat.sh proves the R and Python fronts return byte-identical values from the shared C.
That is agreement, not correctness. These assert the mathematics independently, mirroring
ellstat/tests/testthat/, so the Python package is checked on its own terms.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import ellstat as e  # noqa: E402

HALF_PI = math.pi / 2.0
U = [0.3, 0.9, 1.7, 2.4]


# ------------------------------------------------------------- Jacobi elliptic functions
def test_pythagorean_identity_for_sn_and_cn():
    """sn^2 + cn^2 = 1 for every modulus. This is the defining identity of the family and it
    should hold to machine precision, not to a tolerance."""
    for k in (0.0, 0.3, 0.7, 0.95):
        sn, cn, _ = e.ell_jacobi(U, k)
        for s, c in zip(sn, cn):
            assert abs(s * s + c * c - 1.0) < 1e-14


def test_dn_identity():
    """dn^2 + k^2 sn^2 = 1, the second Jacobi identity."""
    for k in (0.0, 0.3, 0.7, 0.95):
        sn, _, dn = e.ell_jacobi(U, k)
        for s, d in zip(sn, dn):
            assert abs(d * d + k * k * s * s - 1.0) < 1e-14


def test_degenerate_modulus_reduces_to_the_circle():
    """At k = 0 the elliptic functions collapse to the trigonometric ones: sn = sin, cn = cos,
    dn = 1. This is the sense in which elliptic moments generalise trigonometric moments."""
    sn, cn, dn = e.ell_jacobi(U, 0.0)
    for u, s, c, d in zip(U, sn, cn, dn):
        assert abs(s - math.sin(u)) < 1e-12
        assert abs(c - math.cos(u)) < 1e-12
        assert abs(d - 1.0) < 1e-12


def test_complete_integrals_at_zero_modulus():
    """K(0) = E(0) = pi/2, and the nome vanishes."""
    K, E = e.ell_KE(0.0)
    assert abs(K - HALF_PI) < 1e-12
    assert abs(E - HALF_PI) < 1e-12
    assert e.ell_nome(0.0) == 0.0


def test_complete_integrals_are_ordered_and_monotone():
    """For 0 < k < 1, E(k) < pi/2 < K(k); K increases in k and E decreases."""
    Ks, Es = [], []
    for k in (0.1, 0.4, 0.7, 0.9):
        K, E = e.ell_KE(k)
        assert E < HALF_PI < K
        Ks.append(K)
        Es.append(E)
    assert all(y > x for x, y in zip(Ks, Ks[1:]))
    assert all(y < x for x, y in zip(Es, Es[1:]))


# ------------------------------------------------------------- efficiency study
def test_efficiency_gain_is_its_own_variance_ratio():
    """gain is defined as var_trig / var_elliptic, so it must equal that ratio exactly."""
    eta, var_trig, var_ell, gain, solved_t, solved_e = e.ell_efficiency(0.5, 3, 100, 50)
    assert var_trig > 0 and var_ell > 0
    assert abs(gain - var_trig / var_ell) < 1e-12
    # the solved fractions are proportions of the replicates. They exist only because the
    # buffer was widened to the five doubles the C actually writes; before that this front
    # overflowed its allocation on every call.
    assert 0.0 <= solved_t <= 1.0
    assert 0.0 <= solved_e <= 1.0


def test_efficiency_study_is_deterministic():
    a1 = e.ell_efficiency(0.5, 3, 100, 50)
    a2 = e.ell_efficiency(0.5, 3, 100, 50)
    assert a1 == a2


# ------------------------------------------------------------- projected family
def test_projected_family_is_a_probability_measure():
    """mass must be one whatever the axis ratio, or the quadrature under the efficiency is not
    integrating a density."""
    for tau in (0.5, 1.0, 2.0):
        r = e.ell_projfam_eff(tau, 1.0, kgrid=[0.0, 0.5, 0.9])
        assert abs(r["mass"] - 1.0) < 1e-6
        assert r["fisher"] > 0
        assert len(r["efficiency"]) == 3
        assert all(0.0 < x <= 1.0 + 1e-8 for x in r["efficiency"])


def test_optimal_modulus_beats_a_fine_grid():
    """ell_projfam_kopt claims the maximiser; no point on a fine grid may beat it."""
    k, eff = e.ell_projfam_kopt(1.0, 1.0)
    assert 0.0 <= k < 1.0
    grid = [i / 200.0 for i in range(200)]
    curve = e.ell_projfam_eff(1.0, 1.0, kgrid=grid)["efficiency"]
    assert eff >= max(curve) - 1e-6
    assert eff > curve[0]


def test_both_cited_families_are_reachable_and_differ():
    n = e.ell_projfam_eff(1.0, 1.0, family="normal", kgrid=[0.0, 0.9])
    c = e.ell_projfam_eff(1.0, 1.0, family="cauchy", kgrid=[0.0, 0.9])
    assert abs(n["mass"] - 1.0) < 1e-6
    assert abs(c["mass"] - 1.0) < 1e-6
    assert n["fisher"] != c["fisher"]


# ------------------------------------------------------------- projected density and sampling
def test_projected_density_is_non_negative_and_integrates_to_one():
    n = 4000
    step = 2.0 * math.pi / n
    grid = [i * step for i in range(n)]
    dens = e.dprojfam(grid, 1.0, 1.0)
    assert all(d >= 0 for d in dens)
    assert abs(sum(dens) * step - 1.0) < 1e-4


def test_projected_sampler_lands_on_the_circle_and_is_seeded():
    s1 = e.rprojfam(500, 1.0, 1.0, seed=7)
    assert len(s1) == 500
    assert all(-2 * math.pi - 1e-9 <= t <= 2 * math.pi + 1e-9 for t in s1)
    assert e.rprojfam(500, 1.0, 1.0, seed=7) == s1
    assert e.rprojfam(500, 1.0, 1.0, seed=8) != s1
