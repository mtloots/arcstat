"""arcstat: arc-length statistics -- goodness of fit, distributions and a Bayesian test. ctypes
binding to a shared C back-end (compiled on first import from the bundled sources; the same sources
back the R package 'arcstat'). Pure standard library."""
# goodness of fit
from ._core import (al_band_model, al_band_sample, al_moments, al_pvalue, al_scale, al_statistic)
# the arcq distribution family and L-moments
from ._core import (arcq_arclength, arcq_qd, fit_arcq_cf, sample_lmoments)
# the Bayesian arc-length test
from ._core import (bb_evidence, bb_post_disc, bb_ref_disc)
# characteristic-function arc length
from ._core import (cf_arclength, cf_arclength_family)
# the circular arc-length family
from ._core import (arcc_Q3, arcc_Q_fr, arcc_admiss, arcc_c3max, arcc_dens3, arcc_dens_fr, arcc_factorise, arcc_fit3, arcc_fit_fr, arcc_lmom, arcc_qd3, arcc_qd_fr, arcc_qmin_rho, arcc_rand3, arcc_rand_fr, arcc_rho3, arcc_rho_from_lmom, arcc_temper_vm, arcc_temper_vm_trigmom, arcc_trigmom3, arcc_trigmom_fr)
# the beta-companion equivalence family
from ._core import (bc_cdf, bc_pdf, bc_q, eq_E, eq_bstar, eq_readings, eq_readings_vsl, eq_ub_quad)
# the kappa-4 induction curve
from ._core import (k4_band_model, k4_band_sample, k4_cdf, k4_fit_aleq, k4_fit_lmom, k4_fit_nalr, k4_fit_nls, k4_pdf, k4_q, k4_readings, k4_runmed, k4_tau34)

__all__ = [
    "al_band_model",
    "al_band_sample",
    "al_moments",
    "al_pvalue",
    "al_scale",
    "al_statistic",
    "arcc_Q3",
    "arcc_Q_fr",
    "arcc_admiss",
    "arcc_c3max",
    "arcc_dens3",
    "arcc_dens_fr",
    "arcc_factorise",
    "arcc_fit3",
    "arcc_fit_fr",
    "arcc_lmom",
    "arcc_qd3",
    "arcc_qd_fr",
    "arcc_qmin_rho",
    "arcc_rand3",
    "arcc_rand_fr",
    "arcc_rho3",
    "arcc_rho_from_lmom",
    "arcc_temper_vm",
    "arcc_temper_vm_trigmom",
    "arcc_trigmom3",
    "arcc_trigmom_fr",
    "arcq_arclength",
    "arcq_qd",
    "bb_evidence",
    "bb_post_disc",
    "bb_ref_disc",
    "bc_cdf",
    "bc_pdf",
    "bc_q",
    "cf_arclength",
    "cf_arclength_family",
    "eq_E",
    "eq_bstar",
    "eq_readings",
    "eq_readings_vsl",
    "eq_ub_quad",
    "fit_arcq_cf",
    "k4_band_model",
    "k4_band_sample",
    "k4_cdf",
    "k4_fit_aleq",
    "k4_fit_lmom",
    "k4_fit_nalr",
    "k4_fit_nls",
    "k4_pdf",
    "k4_q",
    "k4_readings",
    "k4_runmed",
    "k4_tau34",
    "sample_lmoments",
]
__version__ = "0.1.0"
