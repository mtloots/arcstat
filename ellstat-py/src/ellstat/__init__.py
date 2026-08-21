"""ellstat: the elliptic moment system for circular data.

A modulus-indexed deformation of the trigonometric moment system, built on the
Jacobi elliptic functions, which reduces to the trigonometric system exactly at
modulus zero, characterises every circular law at every modulus, and carries
exact finite-sample first- and second-order theory. ctypes binding to a shared
C back end also used by the R package 'ellstat'. Pure standard library."""
from ._core import (ell_KE, ell_nome, ell_jacobi, ell_basis, ell_moments,
                    ell_moment_var, ell_transfer, ell_dn2_moments,
                    dellvm, rellvm, ell_evm_moments, fit_ellvm,
                    ell_exact_uniform, ell_exact_ci, ell_efficiency,
                    dprojell, rprojell, ell_eff_curve, ell_kopt,
                    dprojfam, rprojfam, ell_projfam_eff, ell_catalogue,
                    ell_joint, ell_joint_sample, ell_joint_uniform)
__all__ = ["ell_KE", "ell_nome", "ell_jacobi", "ell_basis", "ell_moments",
           "ell_moment_var", "ell_transfer", "ell_dn2_moments", "dellvm",
           "rellvm", "ell_evm_moments", "fit_ellvm", "ell_exact_uniform",
           "ell_exact_ci", "ell_efficiency", "dprojell", "rprojell",
           "ell_eff_curve", "ell_kopt", "dprojfam", "rprojfam",
           "ell_projfam_eff", "ell_catalogue", "ell_joint",
           "ell_joint_sample", "ell_joint_uniform"]
__version__ = "0.1.0"
