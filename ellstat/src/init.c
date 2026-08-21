#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <stdlib.h>

extern void C_ell_ke(void *, void *, void *, void *);
extern void C_ell_nome(void *, void *, void *);
extern void C_ell_jac(void *, void *, void *, void *, void *, void *);
extern void C_ell_basis(void *, void *, void *, void *, void *, void *);
extern void C_ell_smom(void *, void *, void *, void *, void *);
extern void C_ell_transfer(void *, void *, void *, void *);
extern void C_ell_dn2mom(void *, void *, void *);
extern void C_ell_evm_d(void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_evm_r(void *, void *, void *, void *, void *, void *);
extern void C_ell_moment_var(void *, void *, void *, void *, void *);
extern void C_ell_popmom_grid(void *, void *, void *, void *, void *);
extern void C_ell_evm_pop(void *, void *, void *, void *, void *, void *);
extern void C_ell_evm_fit(void *, void *, void *, void *, void *, void *);
extern void C_ell_exact_unif(void *, void *, void *, void *, void *, void *);
extern void C_ell_exact_ci(void *, void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *);
extern void C_ell_effstudy(void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_pj_dens(void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_pj_rand(void *, void *, void *, void *, void *, void *);
extern void C_ell_pj_effcurve(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_pj_kopt(void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_projfam_dens(void *, void *, void *, void *, void *, void *);
extern void C_ell_projfam_rand(void *, void *, void *, void *, void *, void *);
extern void C_ell_projfam_eff(void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_projfam_kopt(void *, void *, void *, void *, void *, void *);
extern void C_ell_cat_eff(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_cat_kopt(void *, void *, void *, void *, void *, void *, void *);
extern void C_ell_joint_cov(void *, void *, void *, void *, void *, void *);
extern void C_ell_joint_samp(void *, void *, void *, void *, void *, void *);
extern void C_ell_joint_unif(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"C_ell_ke",          (DL_FUNC) &C_ell_ke,           4},
    {"C_ell_nome",        (DL_FUNC) &C_ell_nome,         3},
    {"C_ell_jac",         (DL_FUNC) &C_ell_jac,          6},
    {"C_ell_basis",       (DL_FUNC) &C_ell_basis,        6},
    {"C_ell_smom",        (DL_FUNC) &C_ell_smom,         5},
    {"C_ell_transfer",    (DL_FUNC) &C_ell_transfer,     4},
    {"C_ell_dn2mom",      (DL_FUNC) &C_ell_dn2mom,       3},
    {"C_ell_evm_d",       (DL_FUNC) &C_ell_evm_d,        7},
    {"C_ell_evm_r",       (DL_FUNC) &C_ell_evm_r,        6},
    {"C_ell_moment_var",  (DL_FUNC) &C_ell_moment_var,   5},
    {"C_ell_popmom_grid", (DL_FUNC) &C_ell_popmom_grid,  5},
    {"C_ell_evm_pop",     (DL_FUNC) &C_ell_evm_pop,      6},
    {"C_ell_evm_fit",     (DL_FUNC) &C_ell_evm_fit,      6},
    {"C_ell_exact_unif",  (DL_FUNC) &C_ell_exact_unif,   6},
    {"C_ell_exact_ci",    (DL_FUNC) &C_ell_exact_ci,    11},
    {"C_ell_effstudy", (DL_FUNC) &C_ell_effstudy,  7},
    {"C_ell_pj_dens",     (DL_FUNC) &C_ell_pj_dens,      7},
    {"C_ell_pj_rand",     (DL_FUNC) &C_ell_pj_rand,      6},
    {"C_ell_pj_effcurve", (DL_FUNC) &C_ell_pj_effcurve,  9},
    {"C_ell_pj_kopt",     (DL_FUNC) &C_ell_pj_kopt,      7},
    {"C_ell_projfam_dens", (DL_FUNC) &C_ell_projfam_dens,  6},
    {"C_ell_projfam_rand", (DL_FUNC) &C_ell_projfam_rand,  6},
    {"C_ell_projfam_eff",  (DL_FUNC) &C_ell_projfam_eff,   8},
    {"C_ell_projfam_kopt", (DL_FUNC) &C_ell_projfam_kopt,  6},
    {"C_ell_cat_eff",  (DL_FUNC) &C_ell_cat_eff,   9},
    {"C_ell_cat_kopt", (DL_FUNC) &C_ell_cat_kopt,  7},
    {"C_ell_joint_cov",  (DL_FUNC) &C_ell_joint_cov,   6},
    {"C_ell_joint_samp", (DL_FUNC) &C_ell_joint_samp,  6},
    {"C_ell_joint_unif", (DL_FUNC) &C_ell_joint_unif,  8},
    {NULL, NULL, 0}
};

void attribute_visible R_init_ellstat(DllInfo *dll){
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
