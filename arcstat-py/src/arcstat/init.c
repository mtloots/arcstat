/* Native-routine registration for arcstat (CRAN-compliant .C interface). */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "arclen.h"
#include "arcdistc.h"
#include "bayesarc.h"
#include "cfarc.h"
#include "arccirc.h"

static const R_CMethodDef CEntries[] = {
  {"al_statistic",      (DL_FUNC) &al_statistic,      3},
  {"al_pvalue",         (DL_FUNC) &al_pvalue,         3},
  {"al_moments",        (DL_FUNC) &al_moments,        2},
  {"arcq_qd",           (DL_FUNC) &arcq_qd,           6},
  {"arcq_arclength",    (DL_FUNC) &arcq_arclength,    5},
  {"al_band_model",     (DL_FUNC) &al_band_model,     5},
  {"al_band_sample",    (DL_FUNC) &al_band_sample,    5},
  {"al_scale",          (DL_FUNC) &al_scale,          5},
  {"sample_lmoments_c", (DL_FUNC) &sample_lmoments_c, 4},
  {"arcq_fit_cf_c",     (DL_FUNC) &arcq_fit_cf_c,     3},
  {"bb_post",           (DL_FUNC) &bb_post,           5},
  {"bb_ref",            (DL_FUNC) &bb_ref,            5},
  {"cf_arclength_emp",  (DL_FUNC) &cf_arclength_emp,  5},
  {"arcc_c3max",          (DL_FUNC) &arcc_c3max,          1},
  {"arcc_qd3",            (DL_FUNC) &arcc_qd3,            4},
  {"arcc_Q3",             (DL_FUNC) &arcc_Q3,             4},
  {"arcc_dens3",          (DL_FUNC) &arcc_dens3,          5},
  {"arcc_trigmom3",       (DL_FUNC) &arcc_trigmom3,       5},
  {"arcc_rand3",          (DL_FUNC) &arcc_rand3,          5},
  {"arcc_fit3",           (DL_FUNC) &arcc_fit3,           4},
  {"arcc_qd_fr",          (DL_FUNC) &arcc_qd_fr,          6},
  {"arcc_trigmom_fr",     (DL_FUNC) &arcc_trigmom_fr,     7},
  {"arcc_Q_fr",           (DL_FUNC) &arcc_Q_fr,           7},
  {"arcc_dens_fr",        (DL_FUNC) &arcc_dens_fr,        8},
  {"arcc_rand_fr",        (DL_FUNC) &arcc_rand_fr,        8},
  {"arcc_lmom",           (DL_FUNC) &arcc_lmom,           4},
  {"arcc_rho_from_lmom",  (DL_FUNC) &arcc_rho_from_lmom,  4},
  {"arcc_qmin_rho",       (DL_FUNC) &arcc_qmin_rho,       5},
  {"arcc_admiss",         (DL_FUNC) &arcc_admiss,         6},
  {"arcc_factorise",      (DL_FUNC) &arcc_factorise,      4},
  {"arcc_fit_fr",         (DL_FUNC) &arcc_fit_fr,         6},
  {"arcc_temper_vm",     (DL_FUNC) &arcc_temper_vm,     7},
  {"arcc_temper_vm_trigmom", (DL_FUNC) &arcc_temper_vm_trigmom, 6},
  {NULL, NULL, 0}
};

void R_init_arcstat(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
