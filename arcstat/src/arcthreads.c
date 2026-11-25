/* arcthreads.c -- the thread setting shared by every parallel region.
 *
 * The default is at most two threads, which is what CRAN allows, and never more than OpenMP itself
 * would use. A user may ask for more (or fewer): arc_set_threads, which the fronts call with the R
 * option arcstat.threads or the Python set_threads value, and otherwise the environment variable
 * ARCSTAT_THREADS. A setting below 1 or above 1024, or one that is not a number, is ignored. */
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "arcthreads.h"

static int arc_nthreads_set = 0;                 /* 0: not set, ask the environment */

int arc_threads(void) {
#ifdef _OPENMP
  int n = arc_nthreads_set;
  if (n < 1) {
    const char *e = getenv("ARCSTAT_THREADS");
    if (e && *e) {
      char *end;
      long v = strtol(e, &end, 10);
      if (*end == '\0' && v >= 1 && v <= 1024) n = (int)v;
    }
  }
  if (n >= 1) return n;
  n = omp_get_max_threads();
  return n > 2 ? 2 : (n < 1 ? 1 : n);
#else
  return 1;
#endif
}

/* *n >= 1 sets the number, *n == 0 clears the setting, *n < 0 only asks. *now is the number in force. */
void arc_set_threads(const int *n, int *now) {
  if (*n >= 0) arc_nthreads_set = (*n <= 1024) ? *n : 0;
  *now = arc_threads();
}
