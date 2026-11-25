#ifndef ARCTHREADS_H
#define ARCTHREADS_H
/* One place decides how many OpenMP threads every parallel region of the back end uses. Results
 * never depend on the number: every parallel loop writes its own slots and selects serially. */
int arc_threads(void);
void arc_set_threads(const int *n, int *now);
#endif
