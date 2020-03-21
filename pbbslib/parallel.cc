#include "parallel.h"

#if defined(CILK)

void set_num_workers(int n) {
  __cilkrts_end_cilk();
  std::stringstream ss;
  ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}

#elif defined(OPENMP)

bool in_par_do = false;

void set_num_workers(int n) { omp_set_num_threads(n); }

#elif defined(HOMEGROWN)

void set_num_workers(int n) { global_scheduler.set_num_workers(n); }

#else

void set_num_workers(int n) { ; }

#endif
