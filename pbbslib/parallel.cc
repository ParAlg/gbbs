#include "parallel.h"

#if defined(CILK)

int num_workers() { return __cilkrts_get_nworkers(); }
int worker_id() { return __cilkrts_get_worker_number(); }
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

int num_workers() { return omp_get_max_threads(); }
int worker_id() { return omp_get_thread_num(); }
void set_num_workers(int n) { omp_set_num_threads(n); }

#elif defined(HOMEGROWN)

int num_workers() { return global_scheduler.num_workers(); }

int worker_id() { return global_scheduler.worker_id(); }

void set_num_workers(int n) { global_scheduler.set_num_workers(n); }

#else

inline int num_workers() { return 1; }
inline int worker_id() { return 0; }
inline void set_num_workers(int n) { ; }

#endif
