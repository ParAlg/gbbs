#include "parallel.h"

#if defined(CILK)

namespace pbbs {
void set_num_workers(int n) {
  __cilkrts_end_cilk();
  std::stringstream ss;
  ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}
}  // namespace pbbs

#elif defined(OPENMP)

namespace pbbs {
bool in_par_do = false;
void set_num_workers(int n) { omp_set_num_threads(n); }
}  // namespace pbbs

#elif defined(HOMEGROWN)

namespace pbbs {
void set_num_workers(int n) { pbbs::global_scheduler.set_num_workers(n); }
}  // namespace pbbs

#else

namespace pbbs {
void set_num_workers(int n) { ; }
}  // namespace pbbs

#endif

