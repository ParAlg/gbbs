#pragma once

//***************************************
// All the pbbs library uses only four functions for
// accessing parallelism.
// These can be implemented on top of any scheduler.
//***************************************
// number of threads available from OS
int num_workers();

// id of running thread, should be numbered from [0...num-workers)
int worker_id();

// the granularity of a simple loop (e.g. adding one to each element
// of an array) to reasonably hide cost of scheduler
// #define PAR_GRANULARITY 2000

// parallel loop from start (inclusive) to end (exclusive) running
// function f.
//    f should map long to void.
//    granularity is the number of iterations to run sequentially
//      if 0 (default) then the scheduler will decide
//    conservative uses a safer scheduler
template <typename F>
static void parallel_for(long start, long end, F f,
			 long granularity = 0,
			 bool conservative = false);

// runs the thunks left and right in parallel.
//    both left and write should map void to void
//    conservative uses a safer scheduler
template <typename Lf, typename Rf>
  inline void par_do(Lf left, Rf right, bool conservative=false);

//***************************************

// cilkplus
#if defined(CILK)
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <sstream>
#define PAR_GRANULARITY 2000

int num_workers() {return __cilkrts_get_nworkers();}
int worker_id() {return __cilkrts_get_worker_number();}
void set_num_workers(int n) {
  __cilkrts_end_cilk();
  std::stringstream ss; ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}

template <typename F>
static void parallel_for(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  if (granularity == 0)
    cilk_for(long i=start; i<end; i++) f(i);
  else if ((end - start) <= granularity)
    for (long i=start; i < end; i++) f(i);
  else {
    long n = end-start;
    long mid = (start + (9*(n+1))/16);
    cilk_spawn parallel_for(start, mid, f, granularity);
    parallel_for(mid, end, f, granularity);
    cilk_sync;
  }
}

template <typename Lf, typename Rf>
  inline void par_do(Lf left, Rf right, bool conservative) {
    cilk_spawn right();
    left();
    cilk_sync;
}

template <typename Job>
static void parallel_run(Job job, int num_threads=0) {
  job();
}

// openmp
#elif defined(OPENMP)
#include <omp.h>
#define PAR_GRANULARITY 200000

int num_workers() { return omp_get_max_threads(); }
int worker_id() { return omp_get_thread_num(); }
void set_num_workers(int n) { omp_set_num_threads(n); }

template <class F>
inline void parallel_for(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  _Pragma("omp parallel for")
    for(long i=start; i<end; i++) f(i);
}

template <typename Lf, typename Rf>
static void par_do(Lf left, Rf right, bool conservative) {
#pragma omp task
    left();
#pragma omp task
    right();
#pragma omp taskwait
}

template <typename Job>
static void parallel_run(Job job, int num_threads=0) {
  job();
}

// Guy's scheduler (ABP)
#elif defined(HOMEGROWN)
#include "scheduler.h"
fork_join_scheduler fj;

// Calls fj.destroy() before the program exits
void destroy_fj() {
  fj.destroy();
}
const int atexit_result = std::atexit(destroy_fj);

#define PAR_GRANULARITY 512

int num_workers() {
  return fj.num_workers();
}

int worker_id() {
  return fj.worker_id();
}

void set_num_workers(int n) {
  fj.set_num_workers(n);
}

template <class F>
inline void parallel_for(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  fj.parfor(start, end, f, granularity, conservative);
}

template <typename Lf, typename Rf>
static void par_do(Lf left, Rf right, bool conservative) {
  return fj.pardo(left, right, conservative);
}

template <typename Job>
static void parallel_run(Job job, int num_threads=0) {
  job();
}

// c++
#else

int num_workers() { return 1;}
int worker_id() { return 0;}
void set_num_workers(int n) { ; }
#define PAR_GRANULARITY 1000

template <class F>
inline void parallel_for(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  for (long i=start; i<end; i++) {
    f(i);
  }
}

template <typename Lf, typename Rf>
static void par_do(Lf left, Rf right, bool conservative) {
  left(); right();
}

template <typename Job>
static void parallel_run(Job job, int num_threads=0) {
  job();
}

#endif
