#pragma once

//***************************************
// All the pbbs library uses only four functions for
// accessing parallelism.
// These can be implemented on top of any scheduler.
//***************************************
// number of threads available from OS
//template <>
static int num_workers();

// id of running thread, should be numbered from [0...num-workers)
static int worker_id();

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
static void parallel_for(long start, long end, F&& f,
			 long granularity = 0,
			 bool conservative = false);

// runs the thunks left and right in parallel.
//    both left and write should map void to void
//    conservative uses a safer scheduler
template <typename Lf, typename Rf>
static void par_do(Lf&& left, Rf&& right, bool conservative=false);

template <typename A, typename Af, typename Df, typename F>
static void parallel_for_alloc(Af&& init_alloc, Df&& finish_alloc, long start, long end, F&& f, long granularity = 0, bool conservative=false);

//***************************************

// cilkplus
#if defined(CILK)
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer.h>
#include <iostream>
#include <sstream>
#define PAR_GRANULARITY 2000

inline int num_workers() {return __cilkrts_get_nworkers();}
inline int worker_id() {return __cilkrts_get_worker_number();}
inline void set_num_workers(int n) {
  __cilkrts_end_cilk();
  std::stringstream ss; ss << n;
  if (0 != __cilkrts_set_param("nworkers", ss.str().c_str())) {
    std::cerr << "failed to set worker count!" << std::endl;
    std::abort();
  }
}

template <typename F>
inline void parallel_for(long start, long end, F&& f,
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

template <typename F>
inline void parallel_for_1(long start, long end, F&& f,
			 long granularity,
			 bool conservative) {
  _Pragma("cilk grainsize = 1") cilk_for(long i=start; i<end; i++) f(i);
}

template <typename Lf, typename Rf>
inline void par_do(Lf&& left, Rf&& right, bool conservative) {
    cilk_spawn std::forward<Rf>(right)();
    std::forward<Lf>(left)();
    cilk_sync;
}

template <typename A>
class alloc_holder
{
   struct Monoid: cilk::monoid_base<A>
   {
     static void reduce (A *left, A *right) {}
   };

public:
  cilk::reducer<Monoid> imp_;
  alloc_holder() : imp_() { }
};

// TODO try parallel_for_1
template <typename A, typename Af, typename Df, typename F>
inline void parallel_for_alloc(Af&& init_alloc, Df&& finish_alloc, long start, long end, F&& f, long granularity, bool conservative) {
  alloc_holder<A> alloc;

  parallel_for_1(start, end, [&](size_t i)
  {
    init_alloc(&alloc.imp_.view());
    f(i, &(alloc.imp_.view()));
    //finish_alloc(&(alloc.imp_.view()));
  }, granularity, conservative);
}

// openmp
#elif defined(OPENMP)
#include <omp.h>
#define PAR_GRANULARITY 200000

inline int num_workers() { return omp_get_max_threads(); }
inline int worker_id() { return omp_get_thread_num(); }
inline void set_num_workers(int n) { omp_set_num_threads(n); }

template <class F>
inline void parallel_for(long start, long end, F&& f,
			 long granularity,
			 bool conservative) {
  _Pragma("omp parallel for")
    for(long i=start; i<end; i++) f(i);
}

template <typename F>
inline void parallel_for_1(long start, long end, F&& f,
			 long granularity,
			 bool conservative) {
  #pragma omp for schedule(dynamic, 1) nowait
    for(long i=start; i<end; i++) f(i);
}

bool in_par_do = false;

template <typename Lf, typename Rf>
inline void par_do(Lf&& left, Rf&& right, bool conservative) {
  if (!in_par_do) {
    in_par_do = true;  // at top level start up tasking
#pragma omp parallel
#pragma omp single
#pragma omp task
    std::forward<Lf>(left)();
#pragma omp task
    std::forward<Rf>(right)();
#pragma omp taskwait
    in_par_do = false;
  } else {   // already started
#pragma omp task
    std::forward<Lf>(left)();
#pragma omp task
    std::forward<Rf>(right)();
  }
}

template <typename Job>
inline void parallel_run(Job&& job, int num_threads=0) {
  std::forward<Job>(job)();
}

template <typename A, typename Af, typename Df, typename F>
inline void parallel_for_alloc(Af&& init_alloc, Df&& finish_alloc, long start, long end, F&& f, long granularity, bool conservative) {
  A* alloc = nullptr;
  #pragma omp parallel private(alloc)
  {
    alloc = new A();
    std::forward<Af>(init_alloc)(alloc);
    parallel_for_1(start, end, [&](size_t i) { f(i, alloc); }, granularity, conservative);
    //#pragma omp for schedule(dynamic, 1) nowait
    //for(long i=start; i<end; i++) f(i, alloc);
    std::forward<Df>(finish_alloc)(alloc);
  }
}

// Guy's scheduler (ABP)
#elif defined(HOMEGROWN)
#include "scheduler.h"

#ifdef NOTMAIN
extern fork_join_scheduler fj;
#else
fork_join_scheduler fj;
#endif

// Calls fj.destroy() before the program exits
inline void destroy_fj() {
  fj.destroy();
}

struct __atexit {__atexit() {std::atexit(destroy_fj);}};
static __atexit __atexit_var;

#define PAR_GRANULARITY 512

inline int num_workers() {
  return fj.num_workers();
}

inline int worker_id() {
  return fj.worker_id();
}

inline void set_num_workers(int n) {
  fj.set_num_workers(n);
}

template <class F>
inline void parallel_for(long start, long end, F&& f,
			 long granularity,
			 bool conservative) {
  fj.parfor(start, end, std::forward<F>(f), granularity, conservative);
}

template <typename Lf, typename Rf>
inline void par_do(Lf&& left, Rf&& right, bool conservative) {
  return fj.pardo(std::forward<Lf>(left), std::forward<Rf>(right), conservative);
}

template <typename Job>
inline void parallel_run(Job&& job, int num_threads=0) {
  std::forward<Job>(job)();
}

template <typename A, typename Af, typename Df, typename F>
inline void parallel_for_alloc(Af&& init_alloc, Df&& finish_alloc, long start, long end, F&& f, long granularity, bool conservative) {

  parallel_for(start, end, [&](long i)
  {
    static thread_local A* alloc = new A();
    init_alloc(alloc);
    f(i, alloc);
  }, granularity, conservative);
  //std::forward<Df>(finish_alloc)(alloc);
}

// c++
#else

inline int num_workers() { return 1;}
inline int worker_id() { return 0;}
inline void set_num_workers(int n) { ; }
#define PAR_GRANULARITY 1000

template <class F>
inline void parallel_for(long start, long end, F&& f,
			 long granularity,
			 bool conservative) {
  for (long i=start; i<end; i++) {
    f(i);
  }
}

template <typename Lf, typename Rf>
inline void par_do(Lf&& left, Rf&& right, bool conservative) {
  std::forward<Lf>(left)();
  std::forward<Rf>(right)();
}

template <typename Job>
inline void parallel_run(Job&& job, int num_threads=0) {
  std::forward<Job>(job)();
}

template <typename A, typename Af, typename Df, typename F>
inline void parallel_for_alloc(Af&& init_alloc, Df&& finish_alloc, long start, long end, F f, long granularity, bool conservative) {
  A* alloc = new A();
  std::forward<Af>(init_alloc)(alloc);
  for (long i=start; i<end; i++) { f(i, alloc); }
  std::forward<Df>(finish_alloc)(alloc);
}

#endif
