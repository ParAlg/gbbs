#pragma once

//***************************************
// All the pbbs library uses only four functions for
// accessing parallelism.
// These can be implemented on top of any scheduler.
//***************************************
// number of threads available from OS
//template <>
int num_workers();

// id of running thread, should be numbered from [0...num-workers)
int worker_id();

void set_num_workers(int n);

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
static void par_do(Lf left, Rf right, bool conservative=false);

template <typename A, typename Af, typename Df, typename F>
static void parallel_for_alloc(Af init_alloc, Df finish_alloc, long start, long end, F f, long granularity = 0, bool conservative=false);

//***************************************

// cilkplus
#if defined(CILK)
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer.h>
#include <iostream>
#include <sstream>
#define PAR_GRANULARITY 2000

template <typename F>
inline void parallel_for(long start, long end, F f,
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
inline void parallel_for_1(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  _Pragma("cilk grainsize = 1") cilk_for(long i=start; i<end; i++) f(i);
}

template <typename Lf, typename Rf>
inline void par_do(Lf left, Rf right, bool conservative) {
    cilk_spawn right();
    left();
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
inline void parallel_for_alloc(Af init_alloc, Df finish_alloc, long start, long end, F f, long granularity, bool conservative) {
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

template <class F>
inline void parallel_for(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  _Pragma("omp parallel for")
    for(long i=start; i<end; i++) f(i);
}

template <typename F>
inline void parallel_for_1(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  #pragma omp for schedule(dynamic, 1) nowait
    for(long i=start; i<end; i++) f(i);
}

bool in_par_do = false;

template <typename Lf, typename Rf>
inline void par_do(Lf left, Rf right, bool conservative) {
  if (!in_par_do) {
    in_par_do = true;  // at top level start up tasking
#pragma omp parallel
#pragma omp single
#pragma omp task
    left();
#pragma omp task
    right();
#pragma omp taskwait
    in_par_do = false;
  } else {   // already started
#pragma omp task
    left();
#pragma omp task
    right();
  }
}

template <typename Job>
inline void parallel_run(Job job, int num_threads=0) {
  job();
}

template <typename A, typename Af, typename Df, typename F>
inline void parallel_for_alloc(Af init_alloc, Df finish_alloc, long start, long end, F f, long granularity, bool conservative) {
  A* alloc = nullptr;
  #pragma omp parallel private(alloc)
  {
    alloc = new A();
    init_alloc(alloc);
    parallel_for_1(start, end, [&](size_t i) { f(i, alloc); }, granularity, conservative);
    //#pragma omp for schedule(dynamic, 1) nowait
    //for(long i=start; i<end; i++) f(i, alloc);
    finish_alloc(alloc);
  }
}

// Guy's scheduler (ABP)
#elif defined(HOMEGROWN)
#include "scheduler.h"

#define PAR_GRANULARITY 512

template <class F>
inline void parallel_for(long start, long end, F f,
			 long granularity,
			 bool conservative) {
  global_scheduler.parfor(start, end, f, granularity, conservative);
}

template <typename Lf, typename Rf>
inline void par_do(Lf left, Rf right, bool conservative) {
  return global_scheduler.pardo(left, right, conservative);
}

template <typename Job>
inline void parallel_run(Job job, int num_threads=0) {
  job();
}

template <typename A, typename Af, typename Df, typename F>
inline void parallel_for_alloc(Af init_alloc, Df finish_alloc, long start, long end, F f, long granularity, bool conservative) {
  
  parallel_for(start, end, [&](long i)
  {
    static thread_local A* alloc = new A();
    init_alloc(alloc);
    f(i, alloc);
  }, granularity, conservative);
  //finish_alloc(alloc);
}

// c++
#else

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
inline void par_do(Lf left, Rf right, bool conservative) {
  left(); right();
}

template <typename Job>
inline void parallel_run(Job job, int num_threads=0) {
  job();
}

template <typename A, typename Af, typename Df, typename F>
inline void parallel_for_alloc(Af init_alloc, Df finish_alloc, long start, long end, F f, long granularity, bool conservative) {
  A* alloc = new A();
  init_alloc(alloc);
  for (long i=start; i<end; i++) { f(i, alloc); }
  finish_alloc(alloc);
}

#endif
