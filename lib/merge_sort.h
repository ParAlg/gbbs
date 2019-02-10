#pragma once
#include "utilities.h"
#include "merge.h"
#include "quicksort.h" // needed for insertion_sort

namespace pbbs {
  // not yet optimized to use moves instead of copies.
  
  // Parallel mergesort
  // This sort is stable
  // if inplace is true then the output is placed in In and Out is just used
  // as temp space.
  template <class SeqA, class SeqB, class F> 
  void merge_sort_(SeqA In, SeqB Out, const F& f, bool inplace=false) {
    size_t n = In.size();
    if (base_case(In.begin(), n)) {
      pbbs::insertion_sort(In.begin(), n, f);
      if (!inplace)
	for (size_t i=0; i < n; i++) Out[i] = In[i];
      return;
    }
    size_t m = n/2;
    par_do_if(n > 64,
	   [&] () {merge_sort_(In.slice(0,m), Out.slice(0,m), f, !inplace);},
	   [&] () {merge_sort_(In.slice(m,n), Out.slice(m,n), f, !inplace);},
	   true);
    if (inplace)
      pbbs::merge_(Out.slice(0,m), Out.slice(m,n), In, f, true);
    else
      pbbs::merge_(In.slice(0,m), In.slice(m,n), Out, f, true);
  }

  template <class SeqA, class F> 
  sequence<typename SeqA::T> merge_sort(const SeqA &In, const F& f) {
    using T = typename SeqA::T;
    sequence<T> A(In);
    sequence<T> B(In.size());
    merge_sort_(A.slice(), B.slice(), f);
    return B;
  }

  template <class T, class F>
  sequence<T> merge_sort(sequence<T> &&In, const F& f) {
    sequence<T> A(std::move(In));
    sequence<T> B(A.size());
    merge_sort_(A.slice(), B.slice(), f, true);
    return A;
  }

}
