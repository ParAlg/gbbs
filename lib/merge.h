#pragma once
#include "utilities.h"
#include "binary_search.h"

namespace pbbs {

  // the following parameter can be tuned
  constexpr const size_t _merge_base = PAR_GRANULARITY; 

  template <class SeqA, class SeqB, class SeqR, class F> 
  void seq_merge(SeqA A, SeqB B, SeqR R, const F& f) {
    using T = typename SeqA::T;
    size_t nA = A.size();
    size_t nB = B.size();
    size_t i = 0;
    size_t j = 0;
    while (true) {
      if (i == nA) {
	while (j < nB) {R.update(i+j, B[j]); j++;}
	break; }
      if (j == nB) {
	while (i < nA) {R.update(i+j, A[i]); i++;}
	break; }
      T a = A[i];
      T b = B[j];
      if (f(a, b)) {
	R.update(i+j, a); 
	i++;
      } else {
	R.update(i+j, b); 
	j++;
      }
    }
  }
    
  template <class SeqA, class SeqB, class SeqR, class F> 
  void par_merge(SeqA A, SeqB B, SeqR R, const F& f, bool cons=false) {
    size_t nA = A.size();
    size_t nB = B.size();
    size_t nR = nA + nB;
    if (nR < _merge_base) seq_merge(A, B, R, f);
    else if (nB > nA) par_merge(B, A, R, f); 
    else {
      size_t mA = nA/2;
      size_t mB = binary_search(B, A[mA], f);
      size_t mR = mA + mB;
      auto left = [&] () {par_merge(A.slice(0, mA), B.slice(0, mB), 
				    R.slice(0, mR), f, cons);};
      auto right = [&] () {par_merge(A.slice(mA, nA), B.slice(mB, nB), 
				     R.slice(mR, nR), f, cons);};
      par_do(left, right, cons);
    }
  }

  template <class SeqA, class SeqB, class SeqR, class F> 
  void merge(SeqA A, SeqB B, SeqR R, const F& f, bool cons=false) {
#if defined(OPENMP)
#pragma omp parallel
#pragma omp single
#endif
    par_merge(A, B, R, f, cons);
  }
}
