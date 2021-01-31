#pragma once
#include "binary_search.h"
#include "seq.h"
#include "utilities.h"

namespace gbbs {
namespace clustering {
namespace utils {

// TODO: not yet optimized to use moves instead of copies.

// the following parameter can be tuned
constexpr const size_t _merge_base = PAR_GRANULARITY;

template <class SeqA, class SeqB, class F>
void seq_merge(SeqA const &A, SeqB const &B,
               range<typename SeqA::value_type *> R, const F &f) {
  using T = typename SeqA::value_type;
  size_t nA = A.size();
  size_t nB = B.size();
  size_t i = 0;
  size_t j = 0;
  while (true) {
    if (i == nA) {
      while (j < nB) {
        R[i + j] = B[j];
        j++;
      }
      break;
    }
    if (j == nB) {
      while (i < nA) {
        R[i + j] = A[i];
        i++;
      }
      break;
    }
    T a = A[i];
    T b = B[j];
    if (f(b, a)) {
      R[i + j] = b;
      j++;
    } else {
      R[i + j] = a;
      i++;
    }
  }
}

// this merge is stable
template <class SeqA, class SeqB, class F>
void merge_(const SeqA &A, const SeqB &B, range<typename SeqA::value_type *> R,
            const F &f, bool cons = false) {
  size_t nA = A.size();
  size_t nB = B.size();
  size_t nR = nA + nB;
  if (nR < _merge_base)
    seq_merge(A, B, R, f);
  else if (nA == 0)
    parallel_for(0, nB, [&](size_t i) { R[i] = B[i]; });
  else if (nB == 0)
    parallel_for(0, nA, [&](size_t i) { R[i] = A[i]; });
  else {
    size_t mA = nA / 2;
    // important for stability that binary search identifies
    // first element in B greater or equal to A[mA]
    size_t mB = binary_search(B, A[mA], f);
    if (mB == 0) mA++;  // ensures at least one on each side
    size_t mR = mA + mB;
    auto left = [&]() {
      merge_(A.slice(0, mA), B.slice(0, mB), R.slice(0, mR), f, cons);
    };
    auto right = [&]() {
      merge_(A.slice(mA, nA), B.slice(mB, nB), R.slice(mR, nR), f, cons);
    };
    par_do(left, right, cons);
  }
}

template <class SeqA, class SeqB, class F>
sequence<typename SeqA::value_type> merge(const SeqA &A, const SeqB &B,
                                          const F &f, bool cons = false) {
  using T = typename SeqA::value_type;
  sequence<T> R(A.size() + B.size());
  merge_(A, B, R.slice(), f, cons);
  return R;
}

}  // namespace utils
}  // namespace clustering
}  // namespace gbbs
