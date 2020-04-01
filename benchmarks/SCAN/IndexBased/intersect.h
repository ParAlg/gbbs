// This function copies `intersection::intersect_f_par` from `ligra/vertex.h`
// and modifies it to get `intersect_f_with_index_par`.
#pragma once

#include <functional>
#include <tuple>

#include "ligra/bridge.h"
#include "ligra/macros.h"
#include "pbbslib/binary_search.h"
#include "pbbslib/seq.h"

namespace indexed_scan {

namespace internal {

constexpr size_t _bs_merge_base = 32;
constexpr size_t _seq_merge_thresh = 2048;

template <class SeqA, class SeqB, class F>
size_t seq_merge_full(const SeqA& A, const SeqB& B, const F& f) {
  using T = typename SeqA::value_type;
  size_t nA = A.size(), nB = B.size();
  size_t i = 0, j = 0;
  size_t ct = 0;
  while (i < nA && j < nB) {
    const T& a = A[i];
    const T& b = B[j];
    if (a == b) {
      f(a, i, j);
      i++;
      j++;
      ct++;
    } else if (a < b) {
      i++;
    } else {
      j++;
    }
  }
  return ct;
}

template <class SeqA, class SeqB, class F>
size_t seq_merge(const SeqA& A, const SeqB& B, const F& f) {
  using T = typename SeqA::value_type;
  size_t nA = A.size();
  size_t ct = 0;
  for (size_t i=0; i < nA; i++) {
    const T& a = A[i];
    size_t mB = pbbslib::binary_search(B, a, std::less<T>());
    const T& b = B[mB];
    if (a == b) {
      f(a, i, mB);
      ct++;
    }
  }
  return ct;
}

template <class SeqA, class SeqB, class F>
size_t merge(const SeqA& A, const SeqB& B, const F& f) {
  using T = typename SeqA::value_type;
  size_t nA = A.size();
  size_t nB = B.size();
  size_t nR = nA + nB;
  if (nR < _seq_merge_thresh) {  // handles (small, small) using linear-merge
    return seq_merge_full(A, B, f);
  } else if (nB < nA) {
    return indexed_scan::internal::merge(B, A, f);
  } else if (nA < _bs_merge_base) {
    return seq_merge(A, B, f);
  } else {
    size_t mA = nA/2;
    size_t mB = pbbslib::binary_search(B, A[mA], std::less<T>());
    size_t m_left = 0;
    size_t m_right = 0;
    par_do(
        [&]() {
          m_left =
            indexed_scan::internal::merge(A.slice(0, mA), B.slice(0, mB), f);
        },
        [&]() {
          m_right =
            indexed_scan::internal::merge(A.slice(mA, nA), B.slice(mB, nB), f);
        });
    return m_left + m_right;
  }
}

// Returns the number of shared neighbors between `A` and `B` and runs `f` for
// each shared neighbor found. This is the same as
// `intersection::intersect_f_par` from `ligra/vertex.h` except that `f` takes
// different arguments.
//
// Arguments:
//   A
//     Vertex 1.
//   B
//     Vertex 2.
//   a
//     ID of vertex 1.
//   b
//     ID of vertex 2.
//   f: (uintE, uintE, uintE) -> void
//     Function to run for each shared neighbor between A and B. It will be run
//     as `f(c, a_to_c_index, b_to_c_index)` where c is the ID of the shared
//     neighbor, `A->getOutNeighbor(a_to_c_index) == c`, and
//     `B->getOutNeighbor(b_to_c_index) == c`. (We could include `a`, `b`, and
//     `a_to_b_index` as arguments, but instead we leave it to the user to
//     capture them into `f`.)
template <template <typename> class vertex, class F>
inline size_t intersect_f_with_index_par(
    vertex<pbbslib::empty>* A,
    vertex<pbbslib::empty>* B,
    uintE a,
    uintE b,
    const F& f) {
  uintT nA = A->getOutDegree(), nB = B->getOutDegree();
  uintE* nghA = reinterpret_cast<uintE*>(A->getOutNeighbors());
  uintE* nghB = reinterpret_cast<uintE*>(B->getOutNeighbors());
  auto seqA = pbbslib::make_sequence<uintE>(nghA, nA);
  auto seqB = pbbslib::make_sequence<uintE>(nghB, nB);
  return indexed_scan::internal::merge(seqA, seqB, f);
}

}  // namespace internal

}  // namespace indexed_scan
