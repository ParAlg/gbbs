// This file copies `intersection::intersect_f_par` from `ligra/vertex.h` and
// modifies it to get `intersect_f_with_index_par`. We copy the function here
// instead of modifying `ligra/vertex.h` because it's probably not broadly
// useful.
#pragma once

#include <functional>
#include <tuple>

#include "gbbs/bridge.h"
#include "gbbs/macros.h"
#include "pbbslib/binary_search.h"
#include "pbbslib/seq.h"

namespace gbbs {
namespace scan {

namespace internal {

constexpr size_t _bs_merge_base = 32;
constexpr size_t _seq_merge_thresh = 2048;

// see comment on `merge` for argument info
template <class Seq, class F>
size_t seq_merge_full(
    const Seq& A,
    const Seq& B,
    const size_t offset_A,
    const size_t offset_B,
    const bool are_sequences_swapped,
    const F& f) {
  using T = typename Seq::value_type;
  const Seq& unswapped_A = are_sequences_swapped ? B : A;
  const Seq& unswapped_B = are_sequences_swapped ? A : B;
  const size_t unswapped_offset_A =
    are_sequences_swapped ? offset_B : offset_A;
  const size_t unswapped_offset_B =
    are_sequences_swapped ? offset_A : offset_B;

  size_t nA = unswapped_A.size(), nB = unswapped_B.size();
  size_t i = 0, j = 0;
  size_t ct = 0;
  while (i < nA && j < nB) {
    const T& a = unswapped_A[i];
    const T& b = unswapped_B[j];
    if (a == b) {
      f(a, unswapped_offset_A + i, unswapped_offset_B + j);
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

// see comment on `merge` for argument info
template <class Seq, class F>
size_t seq_merge(
    const Seq& A,
    const Seq& B,
    const size_t offset_A,
    const size_t offset_B,
    const bool are_sequences_swapped,
    const F& f) {
  using T = typename Seq::value_type;
  size_t nA = A.size();
  size_t ct = 0;
  for (size_t i=0; i < nA; i++) {
    const T& a = A[i];
    size_t mB = pbbslib::binary_search(B, a, std::less<T>());
    if (mB < B.size() && a == B[mB]) {
      if (are_sequences_swapped) {
        f(a, offset_B + mB, offset_A + i);
      } else {
        f(a, offset_A + i, offset_B + mB);
      }
      ct++;
    }
  }
  return ct;
}

// `are_sequences_swapped` is true if sequence A is the neighbor list of vertex
// B and sequence B is the neighbor list of vertex A.
//
// `offset_A` is the index of the vertex's neighbor list at which sequence A
// begins, and likewise for `offset_B`.
template <class Seq, class F>
size_t merge(
    const Seq& A,
    const Seq& B,
    const size_t offset_A,
    const size_t offset_B,
    const bool are_sequences_swapped,
    const F& f) {
  using T = typename Seq::value_type;
  size_t nA = A.size();
  size_t nB = B.size();
  size_t nR = nA + nB;
  if (nR < _seq_merge_thresh) {  // handles (small, small) using linear-merge
    return seq_merge_full(A, B, offset_A, offset_B, are_sequences_swapped, f);
  } else if (nB < nA) {
    return scan::internal::merge(
        B, A, offset_B, offset_A, !are_sequences_swapped, f);
  } else if (nA < _bs_merge_base) {
    return seq_merge(A, B, offset_A, offset_B, are_sequences_swapped, f);
  } else {
    size_t mA = nA/2;
    size_t mB = pbbslib::binary_search(B, A[mA], std::less<T>());
    size_t m_left = 0;
    size_t m_right = 0;
    par_do(
        [&]() {
          m_left = scan::internal::merge(
              A.slice(0, mA),
              B.slice(0, mB),
              offset_A,
              offset_B,
              are_sequences_swapped,
              f);
        },
        [&]() {
          m_right = scan::internal::merge(
              A.slice(mA, nA),
              B.slice(mB, nB),
              offset_A + mA,
              offset_B + mB,
              are_sequences_swapped,
              f);
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
//   f: (uintE, uintE, uintE) -> void
//     Function to run for each shared neighbor between A and B. It will be run
//     as `f(c, a_to_c_index, b_to_c_index)` where c is the ID of the shared
//     neighbor, `A->getOutNeighbor(a_to_c_index) == c`, and
//     `B->getOutNeighbor(b_to_c_index) == c`.
template <template <typename> class VertexTemplate, class F>
size_t intersect_f_with_index_par(
    VertexTemplate<pbbslib::empty>* A,
    VertexTemplate<pbbslib::empty>* B,
    const F& f) {
  uintT nA = A->getOutDegree(), nB = B->getOutDegree();
  uintE* nghA = reinterpret_cast<uintE*>(A->getOutNeighbors());
  uintE* nghB = reinterpret_cast<uintE*>(B->getOutNeighbors());
  auto seqA = pbbslib::make_sequence<uintE>(nghA, nA);
  auto seqB = pbbslib::make_sequence<uintE>(nghB, nB);
  constexpr size_t kOffset{0};
  constexpr bool kAreSeqsSwapped{false};
  return
    scan::internal::merge(seqA, seqB, kOffset, kOffset, kAreSeqsSwapped, f);
}

}  // namespace internal

}  // namespace scan
}  // namespace gbbs
