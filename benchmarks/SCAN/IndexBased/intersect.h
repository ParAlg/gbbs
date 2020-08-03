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

// are_vertices_swapped:  A_swapped = vertex B, B_swapped = vertex A
// !are_vertices_swapped: A_swapped = vertex A, B_swapped = vertex B
template <class Seq, class F>
double seq_merge_full(
    const Seq& A_swapped,
    const Seq& B_swapped,
    const size_t offsetA_swapped,
    const size_t offsetB_swapped,
    const F& f,
    bool are_vertices_swapped) {
  using T = typename Seq::value_type;
  const Seq& A = are_vertices_swapped ? B_swapped : A_swapped;
  const Seq& B = are_vertices_swapped ? A_swapped : B_swapped;
  const size_t offsetA = are_vertices_swapped ? offsetB_swapped : offsetA_swapped;
  const size_t offsetB = are_vertices_swapped ? offsetA_swapped : offsetB_swapped;

  size_t nA = A.size(), nB = B.size();
  size_t i = 0, j = 0;
  double ct = 0;
  while (i < nA && j < nB) {
    const T& a = A[i];
    const T& b = B[j];
    if (std::get<0>(a) == std::get<0>(b)) {
      f(std::get<0>(a), offsetA + i, offsetB + j, std::get<1>(a), std::get<1>(b));
      i++;
      j++;
      ct += std::get<1>(a) * std::get<1>(b);
    } else if (std::get<0>(a) < std::get<0>(b)) {
      i++;
    } else {
      j++;
    }
  }
  return ct;
}

// are_vertices_swapped:  A = original first vertex, B = original second vertex
// !are_vertices_swapped: B = original first vertex, A = original second vertex
template <class Seq, class F>
double seq_merge(
    const Seq& A,
    const Seq& B,
    const size_t offsetA,
    const size_t offsetB,
    const F& f,
    const bool are_vertices_swapped) {
  using T = typename Seq::value_type;
  size_t nA = A.size();
  double ct = 0;
  const auto B0 = pbbslib::make_sequence<uintE>(B.size(), [&](size_t j) { return std::get<0>(B[j]); });
  for (size_t i=0; i < nA; i++) {
    const T& a = A[i];
    size_t mB = pbbslib::binary_search(B0, std::get<0>(a), std::less<uintE>());
    if (mB < B.size()) {
      const T& b = B[mB];
      if (std::get<0>(a) == std::get<0>(b)) {
        if (are_vertices_swapped) {
          f(std::get<0>(a), offsetB + mB, offsetA + i, std::get<1>(b), std::get<1>(a));
        } else {
          f(std::get<0>(a), offsetA + i, offsetB + mB, std::get<1>(a), std::get<1>(b));
        }
        ct += std::get<1>(a) * std::get<1>(b);
      }
    }
  }
  return ct;
}

// are_vertices_swapped:  A = original first vertex, B = original second vertex
// !are_vertices_swapped: B = original first vertex, A = original second vertex
template <class Seq, class F>
double merge(
    const Seq& A,
    const Seq& B,
    const size_t offsetA,
    const size_t offsetB,
    const F& f,
    const bool are_vertices_swapped) {
  size_t nA = A.size();
  size_t nB = B.size();
  size_t nR = nA + nB;
  if (nR < _seq_merge_thresh) {  // handles (small, small) using linear-merge
    return scan::internal::seq_merge_full(A, B, offsetA, offsetB, f, are_vertices_swapped);
  } else if (nB < nA) {
    return scan::internal::merge(B, A, offsetB, offsetA, f, !are_vertices_swapped);
  } else if (nA < _bs_merge_base) {
    return scan::internal::seq_merge(A, B, offsetA, offsetB, f, are_vertices_swapped);
  } else {
    size_t mA = nA/2;
    size_t mB = pbbslib::binary_search(
        pbbslib::make_sequence<uintE>(B.size(), [&](size_t i) { return std::get<0>(B[i]); }),
        std::get<0>(A[mA]),
        std::less<uintE>());
    double m_left = 0;
    double m_right = 0;
    par_do(
        [&]() {
          m_left = scan::internal::merge(
              A.slice(0, mA),
              B.slice(0, mB),
              offsetA,
              offsetB,
              f,
              are_vertices_swapped);
        },
        [&]() {
          m_right = scan::internal::merge(
              A.slice(mA, nA),
              B.slice(mB, nB),
              offsetA + mA,
              offsetB + mB,
              f,
              are_vertices_swapped);
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
//     capture them into `f` if so desired.)
template <class Vertex, class F>
double intersect_f_with_index_par(
    Vertex* A,
    Vertex* B,
    const F& f) {
  using N = typename Vertex::edge_type;
  uintT nA = A->getOutDegree(), nB = B->getOutDegree();
  N* nghA = A->getOutNeighbors();
  N* nghB = B->getOutNeighbors();
  auto seqA = pbbslib::make_sequence<N>(nghA, nA);
  auto seqB = pbbslib::make_sequence<N>(nghB, nB);
  constexpr size_t offsetA = 0;
  constexpr size_t offsetB = 0;
  constexpr bool kAreVerticesSwapped = false;
  return scan::internal::merge(seqA, seqB, offsetA, offsetB, f, kAreVerticesSwapped);
}

}  // namespace internal

}  // namespace scan
}  // namespace gbbs
