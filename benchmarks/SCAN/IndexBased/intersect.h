// This file copies `intersection::intersect_f_par` from `ligra/vertex.h` and
// modifies it to get `intersect_f_with_index_par`. We copy the function here
// instead of modifying `ligra/vertex.h` because it's probably not broadly
// useful.
#pragma once

#include <functional>
#include <tuple>

#include "gbbs/bridge.h"
#include "gbbs/helpers/assert.h"
#include "gbbs/macros.h"
#include "gbbs/vertex.h"

namespace gbbs {
namespace scan {

namespace internal {

constexpr size_t _bs_merge_base = 32;
constexpr size_t _seq_merge_thresh = 2048;

template <typename Weight>
struct IntersectReturn {
  using type = Weight;  // Return sum of weights to shared neighbors
};

template <>
struct IntersectReturn<gbbs::empty> {
  using type = size_t;  // Return number of shared neighbors
};

// Return a sequence consisting of the zero-th tuple attribute of each element
// in `sequence`.
template <class Seq>
auto ProjectSequenceZero(const Seq& sequence) {
  using Element =
      typename std::tuple_element<0, typename Seq::value_type>::type;
  return parlay::delayed_seq<Element>(sequence.size(), [&](const size_t i) {
    return std::get<0>(sequence[i]);
  });
}

// see comment on `merge` for argument info
template <typename Weight, class Seq, class F>
typename IntersectReturn<Weight>::type seq_merge_full(
    const Seq& A, const Seq& B, const size_t offset_A, const size_t offset_B,
    const bool are_sequences_swapped, const F& f) {
  using ReturnType = typename IntersectReturn<Weight>::type;
  const Seq& unswapped_A = are_sequences_swapped ? B : A;
  const Seq& unswapped_B = are_sequences_swapped ? A : B;
  const size_t unswapped_offset_A = are_sequences_swapped ? offset_B : offset_A;
  const size_t unswapped_offset_B = are_sequences_swapped ? offset_A : offset_B;

  size_t nA = unswapped_A.size(), nB = unswapped_B.size();
  size_t i = 0, j = 0;
  ReturnType ct = 0;
  while (i < nA && j < nB) {
    if
      constexpr(std::is_same<Weight, gbbs::empty>::value) {
        // unweighted case
        const uintE a_id = std::get<0>(unswapped_A[i]);
        const uintE b_id = std::get<0>(unswapped_B[j]);
        if (a_id == b_id) {
          f(a_id, unswapped_offset_A + i, unswapped_offset_B + j);
          i++;
          j++;
          ct++;
        } else if (a_id < b_id) {
          i++;
        } else {
          j++;
        }
      }
    else {  // weighted case
      const auto & [ a_id, a_weight ] = unswapped_A[i];
      const auto & [ b_id, b_weight ] = unswapped_B[j];
      if (a_id == b_id) {
        f(a_id, unswapped_offset_A + i, unswapped_offset_B + j, a_weight,
          b_weight);
        i++;
        j++;
        ct += a_weight * b_weight;
      } else if (a_id < b_id) {
        i++;
      } else {
        j++;
      }
    }
  }
  return ct;
}

// see comment on `merge` for argument info
template <typename Weight, class Seq, class F>
typename IntersectReturn<Weight>::type seq_merge(
    const Seq& A, const Seq& B, const size_t offset_A, const size_t offset_B,
    const bool are_sequences_swapped, const F& f) {
  using ReturnType = typename IntersectReturn<Weight>::type;
  size_t nA = A.size();
  ReturnType ct = 0;
  const auto B_ids = ProjectSequenceZero(B);
  for (size_t i = 0; i < nA; i++) {
    if
      constexpr(std::is_same<Weight, gbbs::empty>::value) {
        // unweighted case
        const uintE a_id = std::get<0>(A[i]);
        size_t mB = parlay::binary_search(B_ids, a_id, std::less<uintE>());
        if (mB < B.size() && a_id == std::get<0>(B[mB])) {
          if (are_sequences_swapped) {
            f(a_id, offset_B + mB, offset_A + i);
          } else {
            f(a_id, offset_A + i, offset_B + mB);
          }
          ct++;
        }
      }
    else {  // weighted case
      const auto[a_id, a_weight] = A[i];
      size_t mB = parlay::binary_search(B_ids, a_id, std::less<uintE>());
      if (mB < B.size()) {
        const auto[b_id, b_weight] = B[mB];
        if (a_id == b_id) {
          if (are_sequences_swapped) {
            f(a_id, offset_B + mB, offset_A + i, b_weight, a_weight);
          } else {
            f(a_id, offset_A + i, offset_B + mB, a_weight, b_weight);
          }
          ct += a_weight * b_weight;
        }
      }
    }
  }
  return ct;
}

// `are_sequences_swapped` is true if sequence A is the neighbor list of vertex
// B and sequence B is the neighbor list of vertex A.
//
// `offset_A` is the index of the vertex's neighbor list at which sequence A
// begins, and likewise for `offset_B`.
template <typename Weight, class Seq, class F>
typename IntersectReturn<Weight>::type merge(const Seq& AA, const Seq& BB,
                                             const size_t offset_A,
                                             const size_t offset_B,
                                             const bool are_sequences_swapped,
                                             const F& f) {
  using ReturnType = typename IntersectReturn<Weight>::type;
  auto A = make_slice(AA);
  auto B = make_slice(BB);
  size_t nA = A.size();
  size_t nB = B.size();
  size_t nR = nA + nB;
  if (nR < _seq_merge_thresh) {  // handles (small, small) using linear-merge
    return seq_merge_full<Weight>(A, B, offset_A, offset_B,
                                  are_sequences_swapped, f);
  } else if (nB < nA) {
    return scan::internal::merge<Weight>(B, A, offset_B, offset_A,
                                         !are_sequences_swapped, f);
  } else if (nA < _bs_merge_base) {
    return seq_merge<Weight>(A, B, offset_A, offset_B, are_sequences_swapped,
                             f);
  } else {
    const auto B_ids = ProjectSequenceZero(B);
    size_t mA = nA / 2;
    size_t mB =
        parlay::binary_search(B_ids, std::get<0>(A[mA]), std::less<uintE>());
    ReturnType m_left = 0;
    ReturnType m_right = 0;
    par_do(
        [&]() {
          m_left = scan::internal::merge<Weight>(A.cut(0, mA), B.cut(0, mB),
                                                 offset_A, offset_B,
                                                 are_sequences_swapped, f);
        },
        [&]() {
          m_right = scan::internal::merge<Weight>(A.cut(mA, nA), B.cut(mB, nB),
                                                  offset_A + mA, offset_B + mB,
                                                  are_sequences_swapped, f);
        });
    return m_left + m_right;
  }
}

//////////////////////////////
// For unweighted vertices: //
//////////////////////////////
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
//
////////////////////////////
// For weighted vertices: //
////////////////////////////
// Returns the sum of (weight from A to C) * (weight from B to C) for shared
// neighbors C between A and B.
//
// Arguments:
//   ...
//   f: (uintE, uintE, uintE, Weight, Weight) -> void
//     Same as unweighted case, but will be run as
//     `f(c, a_to_c_index, b_to_c_index, a_to_c_weight, b_to_c_weight)`.
template <template <typename> class VertexTemplate, typename Weight, class F>
typename IntersectReturn<Weight>::type intersect_f_with_index_par(
    VertexTemplate<Weight>* A, VertexTemplate<Weight>* B, const F& f) {
  if
    constexpr(
        std::is_same<VertexTemplate<Weight>, symmetric_vertex<Weight>>::value) {
      using Neighbor = typename VertexTemplate<Weight>::edge_type;
      const auto seqA{
          gbbs::make_slice<Neighbor>(A->neighbors, A->out_degree())};
      const auto seqB{
          gbbs::make_slice<Neighbor>(B->neighbors, B->out_degree())};
      constexpr size_t kOffset{0};
      constexpr bool kAreSeqsSwapped{false};
      return scan::internal::merge<Weight>(seqA, seqB, kOffset, kOffset,
                                           kAreSeqsSwapped, f);
    }
  else {
    // TODO(tomtseng) copy over code for intersected compressed vertices
    ABORT("Not yet implemented for compressed vertices");
    (void)A;
    (void)B;
    (void)f;
  }
}

}  // namespace internal

}  // namespace scan
}  // namespace gbbs
