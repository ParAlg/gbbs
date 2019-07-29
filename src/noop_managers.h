#pragma once

#include "encodings/byte_pd_amortized.h"
#include "macros.h"

template <template <class W> class vertex, class W>
struct sym_noop_manager {
  using WV = vertex<W>;
  using E = typename WV::E;

  uintE vtx_id;
  size_t degree;
  size_t block_size;

#ifndef NVM
  E* e0;
  sym_noop_manager(vertex<W>& V, const uintE vtx_id, const uintE block_size)
      : vtx_id(vtx_id),
        degree(V.getOutDegree()),
        block_size(block_size),
        e0(V.getOutNeighbors()) {}
#else
  E* e0;
  E* e1;
  sym_noop_manager(vertex<W>& V0, vertex<W>& V1, const uintE vtx_id,
                   const uintE block_size)
      : vtx_id(vtx_id),
        degree(V0.getOutDegree()),
        block_size(block_size),
        e0(V0.getOutNeighbors),
        e1(V1.getOutNeighbors) {}
#endif

  __attribute__((always_inline)) inline uintE get_degree() { return degree; }

  __attribute__((always_inline)) inline size_t num_blocks() {
    if (degree > 0) {
      return pbbs::num_blocks(degree, block_size);
    }
    return 0;
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_num) {
    uintE block_start = block_num * block_size;
    uintE block_end = std::min(block_start + block_size, degree);
    return block_end - block_start;
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block(uintE block_num,
                                                          F f) {
    uintE block_start = block_num * block_size;
    uintE block_end = std::min(block_start + block_size, degree);
    E* e = e0;

    for (size_t i = block_start; i < block_end; i++) {
      auto& ee = e[i];
      f(std::get<0>(ee), std::get<1>(ee), i);
    }
  }

  template <class F>
  __attribute__((always_inline)) inline void decode_block_cond(uintE block_num,
                                                               F f) {
    uintE block_start = block_num * block_size;
    uintE block_end = std::min(block_start + block_size, degree);
    E* e = e0;

    for (size_t i = block_start; i < block_end; i++) {
      auto& ee = e[i];
      bool ret = f(std::get<0>(ee), std::get<1>(ee), i);
      if (!ret) break;
    }
  }

  std::tuple<uintE, W> ith_neighbor(size_t i) {
    E* edges = e0;
#ifdef NVM
    if (numanode() == 0)
      edges = e0;
    else
      edges = e1;
#endif
    return edges[i];
  }
};

template <template <class W> class vertex, class W>
struct compressed_sym_noop_manager {
  using WV = vertex<W>;
  using E = typename WV::E;

  uintE vtx_id;
  size_t degree;

#ifndef NVM
  E* e0;
  compressed_sym_noop_manager(vertex<W>& V, const uintE vtx_id)
      : vtx_id(vtx_id), degree(V.getOutDegree()), e0(V.getOutNeighbors()) {}
#else
  E* e0;
  E* e1;
  compressed_sym_noop_manager(vertex<W>& V0, vertex<W>& V1, const uintE vtx_id)
      : vtx_id(vtx_id),
        degree(V0.getOutDegree()),
        e0(V0.getOutNeighbors),
        e1(V1.getOutNeighbors) {}
#endif

  __attribute__((always_inline)) inline uintE get_degree() { return degree; }

  __attribute__((always_inline)) inline size_t num_blocks() {
    return bytepd_amortized::get_num_blocks(e0, degree);  // TODO: nvm
  }

  __attribute__((always_inline)) inline uintE block_degree(uintE block_num) {
    return bytepd_amortized::get_block_degree(e0, degree,
                                              block_num);  // TODO: nvm
  }

  template <class F>
  inline void decode_block(uintE block_num, F f) {
    bytepd_amortized::template decode_block<W>(f, e0, vtx_id, degree,
                                               block_num);
  }

  template <class F>
  inline void decode_block_cond(uintE block_num, F f) {
    bytepd_amortized::template decode_block_cond<W>(f, e0, vtx_id, degree,
                                                    block_num);
  }
};
