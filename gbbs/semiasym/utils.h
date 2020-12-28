#pragma once

#include "gbbs/macros.h"

namespace gbbs {

/* Note that the block degree is computable by differencing two starts. */
struct vtx_info {
  /* vertex's (packed) degree. */
  uintE vtx_degree;
  /* number of blocks associated with v */
  uintE vtx_num_blocks;
  /* pointer into the block structure */
  size_t vtx_block_offset;

  vtx_info(uintE degree, uintE num_blocks, size_t block_offset)
      : vtx_degree(degree),
        vtx_num_blocks(num_blocks),
        vtx_block_offset(block_offset) {}
};

// A sequence which provides value references for [], but accesses these
// values through a provided function. The function returns a pointer
// (indirect) reference to the desired location.
template <typename T, typename F>
struct indirect_value_sequence {
  using value_type = T;
  indirect_value_sequence(size_t n, F _f) : f(_f), s(0), e(n) {};
  indirect_value_sequence(size_t n, value_type v) : f([&] (size_t i) {return v;}), s(0), e(n) {};
  indirect_value_sequence(size_t s, size_t e, F _f) : f(_f), s(s), e(e) {};
  value_type& operator[] (size_t i) const {return *((f)(i+s));}
  indirect_value_sequence<T,F> slice(size_t ss, size_t ee) const {
    return indirect_value_sequence<T,F>(s+ss,s+ee,f); }
  indirect_value_sequence<T,F> slice() const {
    return indirect_value_sequence<T,F>(s,e,f); }
  size_t size() const { return e - s;}
private:
  const F f;
  const size_t s, e;
};

// used so second template argument can be inferred
template <class T, class F>
indirect_value_sequence<T,F> indirect_value_seq (size_t n, F f) {
  return indirect_value_sequence<T,F>(n,f);
}

}  // namespace gbbs
