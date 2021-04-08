#pragma once
#include <optional>
#include <string.h>

// *******************************************
//   Utils
// *******************************************

namespace cpam {
namespace utils {

  //#define LOG(x) (std::cout << x << std::endl)
  #define LOG(x)

  // The block size to use for compressed trees. TODO: we may want to make this
  // a runtime or template argument for maps.
  constexpr const size_t compression_block_size = 60;
  //constexpr const size_t compression_block_size = 35;

  static constexpr size_t kBaseCaseSize = 7*compression_block_size;

  // for granularity control
  // if the size of a node is smaller than this number then it
  // processed sequentially instead of in parallel
  constexpr const size_t node_limit = compression_block_size*4;

  template <class E>
  E* new_array_no_init(size_t n) {
#ifndef PARLAY_USE_STD_ALLOC
    auto allocator = parlay::allocator<E>();
#else
    auto allocator = std::allocator<E>();
#endif
    return allocator.allocate(n);
  }

  template <class E>
  void free_array(E* e, size_t n) {
#ifndef PARLAY_USE_STD_ALLOC
    auto allocator = parlay::allocator<E>();
#else
    auto allocator = std::allocator<E>();
#endif
    allocator.deallocate(e, n);
  }

  // for two input sizes of n and m, should we do a parallel fork
  // assumes work proportional to m log (n/m + 1)
  static bool do_parallel(size_t n, size_t m) {
    if (m > n) std::swap(n,m);
    return (m > 8 && (m * parlay::log2_up(n/m + 1)) > node_limit);
  }

  // fork-join parallel call, returning a pair of values
  template <class RT, class Lf, class Rf>
  static std::pair<RT,RT> fork(bool do_parallel, Lf left, Rf right) {
    if (do_parallel) { //do_parallel) {
      RT r, l;
      auto do_right = [&] () {r = right();};
      auto do_left = [&] () {l = left();};
      parlay::par_do(do_left, do_right);
      return std::pair<RT,RT>(l,r);
    } else {
      RT l = left();
      RT r = right();
      return std::make_pair(l,r);
    }
  }

  // fork-join parallel call, returning nothing
  template <class Lf, class Rf>
  static void fork_no_result(bool do_parallel, Lf left, Rf right) {
    if (do_parallel) parlay::par_do(left, right);
    else {left(); right();}
  }

  template<class V>
  struct get_left {
    V operator () (const V& a, const V& b) const
    { return a; }
  };

  template<class V>
  struct get_right {
    V operator () (const V& a, const V& b) const
    { return b; }
  };

  // should update the following to use more portable implementation
  template <typename ET>
  inline bool atomic_compare_and_swap(ET* a, ET oldval, ET newval) {
    static_assert(sizeof(ET) <= 8, "Bad CAS length");
    if (sizeof(ET) == 1) {
      uint8_t r_oval, r_nval;
      memcpy(&r_oval, &oldval, sizeof(ET));
      memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(a), r_oval, r_nval);
    } else if (sizeof(ET) == 4) {
      uint32_t r_oval, r_nval;
      memcpy(&r_oval, &oldval, sizeof(ET));
      memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(a), r_oval, r_nval);
    } else { // if (sizeof(ET) == 8) {
      uint64_t r_oval, r_nval;
      memcpy(&r_oval, &oldval, sizeof(ET));
      memcpy(&r_nval, &newval, sizeof(ET));
      return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(a), r_oval, r_nval);
    }
  }

  template <typename E, typename EV>
  inline E fetch_and_add(E *a, EV b) {
    volatile E newV, oldV;
    do {oldV = *a; newV = oldV + b;}
    while (!atomic_compare_and_swap(a, oldV, newV));
    return oldV;
  }

  template <typename E, typename EV>
  inline void write_add(E *a, EV b) {
    //volatile E newV, oldV;
    E newV, oldV;
    do {oldV = *a; newV = oldV + b;}
    while (!atomic_compare_and_swap(a, oldV, newV));
  }

}  // namespace utils
}  // namespace cpam
