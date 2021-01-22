// This file is a bridge connecting the "lib interface" gbbs exports and the
// interfact that the current pbbslib exports. We would like to support both
// C++11 users, and the current (C++17) implementation of the lib. Using this
// bridge will hopefully simplify having two separate implementations of the lib
// interface.

#pragma once

#include <type_traits>
#include <utility>

#include "parlay/internal/binary_search.h"
#include "parlay/primitives.h"
#include "parlay/monoid.h"
#include "parlay/parallel.h"
#include "parlay/random.h"
#include "parlay/delayed_sequence.h"
#include "parlay/sequence.h"
#include "parlay/slice.h"
#include "parlay/range.h"
#include "parlay/utilities.h"

namespace gbbs {
  // ================== parallel primitives ===================

  struct empty {};

  constexpr const size_t kSequentialForThreshold = 2048;

  using parlay::parallel_for;
  using parlay::par_do;
  // using pbbs::parallel_for_alloc;  // TODO
  using parlay::num_workers;
  using parlay::worker_id;

  template <class T>
  using slice = parlay::slice<T*, T*>;

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

  template <class E, class T>
  slice<E> make_slice(T* start, T* end) {
    return parlay::make_slice((E*)start, (E*)end);
  }

  // TODO: deprecate
  template <typename F>
  static void par_for(size_t start, size_t end, size_t granularity, F f, bool parallel=true) {
    if (!parallel) {
      for (size_t i=start; i<end; i++) {
        f(i);
      }
    } else {
      parallel_for(start, end, f, granularity);
    }
  }

  // TODO: deprecate
  template <typename F>
  static void par_for(size_t start, size_t end, F f, bool parallel=true, size_t granularity=std::numeric_limits<size_t>::max()) {
    if (!parallel) {
      for (size_t i=start; i<end; i++) {
        f(i);
      }
    } else {
      if (granularity == std::numeric_limits<size_t>::max()) {
        parallel_for(start, end, f);
      } else {
        parallel_for(start, end, f, granularity);
      }
    }
  }

  // Alias template so that sequence is exposed w/o namespacing
  template<typename T>
  using sequence = parlay::sequence<T>;

  inline size_t num_blocks(size_t n, size_t block_size) {
    if (n == 0)
      return 0;
    else
      return (1 + ((n)-1) / (block_size));
  }

  // Combines two hash values.
  inline uint64_t hash_combine(uint64_t hash_value_1, uint64_t hash_value_2) {
    // This is the same as boost's 32-bit `hash_combine` implementation, but with
    // 2 ^ 64 / (golden ratio) chosen as an arbitrary 64-bit additive magic number
    // rather than 2 ^ 32 / (golden ratio).
    return hash_value_1 ^ (hash_value_2 + 0x9e3779b97f4a7c15 + (hash_value_1 << 6)
        + (hash_value_1 >> 2));
  }

}  // namespace gbbs

// Bridge to pbbslib (c++17)
namespace pbbslib {

  // ====================== utilities =======================
  using empty = gbbs::empty;  // TODO: handle defining

  using parlay::parallel_for;
  using parlay::par_do;
  // using parlay::parallel_for_alloc;  // TODO
  using parlay::num_workers;
  using parlay::worker_id;
  using parlay::log2_up;

  using parlay::hash32;
  using parlay::hash32_2;
  using parlay::hash32_3;
  using parlay::hash64;
  using parlay::hash64_2;


template <typename ET>
inline bool atomic_compare_and_swap(ET* a, ET oldval, ET newval) {
  if constexpr (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 16) {
    __int128 r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap_16(reinterpret_cast<__int128*>(a),
                                           r_oval, r_nval);
  } else {
    std::cout << "Bad CAS Length" << sizeof(ET) << std::endl;
    exit(0);
  }
}

template <typename ET>
inline bool atomic_compare_and_swap(volatile ET* a, ET oldval, ET newval) {
  if (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint8_t*>(a),
                                        r_oval, r_nval);
  } else if (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint32_t*>(a),
                                        r_oval, r_nval);
  } else if (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint64_t*>(a),
                                        r_oval, r_nval);
  } else if (sizeof(ET) == 16) {
    __int128 r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap_16(
        reinterpret_cast<volatile __int128*>(a), r_oval, r_nval);
  } else {
    std::cout << "Bad CAS Length" << sizeof(ET) << std::endl;
    exit(0);
  }
}

template <typename E, typename EV>
inline E fetch_and_add(E* a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
  return oldV;
}

template <typename E, typename EV>
inline void write_add(E* a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
}

template <typename E, typename EV>
inline void write_add(std::atomic<E>* a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = a->load();
    newV = oldV + b;
  } while (!std::atomic_compare_exchange_strong(a, &oldV, newV));
}

template <typename ET, typename F>
inline bool write_min(ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_min(volatile ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_min(std::atomic<ET>* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = a->load();
  while (less(b, c) && !(r = std::atomic_compare_exchange_strong(a, &c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(volatile ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(std::atomic<ET>* a, ET b, F less) {
  ET c;
  bool r = 0;
  do
    c = a->load();
  while (less(c, b) && !(r = std::atomic_compare_exchange_strong(a, &c, b)));
  return r;
}


  // Alias template so that sequence is exposed w/o namespacing
  template<typename T>
  using sequence = parlay::sequence<T>;

  template<typename T>
  inline void assign_uninitialized(T& a, const T& b) {
    new (static_cast<void*>(std::addressof(a))) T(b);
  }

  template<typename T>
  inline void move_uninitialized(T& a, const T& b) {
    new (static_cast<void*>(std::addressof(a))) T(std::move(b));
  }

  // Currently unused, but may be useful in the future; including commented out.
  // template <class ET>
  // inline bool CAS128(ET* a, ET b, ET c) {
  //   return __sync_bool_compare_and_swap_16((__int128*)a, *((__int128*)&b),
  //                                          *((__int128*)&c));
  // }

  template <typename ET>
  inline bool CAS(ET* ptr, const ET oldv, const ET newv) {
    return atomic_compare_and_swap(ptr, oldv, newv);
  }

  inline long xaddl(long* variable, long value) {
    asm volatile("lock; xaddl %%eax, %2;"
                 : "=a"(value)                 // Output
                 : "a"(value), "m"(*variable)  // Input
                 : "memory");
    return value;
  }

  inline int xaddi(int* variable, int value) {
    asm volatile("lock; xadd %%eax, %2;"
                 : "=a"(value)                 // Output
                 : "a"(value), "m"(*variable)  // Input
                 : "memory");
    return value;
  }

  // The conditional should be removed by the compiler
  // this should work with pointer types, or pairs of integers
  template <class ET>
  inline ET xadd(ET* variable, ET value) {
    if (sizeof(ET) == 8) {
      return xaddl((long*)variable, (long)value);
    } else if (sizeof(ET) == 4) {
      return xaddi((int*)variable, (int)value);
    } else {
      std::cout << "xadd bad length"
                << "\n";
      abort();
    }
  }

  template <typename ET>
  inline bool write_min(ET *a, ET b) {
    return write_min<ET>(a, b, std::less<ET>());
  }

  template <typename ET>
  inline bool write_max(ET *a, ET b) {
    return write_max<ET>(a, b, std::less<ET>());
  }


  // Scans the input sequence using the addm monoid.
  //
  // This computes in-place an exclusive prefix sum on the input sequence, that is,
  //   Out[i] = In[0] + In[1] + ... + In[i - 1].
  // The return value is the sum over the whole input sequence.
  template <class In_Seq>
  inline auto scan_add_inplace(
      In_Seq&& In)
    -> typename std::remove_reference<In_Seq>::type::value_type {
    using T = typename std::remove_reference<In_Seq>::type::value_type;
    return parlay::scan_inplace(std::forward<In_Seq>(In), parlay::addm<T>());
  }

  template <class Seq>
  inline auto reduce_add(Seq const& I) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return parlay::reduce(parlay::make_slice(I), parlay::addm<T>());
  }

  template <class Seq>
  inline auto reduce_max(Seq const& I) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return parlay::reduce(parlay::make_slice(I), parlay::maxm<T>());
  }

  template <class Seq>
  inline auto reduce_min(Seq const& I) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return parlay::reduce(I, parlay::minm<T>());
  }

  template <class Seq>
  inline auto reduce_xor(Seq const& I) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return parlay::reduce(I, parlay::xorm<T>());
  }

  // Writes the list of indices `i` where `Fl[i] == true` to range `Out`.
  template <class Bool_Seq, class Out_Seq>
  size_t pack_index_out(Bool_Seq const &Fl, Out_Seq&& Out) {
    using Idx_Type = typename std::remove_reference<Out_Seq>::type::value_type;
    auto identity = [] (size_t i) {return (Idx_Type) i;};
    return parlay::pack_into(
        parlay::delayed_seq<Idx_Type>(Fl.size(),identity),
        Fl,
        std::forward<Out_Seq>(Out));
  }


  // ====================== binary search =======================

  using parlay::internal::binary_search;

  // ====================== sample sort =======================

  using parlay::internal::sample_sort;
  using parlay::internal::sample_sort_inplace;

  // ====================== integer sort =======================

  using parlay::integer_sort_inplace;
  using parlay::integer_sort;
  using parlay::internal::count_sort;

  // ====================== random shuffle =======================
  using random = parlay::random;
}


// Other extensions to pbbs used by the graph benchmarks.
namespace pbbslib {

  constexpr size_t _F_BSIZE = 2000;


  template <class Idx_Type, class D, class F>
  inline parlay::sequence<std::tuple<Idx_Type, D> > pack_index_and_data(
      F& f, size_t size) {
    auto id_seq = parlay::delayed_seq<std::tuple<Idx_Type, D> >(size,  [&](size_t i) {
      return std::make_tuple((Idx_Type)i, std::get<1>(f[i]));
    });
    auto flgs_seq = parlay::delayed_seq<bool>(size, [&](size_t i) { return std::get<0>(f[i]); });

    return parlay::pack(id_seq, flgs_seq);
  }

//  template <class T, class Pred>
//  inline size_t filter_seq(T* in, T* out, size_t n, Pred p) {
//    size_t k = 0;
//    for (size_t i = 0; i < n; i++)
//      if (p(in[i])) out[k++] = in[i];
//    return k;
//  }
//
//  // Faster for a small number in output (about 40% or less)
//  // Destroys the input.   Does not need a bool array.
//  template <class T, class PRED>
//  inline size_t filterf(T* In, T* Out, size_t n, PRED p) {
//    size_t b = _F_BSIZE;
//    if (n < b) return filter_seq(In, Out, n, p);
//    size_t l = pbbs::num_blocks(n, b);
//    size_t* Sums = new_array_no_init<size_t>(l + 1);
//    parallel_for(0, l, [&] (size_t i) {
//      size_t s = i * b;
//      size_t e = std::min(s + b, n);
//      size_t k = s;
//      for (size_t j = s; j < e; j++) {
//        if (p(In[j])) In[k++] = In[j];
//      }
//      Sums[i] = k - s;
//    }, 1);
//    auto isums = make_sequence(Sums, l);
//    size_t m = scan_add_inplace(isums.slice());
//    Sums[l] = m;
//    parallel_for(0, l, [&] (size_t i) {
//      T* I = In + i * b;
//      T* O = Out + Sums[i];
//      for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
//        O[j] = I[j];
//      }
//    }, 1);
//    pbbslib::free_array(Sums);
//    return m;
//  }
//
//
//  // Faster for a small number in output (about 40% or less)
//  // Destroys the input.   Does not need a bool array.
//  template <class T, class PRED, class OUT>
//  inline size_t filterf(T* In, size_t n, PRED p, OUT out, size_t out_off) {
//    size_t b = _F_BSIZE;
//    if (n < b) {
//      size_t k = out_off;
//      for (size_t i = 0; i < n; i++) {
//        if (p(In[i])) out(k++, In[i]);
//      }
//      return k - out_off;
//    }
//    size_t l = pbbs::num_blocks(n, b);
//    size_t* Sums = new_array_no_init<size_t>(l + 1);
//    parallel_for(0, l, [&] (size_t i) {
//      size_t s = i * b;
//      size_t e = std::min(s + b, n);
//      size_t k = s;
//      for (size_t j = s; j < e; j++) {
//        if (p(In[j])) In[k++] = In[j];
//      }
//      Sums[i] = k - s;
//    }, 1);
//    auto isums = make_sequence(Sums, l);
//    size_t m = scan_add_inplace(isums.slice());
//    Sums[l] = m;
//    parallel_for(0, l, [&] (size_t i) {
//      T* I = In + i * b;
//      size_t si = out_off + Sums[i];
//      for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
//        out(si + j, I[j]);
//      }
//    }, 1);
//    pbbslib::free_array(Sums);
//    return m;
//  }
//
//  template <class T, class PRED>
//  inline size_t filterf_and_clear(T* In, T* Out, size_t n, PRED p, T& empty) {
//    size_t b = _F_BSIZE;
//    if (n < b) {
//      size_t ret = filter_seq(In, Out, n, p);
//      for (size_t i=0; i<n; i++) {
//        if (p(In[i])) {
//          In[i] = empty;
//        }
//      }
//      return ret;
//    }
//    size_t l = pbbs::num_blocks(n, b);
//    b = pbbs::num_blocks(n, l);
//    size_t* Sums = new_array_no_init<size_t>(l + 1);
//
//    parallel_for(0, l, [&] (size_t i) {
//      size_t s = i * b;
//      size_t e = std::min(s + b, n);
//      size_t k = s;
//      for (size_t j = s; j < e; j++) {
//        if (p(In[j])) {
//          In[k] = In[j];
//          if (k != j) {
//            In[j] = empty;
//          }
//          k++;
//        }
//      }
//      Sums[i] = k - s;
//    }, 1);
//    auto isums = make_sequence(Sums, l);
//    size_t m = scan_add_inplace(isums.slice());
//    Sums[l] = m;
//    parallel_for(0, l, [&] (size_t i) {
//      T* I = In + (i * b);
//      size_t i_off = Sums[i];
//      size_t num_i = Sums[i+1] - i_off;
//      T* O = Out + i_off;
//      for (size_t j = 0; j < num_i; j++) {
//        O[j] = I[j];
//        I[j] = empty;
//      }
//    }, 1);
//    pbbslib::free_array(Sums);
//    return m;
//  }

  template <class E, class I, class P>
  struct filter_iter {
    I& iter;
    P& pred;
    E cur_val;

    filter_iter(I& _it, P& _pr) : iter(_it), pred(_pr) {
      cur_val = iter.cur();
      while (!pred(cur_val) && iter.has_next()) {
        cur_val = iter.next();
      }
    }

    E cur() { return cur_val; }

    E next() {
      while (iter.has_next()) {
        cur_val = iter.next();
        if (pred(cur_val)) {
          break;
        }
      }
      return cur_val;
    }

    // has_next
  };

  template <class E, class I, class P>
  inline filter_iter<E, I, P> make_filter_iter(I& _it, P& _pr) {
    return filter_iter<E, I, P>(_it, _pr);
  }

  int t_to_stringlen(long a);
  void type_to_string(char* s, long a);

  int t_to_stringlen(unsigned long a);
  void type_to_string(char* s, unsigned long a);

  uint t_to_stringlen(uint a);
  void type_to_string(char* s, uint a);

  int t_to_stringlen(int a);
  void type_to_string(char* s, int a);

  int t_to_stringlen(double a);

  int t_to_stringlen(char* a);
  void type_to_string(char* s, char* a);

  void type_to_string(char* s, double a);

  template <class A, class B>
  inline int t_to_stringlen(std::pair<A, B> a) {
    return t_to_stringlen(a.first) + t_to_stringlen(a.second) + 1;
  }

  template <class A, class B>
  inline int t_to_stringlen(std::tuple<A, B> a) {
    return t_to_stringlen(std::get<0>(a)) + t_to_stringlen(std::get<1>(a)) + 1;
  }

  template <class A, class B, class C>
  inline int t_to_stringlen(std::tuple<A, B, C> a) {
    return t_to_stringlen(std::get<0>(a)) + t_to_stringlen(std::get<1>(a)) + t_to_stringlen(std::get<2>(a)) + 2;
  }

  template <class A, class B>
  inline void type_to_string(char* s, std::pair<A, B> a) {
    int l = t_to_stringlen(a.first);
    type_to_string(s, a.first);
    s[l] = ' ';
    type_to_string(s + l + 1, a.second);
  }

  template <class A, class B>
  inline void type_to_string(char* s, std::tuple<A, B> a) {
    int l = t_to_stringlen(std::get<0>(a));
    type_to_string(s, std::get<0>(a));
    s[l] = ' ';
    type_to_string(s + l + 1, std::get<1>(a));
  }

  template <class A, class B, class C>
  inline void type_to_string(char* s, std::tuple<A, B, C> a) {
    int l = t_to_stringlen(std::get<0>(a));
    type_to_string(s, std::get<0>(a));
    s[l] = ' ';
    int l1 = t_to_stringlen(std::get<1>(a));
    type_to_string(s + l + 1, std::get<1>(a));
    s[l + l1 + 1] = ' ';
    type_to_string(s + l + l1 + 2, std::get<2>(a));
  }

  template <class TSeq>
  parlay::sequence<char> sequence_to_string(TSeq const &T) {
    size_t n = T.size();
    auto S = parlay::sequence<size_t>(n, [&] (size_t i) {
      return t_to_stringlen(T[i])+1; // +1 for \n
    });
    size_t m = parlay::scan_inplace(parlay::make_slice(S), parlay::addm<size_t>());

    auto C = parlay::sequence<char>(m, [&] (size_t i) { return (char)0; });
    parallel_for(0, n-1, [&] (size_t i) {
      type_to_string(C.begin() + S[i], T[i]);
      C[S[i + 1] - 1] = '\n';
    });
    type_to_string(C.begin() + S[n - 1], T[n - 1]);
    C[m - 1] = '\n';

    return parlay::filter(parlay::make_slice(C), [&] (char A) { return A > 0; });
  }
}
