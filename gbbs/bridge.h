// This file is a bridge connecting the "lib interface" gbbs exports and the
// interfact that the current pbbslib exports. We would like to support both
// C++11 users, and the current (C++17) implementation of the lib. Using this
// bridge will hopefully simplify having two separate implementations of the lib
// interface.

#pragma once

#include <type_traits>
#include <utility>

#include "pbbslib/binary_search.h"
#include "pbbslib/counting_sort.h"
#include "pbbslib/integer_sort.h"
#include "pbbslib/monoid.h"
#include "pbbslib/parallel.h"
#include "pbbslib/random.h"
#include "pbbslib/random_shuffle.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/seq.h"
#include "pbbslib/sequence_ops.h"
#include "pbbslib/utilities.h"

// ================== parallel primitives ===================
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

//template <typename F>
//static void par_for(size_t start, size_t end, F f, bool parallel=true) {
//  size_t n = end - start;
//  size_t granularity = (n > 100) ? ceil(sqrt(n)) : 100;
//  par_for<F>(start, end, granularity, f, parallel);
//}

// Alias template so that sequence is exposed w/o namespacing
template<typename T>
using sequence = pbbs::sequence<T>;

// Alias template so that range is exposed w/o namespacing
template<typename T>
using range = pbbs::range<T>;


// Bridge to pbbslib (c++17)
namespace pbbslib {


  constexpr const size_t kSequentialForThreshold = 2048;

  // ====================== utilities =======================
  using empty = pbbs::empty;
  using flags = pbbs::flags;
  const flags no_flag = pbbs::no_flag;
  const flags fl_sequential = pbbs::fl_sequential;
  const flags fl_debug = pbbs::fl_debug;
  const flags fl_time = pbbs::fl_time;
  const flags fl_conservative = pbbs::fl_conservative;
  const flags fl_inplace = pbbs::fl_inplace;
  const flags fl_scan_inclusive = pbbs::fl_scan_inclusive;

  using pbbs::free_array;
  using pbbs::delete_array;
  using pbbs::new_array_no_init;
  using pbbs::new_array;
  using pbbs::hash32;
  using pbbs::hash32_2;
  using pbbs::hash32_3;
  using pbbs::hash64;
  using pbbs::hash64_2;
  using pbbs::atomic_compare_and_swap;
  using pbbs::fetch_and_add;
  using pbbs::write_add;
  using pbbs::write_max;
  using pbbs::write_min;
  using pbbs::log2_up;
  using pbbs::granularity;

  template<typename T>
  inline void assign_uninitialized(T& a, const T& b) {
    new (static_cast<void*>(std::addressof(a))) T(b);
  }

  template<typename T>
  inline void move_uninitialized(T& a, const T& b) {
    new (static_cast<void*>(std::addressof(a))) T(std::move(b));
  }

  // Currently unused; including commented out.
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
    return pbbs::write_min<ET>(a, b, std::less<ET>());
  }

  template <typename ET>
  inline bool write_max(ET *a, ET b) {
    return pbbs::write_max<ET>(a, b, std::less<ET>());
  }

  // ========================= monoid ==========================

  using pbbs::make_monoid;

  template <class T>
  using minm = pbbs::minm<T>;

  template <class T>
  using maxm = pbbs::maxm<T>;

  template <class T>
  using addm = pbbs::addm<T>;


  // ====================== sequence ops =======================

  using pbbs::scan_inplace;
  using pbbs::scan;
  using pbbs::reduce;
  using pbbs::pack;
  using pbbs::pack_index;
  using pbbs::pack_out;
  using pbbs::filter;
  using pbbs::filter_out;

  // used so second template argument can be inferred
  template <class T, class F>
  inline pbbs::delayed_sequence<T,F> make_sequence (size_t n, F f) {
    return pbbs::delayed_sequence<T,F>(n,f);
  }

  // TODO: call this make_range. make_sequence is bogus.
  template <class T>
  inline pbbs::range<T*> make_sequence (T* A, size_t n) {
    return pbbs::range<T*>(A, A+n);
  }

  // Scans the input sequence using the addm monoid.
  //
  // This computes in-place an exclusive prefix sum on the input sequence, that is,
  //   Out[i] = In[0] + In[1] + ... + In[i - 1].
  // The return value is the sum over the whole input sequence.
  template <RANGE In_Seq>
  inline auto scan_add_inplace(
      In_Seq&& In,
      flags fl = no_flag,
      typename std::remove_reference<In_Seq>::type::value_type* tmp = nullptr)
    -> typename std::remove_reference<In_Seq>::type::value_type {
    using T = typename std::remove_reference<In_Seq>::type::value_type;
    return pbbs::scan_inplace(
        std::forward<In_Seq>(In), pbbs::addm<T>(), fl, tmp);
  }

  template <class Seq>
  inline auto reduce_add(Seq const& I, flags fl = no_flag) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return pbbs::reduce(I, pbbs::addm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_max(Seq const& I, flags fl = no_flag) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return pbbs::reduce(I, pbbs::maxm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_min(Seq const& I, flags fl = no_flag) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return pbbs::reduce(I, pbbs::minm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_xor(Seq const& I, flags fl = no_flag) -> typename Seq::value_type {
    using T = typename Seq::value_type;
    return pbbs::reduce(I, pbbs::xorm<T>(), fl);
  }

  // Writes the list of indices `i` where `Fl[i] == true` to range `Out`.
  template <SEQ Bool_Seq, RANGE Out_Seq>
  size_t pack_index_out(Bool_Seq const &Fl, Out_Seq&& Out,
                flags fl = no_flag) {
    using Idx_Type = typename std::remove_reference<Out_Seq>::type::value_type;
    auto identity = [] (size_t i) {return (Idx_Type) i;};
    return pbbs::pack_out(
        pbbs::delayed_seq<Idx_Type>(Fl.size(),identity),
        Fl,
        std::forward<Out_Seq>(Out),
        fl);
  }

  // ====================== binary search =======================

  using pbbs::binary_search;

  // ====================== sample sort =======================

  using pbbs::sample_sort;
  using pbbs::sample_sort_inplace;

  // ====================== integer sort =======================

  using pbbs::integer_sort_inplace;
  using pbbs::integer_sort;
  using pbbs::count_sort;

  // ====================== random shuffle =======================
  using random = pbbs::random;
  using pbbs::random_permutation;
  using pbbs::random_shuffle;

}


// Other extensions to pbbs used by the graph benchmarks.
namespace pbbslib {

  constexpr size_t _F_BSIZE = 2000;

  template <class T>
  void free_arrays(T* first) {
    pbbs::free_array(first);
  }

  template <class T, typename... Args>
  void free_arrays(T* first, Args... args) {
    pbbs::free_array(first);
    free_arrays(args...);
  }

  template <class Idx_Type, class D, class F>
  inline pbbs::sequence<std::tuple<Idx_Type, D> > pack_index_and_data(
      F& f, size_t size, flags fl = no_flag) {
    auto id_seq = pbbslib::make_sequence<std::tuple<Idx_Type, D> >(size,  [&](size_t i) {
      return std::make_tuple((Idx_Type)i, std::get<1>(f[i]));
    });
    auto flgs_seq = pbbslib::make_sequence<bool>(size, [&](size_t i) { return std::get<0>(f[i]); });

    return pbbs::pack(id_seq, flgs_seq, fl);
  }

  template <class T, class Pred>
  inline size_t filter_seq(T* in, T* out, size_t n, Pred p) {
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
      if (p(in[i])) out[k++] = in[i];
    return k;
  }

  // Faster for a small number in output (about 40% or less)
  // Destroys the input.   Does not need a bool array.
  template <class T, class PRED>
  inline size_t filterf(T* In, T* Out, size_t n, PRED p) {
    size_t b = _F_BSIZE;
    if (n < b) return filter_seq(In, Out, n, p);
    size_t l = pbbs::num_blocks(n, b);
    size_t* Sums = new_array_no_init<size_t>(l + 1);
    par_for(0, l, 1, [&] (size_t i) {
      size_t s = i * b;
      size_t e = std::min(s + b, n);
      size_t k = s;
      for (size_t j = s; j < e; j++) {
        if (p(In[j])) In[k++] = In[j];
      }
      Sums[i] = k - s;
    });
    auto isums = make_sequence(Sums, l);
    size_t m = scan_add_inplace(isums.slice());
    Sums[l] = m;
    par_for(0, l, 1, [&] (size_t i) {
      T* I = In + i * b;
      T* O = Out + Sums[i];
      for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
        O[j] = I[j];
      }
    });
    pbbslib::free_array(Sums);
    return m;
  }


  // Faster for a small number in output (about 40% or less)
  // Destroys the input.   Does not need a bool array.
  template <class T, class PRED, class OUT>
  inline size_t filterf(T* In, size_t n, PRED p, OUT out, size_t out_off) {
    size_t b = _F_BSIZE;
    if (n < b) {
      size_t k = out_off;
      for (size_t i = 0; i < n; i++) {
        if (p(In[i])) out(k++, In[i]);
      }
      return k - out_off;
    }
    size_t l = pbbs::num_blocks(n, b);
    size_t* Sums = new_array_no_init<size_t>(l + 1);
    par_for(0, l, 1, [&] (size_t i) {
      size_t s = i * b;
      size_t e = std::min(s + b, n);
      size_t k = s;
      for (size_t j = s; j < e; j++) {
        if (p(In[j])) In[k++] = In[j];
      }
      Sums[i] = k - s;
    });
    auto isums = make_sequence(Sums, l);
    size_t m = scan_add_inplace(isums.slice());
    Sums[l] = m;
    par_for(0, l, 1, [&] (size_t i) {
      T* I = In + i * b;
      size_t si = out_off + Sums[i];
      for (size_t j = 0; j < Sums[i + 1] - Sums[i]; j++) {
        out(si + j, I[j]);
      }
    });
    pbbslib::free_array(Sums);
    return m;
  }

  template <class T, class PRED>
  inline size_t filterf_and_clear(T* In, T* Out, size_t n, PRED p, T& empty) {
    size_t b = _F_BSIZE;
    if (n < b) {
      size_t ret = filter_seq(In, Out, n, p);
      for (size_t i=0; i<n; i++) {
        if (p(In[i])) {
          In[i] = empty;
        }
      }
      return ret;
    }
    size_t l = pbbs::num_blocks(n, b);
    b = pbbs::num_blocks(n, l);
    size_t* Sums = new_array_no_init<size_t>(l + 1);

    par_for(0, l, 1, [&] (size_t i) {
      size_t s = i * b;
      size_t e = std::min(s + b, n);
      size_t k = s;
      for (size_t j = s; j < e; j++) {
        if (p(In[j])) {
          In[k] = In[j];
          if (k != j) {
            In[j] = empty;
          }
          k++;
        }
      }
      Sums[i] = k - s;
    });
    auto isums = make_sequence(Sums, l);
    size_t m = scan_add_inplace(isums.slice());
    Sums[l] = m;
    par_for(0, l, 1, [&] (size_t i) {
      T* I = In + (i * b);
      size_t i_off = Sums[i];
      size_t num_i = Sums[i+1] - i_off;
      T* O = Out + i_off;
      for (size_t j = 0; j < num_i; j++) {
        O[j] = I[j];
        I[j] = empty;
      }
    });
    pbbslib::free_array(Sums);
    return m;
  }

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
  sequence<char> sequence_to_string(TSeq const &T) {
    size_t n = T.size();
    auto S = sequence<size_t>(n, [&] (size_t i) {
      return t_to_stringlen(T[i])+1; // +1 for \n
    });
    size_t m = pbbslib::scan_inplace(S.slice(), addm<size_t>());

    auto C = sequence<char>(m, [&] (size_t i) { return (char)0; });
    parallel_for(0, n-1, [&] (size_t i) {
      type_to_string(C.begin() + S[i], T[i]);
      C[S[i + 1] - 1] = '\n';
    });
    type_to_string(C.begin() + S[n - 1], T[n - 1]);
    C[m - 1] = '\n';

    return pbbslib::filter(C, [&] (char A) { return A > 0; });
  }

}
