// This file is a bridge connecting the "lib interface" gbbs exports and the
// interfact that the current pbbslib exports. We would like to support both
// C++11 users, and the current (C++17) implementation of the lib. Using this
// bridge will hopefully simplify having two separate implementations of the lib
// interface.

#pragma once

#include "pbbslib/binary_search.h"
#include "pbbslib/monoid.h"
#include "pbbslib/parallel.h"
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
static void par_for(size_t start, size_t end, F f, bool parallel=true) {
  size_t n = end - start;
  size_t granularity = (n > 100) ? ceil(sqrt(n)) : 100;
  par_for<F>(start, end, granularity, f, parallel);
}

// Bridge to pbbslib (c++17)
namespace pbbslib {

  // Open pbbs namespace internally.
  using namespace pbbs;

  constexpr const size_t kSequentialForThreshold = 4000;

  // ====================== utilities =======================
  using empty = pbbs::empty;
  using flags = pbbs::flags;
  const flags no_flag = pbbs::no_flag;
  const flags fl_sequential = pbbs::fl_sequential;
  const flags fl_debug = pbbs::fl_debug;
  const flags fl_time = pbbs::fl_time;
  const flags fl_conservative = pbbs::fl_conservative;
  const flags fl_inplace = pbbs::fl_inplace;

  inline void free_array(void* a) {
    return pbbs::free_array(a);
  }

  // Destructs in parallel
  template<typename E>
  void delete_array(E* A, size_t n) {
    return pbbs::delete_array<E>(A, n);
  }

  // Does not initialize the array
  template<typename E>
  E* new_array_no_init(size_t n, bool touch_pages=false) {
    return pbbs::new_array_no_init<E>(n, touch_pages);
  }

  // Initializes in parallel
  template<typename E>
  E* new_array(size_t n) {
    return pbbs::new_array<E>(n);
  }

  // a 32-bit hash function
  inline uint32_t hash32(uint32_t a) {
    return pbbs::hash32(a);
  }

  // cheaper, likely more sketchy versions
  inline uint32_t hash32_2(uint32_t a) {
    return pbbs::hash32_2(a);
  }

  // cheaper, likely more sketchy versions
  inline uint32_t hash32_3(uint32_t a) {
    return pbbs::hash32_3(a);
  }

  // from numerical recipes
  inline uint64_t hash64(uint64_t u )
  {
    return pbbs::hash64(u);
  }

  // a slightly cheaper, but possibly not as good version
  // based on splitmix64
  inline uint64_t hash64_2(uint64_t x) {
    return pbbs::hash64_2(x);
  }

  template <typename ET>
  inline bool CAS_GCC(ET* ptr, const ET oldv, const ET newv) {
    return pbbs::CAS_GCC(ptr, oldv, newv);
  }

  // Currently unused; including commented out.
  // template <class ET>
  // inline bool CAS128(ET* a, ET b, ET c) {
  //   return __sync_bool_compare_and_swap_16((__int128*)a, *((__int128*)&b),
  //                                          *((__int128*)&c));
  // }

  template <typename ET>
  inline bool atomic_compare_and_swap(ET* ptr, ET oldv, ET newv) {
    return pbbs::atomic_compare_and_swap(ptr, oldv, newv);
  }

  inline bool atomic_compare_and_swap(double* a, const double &oldval, const double &newval) {
    return pbbs::atomic_compare_and_swap(a, oldval, newval);
  };

  inline bool atomic_compare_and_swap(float* a, float &oldval, float &newval) {
    return pbbs::atomic_compare_and_swap(a, oldval, newval);
  };

  template <typename E, typename EV>
  inline E fetch_and_add(E *a, EV b) {
    return pbbs::fetch_and_add<E, EV>(a, b);
  }

  template <typename E, typename EV>
  inline void write_add(E *a, EV b) {
    pbbs::write_add<E, EV>(a, b);
  }

  template <typename E, typename EV>
  inline void write_add(std::atomic<E> *a, EV b) {
    pbbs::write_add<E, EV>(a, b);
  }

  template <typename ET, typename F>
  inline bool write_min(ET *a, ET b, F less) {
    return pbbs::write_min<ET, F>(a, b, less);
  }

  template <typename ET, typename F>
  inline bool write_min(std::atomic<ET> *a, ET b, F less) {
    return pbbs::write_min<ET, F>(a, b, less);
  }

  // returns the log base 2 rounded up (works on ints or longs or unsigned versions)
  template <class T>
  inline size_t log2_up(T i) {
    size_t a=0;
    T b=i-1;
    while (b > 0) {b = b >> 1; a++;}
    return a;
  }

  inline size_t granularity(size_t n) {
    return (n > 100) ? ceil(pow(n,0.5)) : 100;
  }

  inline void assert_str(int cond, std::string s) {
    return pbbs::assert_str(cond, s);
  }


  // ====================== sequence ops =======================
  // used so second template argument can be inferred
  template <class T, class F>
  inline delayed_sequence<T,F> make_sequence (size_t n, F f) {
    return delayed_sequence<T,F>(n,f);
  }

  template <class T>
  inline range<T*> make_sequence (T* A, size_t n) {
    return range<T*>(A, A+n);
  }

  template <RANGE Range, class Monoid>
  inline auto scan_inplace(Range In, Monoid m, flags fl = no_flag)
    -> typename Range::value_type {
    return pbbs::scan_inplace<Range, Monoid>(In, m, fl);
  }

  template <SEQ In_Seq, class Monoid>
  inline auto scan(In_Seq const &In, Monoid m, flags fl = no_flag)
    ->  std::pair<sequence<typename In_Seq::value_type>, typename In_Seq::value_type>
  {
    return pbbs::scan<In_Seq, Monoid>(In, m, fl);
  }

  // do in place if rvalue reference to a sequence<T>
  template <class T, class Monoid>
  auto scan(sequence<T> &&In, Monoid m, flags fl = no_flag)
    ->  std::pair<sequence<T>, T> {
    return pbbs::scan(std::move(In), m, fl);
  }

  // Scans the input sequence using the addm monoid.
  template <class In_Seq>
  inline auto scan_add_inplace(In_Seq const& In, flags fl = no_flag) -> typename In_Seq::value_type {
    using T = typename In_Seq::T;
    return scan_inplace(In, pbbslib::addm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_add(Seq const& I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    return reduce(I, addm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_max(Seq const& I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    return reduce(I, maxm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_min(Seq const& I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    return reduce(I, minm<T>(), fl);
  }

  template <class Seq>
  inline auto reduce_xor(Seq const& I, flags fl = no_flag) -> typename Seq::T {
    using T = typename Seq::T;
    return reduce(I, xorm<T>(), fl);
  }


  // ====================== binary search =======================
  // return index to first key greater or equal to v
  template <typename Seq, typename F>
  inline size_t binary_search(Seq const &I, typename Seq::value_type v,
		       const F& less) {
    return pbbs::binary_search<Seq, F>(I, v, less);
  }

  // return index to first key where less is false
  template <typename Seq, typename F>
  inline size_t binary_search(Seq const &I, const F& less) {
    return pbbs::binary_search<Seq, F>(I, less);
  }


  // ====================== sample sort =======================
  template<class Seq, typename BinPred>
  auto sample_sort (Seq const &A, const BinPred& f, bool stable = false)
    -> sequence<typename Seq::value_type> {
      return pbbs::sample_sort(A, f, stable);
  }

  template<class Iter, typename BinPred>
  void sample_sort_inplace (range<Iter> A, const BinPred& f, bool stable = false) {
    return pbbs::sample_sort_inplace(A, f, stable);
  }

  template<typename E, typename BinPred, typename s_size_t>
  void sample_sort (E* A, s_size_t n, const BinPred& f, bool stable) {
    return pbbs::sample_sort(A, n, f, stable);
  }

  // ====================== random shuffle =======================
  template <class intT>
  sequence<intT> random_permutation(size_t n, pbbs::random r = default_random) {
    return pbbs::random_permutation<intT>(n, r);
  }
}




// Other extensions to pbbs used by the graph benchmarks.
namespace pbbslib {

  constexpr size_t _F_BSIZE = 2000;

  template <class Idx_Type, class D, class F>
  inline sequence<std::tuple<Idx_Type, D> > pack_index_and_data(
      F& f, size_t size, flags fl = no_flag) {
    auto identity = [&](size_t i) {
      return std::make_tuple((Idx_Type)i, std::get<1>(f(i)));
    };
    auto flgs_f = [&](size_t i) { return std::get<0>(f(i)); };
    auto flgs_in =
        pbbslib::make_sequence<bool>(size, flgs_f);
    return pack(pbbslib::make_sequence<std::tuple<Idx_Type, D> >(size, identity), flgs_in,
                fl);
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
    size_t l = num_blocks(n, b);
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
    auto isums = sequence<size_t>(Sums, l);
    size_t m = scan_add_inplace(isums);
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
    size_t l = num_blocks(n, b);
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
    auto isums = sequence<size_t>(Sums, l);
    size_t m = scan_add_inplace(isums);
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
    size_t l = num_blocks(n, b);
    b = num_blocks(n, l);
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
    auto isums = sequence<size_t>(Sums, l);
    size_t m = scan_add_inplace(isums);
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

}
