// This file is a bridge connecting the "lib interface" gbbs exports and the
// interfact that the current pbbslib exports. We would like to support both
// C++11 users, and the current (C++17) implementation of the lib. Using this
// bridge will hopefully simplify having two separate implementations of the lib
// interface.

#pragma once

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

  inline void free_array(void* a) {
    return pbbs::free_array(a);
  }

  // Destructs in parallel
  template<typename E>
  void delete_array(E* A, size_t n) {
    return pbbs::delete_array<E>(A, n);
  }

  template<typename T>
  inline void assign_uninitialized(T& a, const T& b) {
    new (static_cast<void*>(std::addressof(a))) T(b);
  }

  template<typename T>
  inline void move_uninitialized(T& a, const T& b) {
    new (static_cast<void*>(std::addressof(a))) T(std::move(b));
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
    return pbbs::atomic_compare_and_swap(ptr, oldv, newv);
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

  template <typename ET>
  inline bool CAS(ET* ptr, const ET oldv, const ET newv) {
    return atomic_compare_and_swap(ptr, oldv, newv);
  }

  template <typename E, typename EV>
  inline E fetch_and_add(E *a, EV b) {
    return pbbs::fetch_and_add<E, EV>(a, b);
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

  template <typename E, typename EV>
  inline void write_add(E *a, EV b) {
    pbbs::write_add<E, EV>(a, b);
  }

  template <typename E, typename EV>
  inline void write_add(std::atomic<E> *a, EV b) {
    pbbs::write_add<E, EV>(a, b);
  }

  template <typename ET>
  inline bool write_min(ET *a, ET b) {
    return pbbs::write_min<ET>(a, b, std::less<ET>());
  }

  template <typename ET, typename F>
  inline bool write_min(ET *a, ET b, F less) {
    return pbbs::write_min<ET, F>(a, b, less);
  }

  template <typename ET>
  inline bool write_max(ET *a, ET b) {
    return pbbs::write_max<ET>(a, b, std::less<ET>());
  }

  template <typename ET, typename F>
  inline bool write_max(ET *a, ET b, F less) {
    return pbbs::write_max<ET, F>(a, b, less);
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


  // ========================= monoid ==========================
  template <class F, class T>
  pbbs::monoid<F,T> make_monoid (F f, T id) {
    return pbbs::monoid<F,T>(f, id);
  }

  template <class T>
  using minm = pbbs::minm<T>;

  template <class T>
  using maxm = pbbs::maxm<T>;

  template <class T>
  using addm = pbbs::addm<T>;


  // ====================== sequence ops =======================
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

  template <RANGE Range, class Monoid>
  inline auto scan_inplace(Range In, Monoid m, flags fl = no_flag)
    -> typename Range::value_type {
    return pbbs::scan_inplace<Range, Monoid>(In, m, fl);
  }

  template <SEQ In_Seq, class Monoid>
  inline auto scan(In_Seq const &In, Monoid m, flags fl = no_flag)
    ->  std::pair<pbbs::sequence<typename In_Seq::value_type>, typename In_Seq::value_type>
  {
    return pbbs::scan<In_Seq, Monoid>(In, m, fl);
  }

  // do in place if rvalue reference to a sequence<T>
  template <class T, class Monoid>
  auto scan(pbbs::sequence<T> &&In, Monoid m, flags fl = no_flag)
    ->  std::pair<pbbs::sequence<T>, T> {
    return pbbs::scan(std::move(In), m, fl);
  }

  // Scans the input sequence using the addm monoid.
  template <RANGE In_Seq>
  inline auto scan_add_inplace(In_Seq const& In, flags fl = no_flag) -> typename In_Seq::value_type {
    using T = typename In_Seq::value_type;
    return pbbs::scan_inplace(In, pbbs::addm<T>(), fl);
  }

  // Scans the input sequence using the addm monoid.
  template <class T>
  inline auto scan_add_inplace(sequence<T> const& In, flags fl = no_flag) -> T {
    return pbbs::scan_inplace(In.slice(), pbbs::addm<T>(), fl);
  }


  template <SEQ Seq, class Monoid>
  auto reduce(Seq const &A, Monoid m, flags fl = no_flag)
    -> typename Seq::value_type {
      return pbbs::reduce(A, m, fl);
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

  template <SEQ In_Seq, SEQ Bool_Seq>
  auto pack(In_Seq const &In, Bool_Seq const &Fl, flags fl = no_flag)
      -> pbbs::sequence<typename In_Seq::value_type> {
    return pbbs::pack(In, Fl, fl);
  }

  // Pack the output to the output range.
  template <SEQ In_Seq, SEQ Bool_Seq, RANGE Out_Seq>
  size_t pack_out(In_Seq const &In, Bool_Seq const &Fl, Out_Seq Out,
                flags fl = no_flag) {
    return pbbs::pack_out(In, Fl, Out, fl);
  }

  // Pack the output to the output range.
  template <SEQ Bool_Seq, RANGE Out_Seq>
  size_t pack_index_out(Bool_Seq const &Fl, Out_Seq Out,
                flags fl = no_flag) {
    using Idx_Type = typename Out_Seq::value_type;
    auto identity = [] (size_t i) {return (Idx_Type) i;};
    return pbbs::pack_out(pbbs::delayed_seq<Idx_Type>(Fl.size(),identity), Fl, Out, fl);
  }

  template <SEQ In_Seq, class F>
  auto filter(In_Seq const &In, F f, flags fl = no_flag)
    -> pbbs::sequence<typename In_Seq::value_type> {
      return pbbs::filter(In, f, fl);
  }

  template <SEQ In_Seq, RANGE Out_Seq, class F>
  size_t filter_out(In_Seq const &In, Out_Seq Out, F f, flags fl = no_flag) {
    return pbbs::filter_out(In, Out, f, fl);
  }

  template <class Idx_Type, SEQ Bool_Seq>
  pbbs::sequence<Idx_Type> pack_index(Bool_Seq const &Fl, flags fl = no_flag) {
    return pbbs::pack_index<Idx_Type>(Fl, fl);
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
    -> pbbs::sequence<typename Seq::value_type> {
      return pbbs::sample_sort(A, f, stable);
  }

  template<class Iter, typename BinPred>
  void sample_sort_inplace (pbbs::range<Iter> A, const BinPred& f, bool stable = false) {
    return pbbs::sample_sort_inplace(A, f, stable);
  }

  template<typename E, typename BinPred, typename s_size_t>
  void sample_sort (E* A, s_size_t n, const BinPred& f, bool stable) {
    return pbbs::sample_sort(A, n, f, stable);
  }


  // ====================== integer sort =======================
  template <typename T, typename Get_Key>
  void integer_sort_inplace(pbbs::range<T*> In,
          Get_Key const &g,
          size_t key_bits=0) {
    return pbbs::integer_sort_inplace(In, g, key_bits);
  }

  template <typename Seq, typename Get_Key>
  pbbs::sequence<typename Seq::value_type> integer_sort(Seq const &In, Get_Key
      const &g, size_t key_bits=0) {
    return pbbs::integer_sort(In, g, key_bits);
  }

  // Parallel version
  template <typename InS, typename KeyS>
  sequence<size_t> count_sort(InS const &In,
			      range<typename InS::value_type*> Out,
			      KeyS const &Keys,
			      size_t num_buckets,
			      float parallelism = 1.0,
			      bool skip_if_in_one=false) {
    return pbbs::count_sort(In, Out, Keys, num_buckets, parallelism, skip_if_in_one);
  }

  // ====================== random shuffle =======================
  using random = pbbs::random;

  template <class intT>
  pbbs::sequence<intT> random_permutation(size_t n, pbbs::random r = pbbs::random()) {
    return pbbs::random_permutation<intT>(n, r);
  }

  template <typename Seq>
  sequence<typename Seq::value_type>
  random_shuffle(Seq const &In, random r = random()) {
    return pbbs::random_shuffle<Seq>(In, r);
  }

}



// Other extensions to pbbs used by the graph benchmarks.
namespace pbbslib {

  constexpr size_t _F_BSIZE = 2000;

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

  inline int t_to_stringlen(long a) { return 21; }
  inline void type_to_string(char* s, long a) { sprintf(s, "%ld", a); }

  inline int t_to_stringlen(unsigned long a) { return 21; }
  inline void type_to_string(char* s, unsigned long a) { sprintf(s, "%lu", a); }

  inline uint t_to_stringlen(uint a) { return 12; }
  inline void type_to_string(char* s, uint a) { sprintf(s, "%u", a); }

  inline int t_to_stringlen(int a) { return 12; }
  inline void type_to_string(char* s, int a) { sprintf(s, "%d", a); }

  inline int t_to_stringlen(double a) { return 18; }

  inline int t_to_stringlen(char* a) { return strlen(a) + 1; }
  inline void type_to_string(char* s, char* a) { sprintf(s, "%s", a); }

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

  inline void type_to_string(char* s, double a) { sprintf(s, "%.11le", a); }

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
