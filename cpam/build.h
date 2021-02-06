#pragma once
#include "parlay/primitives.h"
#include "get_time.h"

namespace cpam {

template <class Entry>
struct build {
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;
  using ET = typename Entry::entry_t;
  constexpr static auto less = [] (ET a, ET b) {
    return Entry::comp(Entry::get_key(a), Entry::get_key(b));};

  // sorts a sequence, then removes all but first element with equal keys
  // the sort is not necessarily stable, so any element could be kept
  template <class Seq>
  static parlay::sequence<ET> sort_remove_duplicates(Seq const &A) { // ?? const
    if (A.size() == 0) return parlay::sequence<ET>(0);
    auto B = parlay::internal::sample_sort(parlay::make_slice(A.begin(),A.end()), less);

    auto Fl = parlay::delayed_seq<bool>(B.size(), [&] (size_t i) {
	return (i==0) || less(B[i-1], B[i]); });

    return parlay::pack(B, Fl);
  }

  // // this version does the sort inplace if passed as an rvalue reference
  // template <class Seq>
  // static parlay::sequence<ET>
  // sort_remove_duplicates(Seq &&A, bool sequential=false) {
  //   if (A.size() == 0) return parlay::sequence<ET>(0);
  //   if (sequential) {
  //     parlay::internal::quicksort(A.begin(), A.size(), less);

  //     // remove duplicates
  //     size_t  j = 1;
  //     for (size_t i=1; i<A.size(); i++)
  // 	if (less(A[j-1], A[i])) A[j++] = A[i];
  //     return parlay::to_sequence(A.cut(0,j));
  //   } else {
  //     parlay::internal::sample_sort_inplace(parlay::make_slice(A), less);
  //     auto g = [&] (size_t i) {return (i==0) || less(A[i-1], A[i]); };
  //     auto C = parlay::pack(A, parlay::delayed_seq<bool>(A.size(),g));
  //     return C;
  //   }
  // }

  // // not sure what where this one is used??? remove?
  // parlay::sequence<ET>
  // static sort_remove_duplicates(ET* A, size_t n) {
  //   auto lessE = [&] (ET a, ET b) {
  //     return Entry::comp(a.first, b.first);};

  //   parlay::internal::sample_sort_inplace(parlay::make_slice(A,A+n), lessE);

  //   auto f = [&] (size_t i) -> ET {return A[i];};
  //   auto g = [&] (size_t i) {return (i==0) || lessE(A[i-1], A[i]); };
  //   auto B = parlay::pack(parlay::delayed_seq<ET>(n,f),
  // 			  parlay::delayed_seq<bool>(n,g));

  //   return B;
  // }

  template<class Seq, class Reduce>
  static parlay::sequence<ET>
  sort_reduce_duplicates(Seq const &A, const Reduce& reduce) {
    using E = typename Seq::value_type;
    using Vi = typename E::second_type;

    timer t("sort_reduce_duplicates", false);
    size_t n = A.size();
    if (n == 0) return parlay::sequence<ET>(0);
    auto lessE = [] (E const &a, E const &b) {
      return Entry::comp(a.first, b.first);};

    auto B = parlay::internal::sample_sort(parlay::make_slice(A), lessE);
    t.next("sort");

    // determines the index of start of each block of equal keys
    // and copies values into vals
    parlay::sequence<bool> Fl(n);
    auto Vals = parlay::tabulate(n, [&] (size_t i) -> Vi {
	Fl[i] = (i==0) || lessE(B[i-1], B[i]);
	return B[i].second;
      });
    t.next("copy, set flags");

    auto I = parlay::pack_index<node_size_t>(Fl);
    t.next("pack index");

    // combines over each block of equal keys using function reduce
    auto a = parlay::tabulate(I.size(), [&] (size_t i) -> ET {
	size_t start = I[i];
	size_t end = (i==I.size()-1) ? n : I[i+1];
	return ET(B[start].first, reduce(Vals.cut(start,end)));
      });
    t.next("reductions");
    // tabulate set over all entries of i
    return a;
  }

  template<class Seq, class Reduce>
  static parlay::sequence<ET>
  sort_reduce_duplicates_a(Seq const &A, const Reduce& reduce) {
    using E = typename Seq::value_type;
    using Vi = typename E::second_type;

    timer t("sort_reduce_duplicates", false);
    size_t n = A.size();
    if (n == 0) return parlay::sequence<ET>(0);
    auto lessE = [] (E const &a, E const &b) {
      return Entry::comp(a.first, b.first);};

    auto B = parlay::internal::sample_sort(parlay::make_slice(A), lessE);
    t.next("sort");

    // determines the index of start of each block of equal keys
    // and copies values into vals
    parlay::sequence<bool> Fl(n);  // ?? should it be unitialized
    auto Vals = parlay::tabulate(n, [&] (size_t i) -> Vi {
	Fl[i] = (i==0) || lessE(B[i-1], B[i]);
	return B[i].second;
      });
    t.next("copy, set flags");

    auto I = parlay::pack_index<node_size_t>(Fl);
    t.next("pack index");

    // combines over each block of equal keys using function reduce
    auto a = parlay::tabulate(I.size(), [&] (size_t i) -> ET {
	size_t start = I[i];
	size_t end = (i==I.size()-1) ? n : I[i+1];
	return ET(B[start].first, reduce(Vals.cut(start,end)));
      });
    t.next("reductions");
    // tabulate set over all entries of i
    return a;
  }

  template<class Seq, class Bin_Op>
  static parlay::sequence<ET>
  sort_combine_duplicates(Seq const &A, Bin_Op& f) {
    auto mon = parlay::make_monoid(f,V());
    auto reduce_op = [&] (parlay::slice<V*,V*> S) { // ?? should it be &
      return parlay::reduce(S, mon);};
    return sort_reduce_duplicates_a(A, reduce_op);
  }

  template<class Seq, class Bin_Op>
  static parlay::slice<ET*,ET*>
  sort_combine_duplicates_inplace(Seq const &A,  Bin_Op& f) {
    auto less = [&] (ET a, ET b) {return Entry::comp(a.first, b.first);};
    parlay::internal::quicksort(A.begin(), A.size(), less);
    size_t j = 0;
    for (size_t i=1; i < A.size(); i++) {
      if (less(A[j], A[i])) A[++j] = A[i];
      else A[j].second = f(A[j].second,A[i].second);
    }
    return A.cut(0,j+1);
  }

};

}  // namespace cpam
