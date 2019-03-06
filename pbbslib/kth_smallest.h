#pragma once
#include "sequence_ops.h"
#include "random.h"

namespace pbbs {

  template <class Seq, class Compare>
  typename Seq::value_type kth_smallest(Seq const &s, size_t k, Compare less,
					random r = random()) {
    using T = typename Seq::value_type;
    size_t n = s.size();
    T pivot = s[r[0]%n];
    sequence<T> smaller = filter(s, [&] (T a) {return less(a, pivot);});
    if (k < smaller.size())
      return kth_smallest(smaller, k, less, r.next());
    else {
      sequence<T> larger = filter(s, [&] (T a) {return less(pivot, a);});
      if (k >= n - larger.size())
	return kth_smallest(larger, k - n + larger.size(), less, r.next());
      else return pivot;
    }
  }

  template <class Seq, class Compare>
  typename Seq::value_type approximate_kth_smallest(Seq S, size_t k, Compare less) {
    // raise exception if empty sequence?
    using T = typename Seq::value_type;
    pbbs::random r;
    size_t n = S.size();
    size_t samples = n/sqrt(n);
    pbbs::sequence<T> sample(samples, [&] (size_t i) -> T {
	return S[i];});
    return kth_smallest(sample, k * samples / n, less);
  }  
}
