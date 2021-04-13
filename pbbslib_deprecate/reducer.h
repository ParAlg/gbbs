#pragma once

#include <cilk/cilk.h>
#include <cilk/reducer.h>

namespace pbbs {

template <class T, int n>
struct histogram_view {
  using value_type = std::array<T, n>;
  value_type hist;

  histogram_view() {
    for (size_t i = 0; i < 6; i++) hist[i] = 0;
  }

  void reduce(histogram_view* right) {
    for (size_t i = 0; i < 6; i++) hist[i] += right->hist[i];
  }

  void add_value(size_t i) { hist[i] += 1; }
  value_type view_get_value() const { return hist; }
};

template <class T, int n>
using histogram_monoid = cilk::monoid_with_view<histogram_view<T, n>>;

template <class T, int n>
using histogram_reducer = cilk::reducer<histogram_monoid<T, n>>;

// use "histogram_reducer<int,n> r;" to define reducer with n int buckets
// use "r->add_value(i)" to increment bucket i
// use "r.get_value()[i]" to get bucket i

}  // namespace pbbs
