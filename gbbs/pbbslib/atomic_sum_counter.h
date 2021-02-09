#pragma once

#include "gbbs/macros.h"

namespace pbbslib {

/* atomic max object for numeric type T */
/* Note the temporary bad hack around usage---see note on reset() */
template <class T>
struct atomic_sum_counter {
  T* entries;
  size_t stride;
  size_t num_elms;
  size_t num_workers_;
  atomic_sum_counter() {
    initialize();
//    num_workers_ = 0;
  }

  void initialize() {
    stride = 128/sizeof(T);
    stride = parlay::log2_up(stride);
    num_workers_ = num_workers();
    num_elms = num_workers_ << stride;
    entries = gbbs::new_array_no_init<T>(num_elms);
    for (size_t i=0; i<num_workers_; i++) {
      entries[i << stride] = (T)0;
    }
  }

  ~atomic_sum_counter() {
    if (entries != nullptr) {
      gbbs::free_array(entries, num_elms);
      entries = nullptr;
    }
  }

  void reset() {
    // bad hack; must call reset() before using.
    if (num_workers_ == 0) {
      initialize();
    } // else already initialized
    for (size_t i=0; i<num_workers_; i++) {
      entries[i << stride] = (T)0;
    }
  }

  T get_value() {
    T sum = 0;
    for (size_t i=0; i<num_workers_; i++) {
      sum += entries[i << stride];
    }
    return sum;
  }

  void update_value(T new_val) {
    size_t id = worker_id();
    entries[id << stride] += new_val;
  }
};

}  // namespace pbbslib
