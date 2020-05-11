#pragma once

#include "pbbslib/utilities.h"

/* atomic max object for numeric type T */
template <class T>
struct atomic_max_counter {
  T entry;
  atomic_max_counter() {
    entry = (T)0;
  }

  void reset() {
    entry = (T)0;
  }

  T get_value() {
    return entry;
  }

  void update_value(T new_val) {
    pbbs::write_min(&entry, new_val, std::greater<T>());
  }
};

