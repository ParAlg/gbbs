#pragma once

#include <stddef.h>

namespace pbbs {
// the following parameter can be tuned
constexpr const size_t _binary_search_base = 16;

template <typename Seq, typename F>
size_t linear_search(Seq const& I, typename Seq::value_type v, const F& less) {
  for (size_t i = 0; i < I.size(); i++)
    if (!less(I[i], v)) return i;
  return I.size();
}

// return index to first key greater or equal to v
template <typename Seq, typename F>
size_t binary_search(Seq const& I, typename Seq::value_type v, const F& less) {
  size_t start = 0;
  size_t end = I.size();
  while (end - start > _binary_search_base) {
    size_t mid = start + (end - start) / 2;
    if (!less(I[mid], v))
      end = mid;
    else
      start = mid + 1;
  }
  return start + linear_search(I.slice(start, end), v, less);
}

template <typename Seq, typename F>
size_t linear_search(Seq const& I, const F& less) {
  for (size_t i = 0; i < I.size(); i++)
    if (!less(I[i])) return i;
  return I.size();
}

// return index to first key where less is false
template <typename Seq, typename F>
size_t binary_search(Seq const& I, const F& less) {
  size_t start = 0;
  size_t end = I.size();
  while (end - start > _binary_search_base) {
    size_t mid = start + (end - start) / 2;
    if (!less(I[mid]))
      end = mid;
    else
      start = mid + 1;
  }
  return start + linear_search(I.slice(start, end), less);
}

template <typename Seq, typename F, typename G>
size_t binary_search_eq(Seq const& I, const F& less, const G& eq) {
  size_t start = 0;
  size_t end = I.size();
  while (end - start > _binary_search_base) {
    size_t mid = start + (end - start) / 2;
    if (eq(I[mid]))
      return mid;
    else if (!less(I[mid]))
      end = mid;
    else
      start = mid + 1;
  }
  return start + linear_search(I.slice(start, end), less);
}

}  // namespace pbbs
