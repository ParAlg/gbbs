// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2016 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include <stddef.h>

namespace pbbs {
// the following parameter can be tuned
constexpr const size_t _binary_search_base = 16;

template <typename Sequence, typename F>
size_t linear_search(Sequence I, typename Sequence::T v, const F& less) {
  for (size_t i = 0; i < I.size(); i++)
    if (!less(I[i], v)) return i;
  return I.size();
}

// return index to first key greater or equal to v
template <typename Sequence, typename F>
size_t binary_search(Sequence I, typename Sequence::T v, const F& less) {
  size_t start = 0;
  size_t end = I.size();
  while (end - start > _binary_search_base) {
    size_t mid = (end + start) / 2;
    if (!less(I[mid], v))
      end = mid;
    else
      start = mid + 1;
  }
  return start + linear_search(I.slice(start, end), v, less);
}
}  // namespace pbbs
