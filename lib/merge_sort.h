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

#include "quicksort.h"  // needed for insertion_sort
#include "utilities.h"

namespace pbbs {

// If swap=0, reads from In and writes to Out,
// If swap=1, reads from In and writes to In (Out is temp space)
template <class SeqA, class SeqB, class F>
inline void merge_sort(SeqA Out, SeqB In, const F& f, bool swap = 0) {
  size_t n = In.size();
  if (n < 24) {
    pbbs::insertion_sort(In.get_array(), n, f);
    if (!swap)
      for (size_t i = 0; i < n; i++) Out[i] = In[i];  //.update(i, In[i]);
    return;
  }
  size_t m = n / 2;
  par_do(true, [&]() { merge_sort(Out.slice(0, m), In.slice(0, m), f, !swap); },
         [&]() { merge_sort(Out.slice(m, n), In.slice(m, n), f, !swap); });
  if (swap)
    pbbs::merge(Out.slice(0, m), Out.slice(m, n), In.slice(0, n), f);
  else
    pbbs::merge(In.slice(0, m), In.slice(m, n), Out.slice(0, n), f);
}
}  // namespace pbbs
