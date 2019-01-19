// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010-2017 Guy Blelloch and the PBBS team
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
#include <math.h>
#include <stdio.h>
#include <cstdint>
#include "utilities.h"
#include "counting_sort.h"
#include "quicksort.h"

namespace pbbs {
  
  // a top down recursive radix sort
  // g extracts the integer keys from In
  // val_bits specifies how many bits there are in the integer
  template <typename KT, typename InS, typename OutS, typename Get_Key>
  void integer_sort(InS In,  OutS Out, Get_Key& g, 
		    size_t val_bits, size_t depth=0) {
    using T = typename InS::T;
    size_t n = In.size();

    // for small inputs call a comparison sort
    if (n < (1 << 12) || depth > 2) {
      auto cmp = [&] (T a, T b) {return g(a) < g(b);};
      quicksort(In.start(),n,cmp);
      //std::sort(In.start(),In.end(),cmp);
      return;
    }
    size_t bits;  // 2^bits = number of buckets for first round
    if (depth == 0) bits = log2_up(n)/3; // cube root
    else bits = log2_up(n)/2;            // square root
    if (sizeof(T) <= 4) bits += 1;

    // for larger elements only one level of counting sort
    if (sizeof(T) > 8) {bits = bits+2; depth = 3;}

    // never do more bits that in the key
    bits = std::min(bits, val_bits);

    size_t num_buckets = (1<<bits);
    size_t mask = num_buckets - 1;
    size_t shift_bits = val_bits - bits;

    auto high_bits = [&] (size_t i) {
      return (g(In[i]) >> shift_bits) & mask;};

    auto get_bits = make_sequence<size_t>(n, high_bits);
    
    // sort into buckets using a counting sort
    sequence<size_t> bucket_offsets 
      = count_sort(In, Out, get_bits, num_buckets);

    // recurse on each bucket
    if (shift_bits > 0) {
      //parallel_for(size_t i = 0; i < num_buckets; i++) {
      auto f = [&] (size_t i) {
	auto out_slice = Out.slice(bucket_offsets[i],bucket_offsets[i+1]);
	integer_sort<KT>(out_slice, out_slice, g, shift_bits, depth+1);
      };
      par_for(0, num_buckets, 1, f);
    }
  }
}
