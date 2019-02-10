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

  constexpr size_t radix = 8;
  constexpr size_t max_buckets = 1 << radix;

  template <class T, class F>
  void radix_step(T* A, T* B, size_t* counts, size_t n, size_t m, F extract) {
    for (size_t i = 0; i < m; i++)  counts[i] = 0;
    for (size_t j = 0; j < n; j++) {
      size_t k = extract(A[j]);
      counts[k]++;
    }

    size_t s = 0;
    for (size_t i = 0; i < m; i++) {
      s += counts[i];
      counts[i] = s;
    }
    
    for (long j = n-1; j >= 0; j--) {
      long x = --counts[extract(A[j])];
      B[x] = A[j];
    }
  }

  template <class T, class GetKey>
  void seq_radix_sort(sequence<T> const &In, sequence<T> &Out,
		      GetKey const &g,
		      size_t bits, bool inplace=true) {
    size_t n = In.size();
    if (n == 0) return;
    size_t counts[max_buckets+1];
    T* InA = In.begin();
    T* OutA = Out.begin();
    bool swapped = false;
    int bit_offset = 0;
    while (bits > 0) {
      size_t round_bits = std::min(radix, bits);
      size_t num_buckets = (1 << round_bits);
      size_t mask = num_buckets-1;
      auto get_key = [&] (T k) -> size_t {
	return (g(k) >> bit_offset) & mask;};
      radix_step(InA, OutA, counts, n, num_buckets, get_key);
      std::swap(InA,OutA);
      bits = bits - round_bits;
      bit_offset = bit_offset + round_bits;
      swapped = !swapped;
    }
    if ((inplace && swapped) || (!inplace && !swapped)) {
      for (size_t i=0; i < n; i++) 
	move_uninitialized(OutA[i], InA[i]);
    }
  }
  
  // a top down recursive radix sort
  // g extracts the integer keys from In
  // key_bits specifies how many bits there are left
  // if inplace is true, then result will be in In, otherwise in Out
  // In and Out cannot be the same
  template <typename T, typename Get_Key>
  void integer_sort_r(sequence<T> &In,  sequence<T> &Out,
		      Get_Key const &g, 
		      size_t key_bits, bool inplace) {
    size_t n = In.size();
    timer t;

    if (key_bits == 0) {
      if (!inplace)
	parallel_for(0, In.size(), [&] (size_t i) {Out[i] = In[i];});
      
    // for small inputs use sequential radix sort
	  } else if (n < (1 << 15)) {
      seq_radix_sort<T>(In, Out, g, key_bits, inplace);
    
    // few bits, just do a single parallel count sort
    } else if (key_bits <= radix) {
      t.start();
      size_t num_buckets = (1 << key_bits);
      size_t mask = num_buckets - 1;
      auto f = [&] (size_t i) {return g(In[i]) & mask;};
      auto get_bits = make_sequence<size_t>(n, f);
      if (inplace) count_sort(In, In, get_bits, num_buckets);
      else count_sort(In, Out, get_bits, num_buckets);
      
    // recursive case  
    } else {
      size_t bits = 8;
      size_t shift_bits = key_bits - bits;
      size_t buckets = (1 << bits);
      size_t mask = buckets - 1;
      auto f = [&] (size_t i) {return (g(In[i]) >> shift_bits) & mask;};
      auto get_bits = make_sequence<size_t>(n, f);

      // divide into buckets
      sequence<size_t> offsets = count_sort(In, Out, get_bits, buckets);

      // recursively sort each bucket
      parallel_for(0, buckets, [&] (size_t i) {
	size_t start = offsets[i];
	size_t end = offsets[i+1];
	auto a = Out.slice(start, end);
	auto b = In.slice(start, end);
	integer_sort_r(a, b, g, shift_bits, !inplace);
	}, 1);
    }
  }

  // a top down recursive radix sort
  // g extracts the integer keys from In
  // result will be placed in out, 
  //    but if inplace is true, then result will be put back into In
  // val_bits specifies how many bits there are in the key
  //    if set to 0, then a max is taken over the keys to determine
  template <typename T, typename Get_Key>
  void integer_sort(sequence<T> &In, sequence<T> &Out,
		    Get_Key const &g, 
		    size_t key_bits=0, bool inplace=false) {
    if (In.begin() == Out.begin()) {
      cout << "in integer_sort : input and output must be different locations" << endl;
      abort();}
    if (key_bits == 0) {
      using P = std::pair<size_t,size_t>;
      auto get_key = [&] (size_t i) -> P {
	size_t k =g(In[i]);
	return P(k,k);};
      auto keys = make_sequence<P>(In.size(), get_key);
      size_t min_val, max_val;
      std::tie(min_val,max_val) = reduce(keys, minmaxm<size_t>());
      key_bits = log2_up(max_val - min_val + 1);
      cout << key_bits << endl;
      if (min_val > max_val / 4) {
	auto h = [&] (T a) {return g(a) - min_val;};
	integer_sort_r(In, Out, h, key_bits, inplace);
	return;
      }
    }
    integer_sort_r(In, Out, g, key_bits, inplace);
  }

  template <typename T, typename Get_Key>
  void integer_sort(sequence<T> &In,
		    Get_Key const &g,
		    size_t key_bits=0) {
    sequence<T> Tmp = sequence<T>::alloc_no_init(In.size());
    integer_sort(In, Tmp, g, key_bits, true);
  }
}
