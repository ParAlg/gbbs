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

  // a bottom up radix sort
  template <class Slice, class GetKey>
  void seq_radix_sort_(Slice In, Slice Out, GetKey const &g,
		       size_t bits, bool inplace=true) {
    size_t n = In.size();
    if (n == 0) return;
    size_t counts[max_buckets+1];
    bool swapped = false;
    int bit_offset = 0;
    while (bits > 0) {
      size_t round_bits = std::min(radix, bits);
      size_t num_buckets = (1 << round_bits);
      size_t mask = num_buckets-1;
      auto get_key = [&] (size_t i) -> size_t {
	return (g(In[i]) >> bit_offset) & mask;};
      seq_count_sort_(In, Out, delayed_seq<size_t>(n, get_key),
		      counts, num_buckets);
      std::swap(In,Out);
      bits = bits - round_bits;
      bit_offset = bit_offset + round_bits;
      swapped = !swapped;
    }
    if ((inplace && swapped) || (!inplace && !swapped)) {
      for (size_t i=0; i < n; i++) 
	move_uninitialized(Out[i], In[i]);
    }
  }

  // wrapper to reduce copies and avoid modifying In when not inplace
  // In and Tmp can be the same, but Out must be different
  template <class SeqIn, class Slice, class GetKey>
  void seq_radix_sort(SeqIn const &In, Slice Out, Slice Tmp, GetKey const &g,
		      size_t key_bits, bool inplace=true) {
    bool odd = ((key_bits-1)/radix) & 1;
    size_t n = In.size();
    if (slice_eq(In.slice(), Tmp)) { // not inplace
      if (odd) {
	for (size_t i=0; i < n; i++) 
	  move_uninitialized(Tmp[i], In[i]);
	seq_radix_sort_(Tmp, Out, g, key_bits, false);
      } else {
	for (size_t i=0; i < n; i++) 
	  move_uninitialized(Out[i], In[i]);
	seq_radix_sort_(Out, Tmp, g, key_bits, true);
      }
    }	else
      seq_radix_sort_(In.slice(), Out, g, key_bits, inplace);
  }


  // a top down recursive radix sort
  // g extracts the integer keys from In
  // key_bits specifies how many bits there are left
  // if inplace is true, then result will be in Tmp, otherwise in Out
  // In and Out cannot be the same, but In and Tmp should be same if inplace
  template <typename SeqIn, typename Slice, typename Get_Key>
  void integer_sort_r(SeqIn const &In, Slice Out, Slice Tmp, Get_Key const &g, 
		      size_t key_bits, bool inplace, float parallelism=1.0) {
    size_t n = In.size();
    timer t("integer sort",false);

    if (key_bits == 0) {
      if (!inplace)
	parallel_for(0, In.size(), [&] (size_t i) {Out[i] = In[i];});

      // for small inputs use sequential radix sort
    } else if (n < (1 << 15)) {
      seq_radix_sort(In, Out, Tmp, g, key_bits, inplace);
      
      // few bits, just do a single parallel count sort
    } else if (key_bits <= radix) {
      size_t num_buckets = (1 << key_bits);
      size_t mask = num_buckets - 1;
      auto f = [&] (size_t i) {return g(In[i]) & mask;};
      auto get_bits = delayed_seq<size_t>(n, f);
      sequence<size_t> cnts = count_sort(In.slice(), Out, get_bits, num_buckets,
					 parallelism, true);
      if (inplace != (cnts.size() == 0))
	parallel_for(0, n, [&] (size_t i) {
	    move_uninitialized(In[i], Out[i]);});
      
      // recursive case  
    } else {
      size_t bits = 8;
      size_t shift_bits = key_bits - bits;
      size_t buckets = (1 << bits);
      size_t mask = buckets - 1;
      auto f = [&] (size_t i) {return (g(In[i]) >> shift_bits) & mask;};
      auto get_bits = delayed_seq<size_t>(n, f);

      // divide into buckets
      sequence<size_t> offsets = count_sort(In.slice(), Out, get_bits, buckets,
					    parallelism, true);

      // if all but one bucket are empty, try again on lower bits
      if (offsets.size() == 0) {
	integer_sort_r(In, Out, Tmp, g, shift_bits, inplace,
		       parallelism);
	if (n > 10000000) t.next("isort");
	return;
      }
      if (n > 10000000) t.next("first");
      // recursively sort each bucket
      parallel_for(0, buckets, [&] (size_t i) {
	  size_t start = offsets[i];
	  size_t end = offsets[i+1];
	  auto a = Out.slice(start, end);
	  auto b = Tmp.slice(start, end);
	  integer_sort_r(a, b, a, g, shift_bits, !inplace,
			 (parallelism * (end - start)) / n);
	}, 1);
      if (n > 10000000) t.next("second");
    }
  }

  // a top down recursive radix sort
  // g extracts the integer keys from In
  // if inplace is false then result will be placed in Out, 
  //    otherwise they are placed in Tmp
  //    Tmp and In can be the same (i.e. to do inplace set them equal)
  // In is not directly modified, but can be indirectly if equal to Tmp
  // val_bits specifies how many bits there are in the key
  //    if set to 0, then a max is taken over the keys to determine
  template <typename SeqIn, typename IterOut, typename Get_Key>
    size_t integer_sort_(SeqIn const &In,
		     range<IterOut> Out,
		     range<IterOut> Tmp,
		     Get_Key const &g, 
		     size_t key_bits=0,
		     bool inplace=false) {
    if (slice_eq(In.slice(), Out)) {
      cout << "in integer_sort : input and output must be different locations"
	   << endl;
      abort();}
    size_t num_buckets = (1 << key_bits);
    if (key_bits == 0) {
      auto get_key = [&] (size_t i) {return g(In[i]);};
      auto keys = delayed_seq<size_t>(In.size(), get_key);
      num_buckets = reduce(keys, maxm<size_t>()) + 1;
      key_bits = log2_up(num_buckets);
    }
    integer_sort_r(In, Out, Tmp, g, key_bits, inplace);
    return num_buckets;
  }

  template <typename T, typename Get_Key>
  void integer_sort_inplace(range<T*> In,
			    Get_Key const &g,
			    size_t key_bits=0) {
    sequence<T> Tmp = sequence<T>::no_init(In.size());
    integer_sort_(In, Tmp.slice(), In, g, key_bits, true);
  }

  template <typename Seq, typename Get_Key>
  sequence<typename Seq::value_type> integer_sort(Seq const &In, Get_Key const &g,
						  size_t key_bits=0) {
    using T = typename Seq::value_type;
    sequence<T> Out = sequence<T>::no_init(In.size());
    sequence<T> Tmp = sequence<T>::no_init(In.size());
    integer_sort_(In, Out.slice(), Tmp.slice(), g, key_bits, false);
    return Out;
  }

  // Given a sorted sequence of integers in the range [0,..,num_buckets)
  // returns a sequence of length num_buckets+1 with the offset for the
  // start of each integer.   If an integer does not appear, its offset
  // will be the same as the next (i.e. offset[i+1]-offset[i] specifies
  // how many i there are.
  // The last element contains the size of the input.
  template <typename Seq, typename Get_Key>
  sequence<size_t>
  get_counts(Seq const &In, Get_Key const &g, size_t num_buckets) {
    size_t n = In.size();
    sequence<size_t> starts(n, (size_t) 0);
    sequence<size_t> ends(n, (size_t) 0);
    parallel_for (0, n-1, [&] (size_t i) {
	if (g(In[i]) != g(In[i+1])) {
	  starts[g(In[i+1])] = i+1;
	  ends[g(In[i])] = i+1;
	};});
    ends[g(In[n-1])] = n;
    return sequence<size_t>(n, [&] (size_t i) {
	return ends[i] - starts[i];});
  }

  template <typename Seq, typename Get_Key>
  std::pair<sequence<typename Seq::value_type>,sequence<size_t>>
  integer_sort_with_counts(Seq const &In, Get_Key const &g,
			    size_t num_buckets=0) {
    size_t key_bits = log2_up(num_buckets);
    auto R = integer_sort(In, g, key_bits);
    return std::make_pair(std::move(R), get_counts(R, g, num_buckets));
  }

}
