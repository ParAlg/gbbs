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
#include "counting_sort.h"
#include "quicksort.h"
#include "utilities.h"

namespace pbbs {

constexpr size_t radix = 8;
constexpr size_t max_buckets = 1 << radix;

// a bottom up radix sort
template <class Slice, class GetKey>
void seq_radix_sort_(Slice In, Slice Out, GetKey const &g, size_t bits,
                     bool inplace = true) {
  size_t n = In.size();
  if (n == 0) return;
  size_t counts[max_buckets + 1];
  bool swapped = false;
  int bit_offset = 0;
  while (bits > 0) {
    size_t round_bits = std::min(radix, bits);
    size_t num_buckets = (1 << round_bits);
    size_t mask = num_buckets - 1;
    auto get_key = [&](size_t i) -> size_t {
      return (g(In[i]) >> bit_offset) & mask;
    };
    seq_count_sort_(In, Out, delayed_seq<size_t>(n, get_key), counts,
                    num_buckets);
    std::swap(In, Out);
    bits = bits - round_bits;
    bit_offset = bit_offset + round_bits;
    swapped = !swapped;
  }
  if ((inplace && swapped) || (!inplace && !swapped)) {
    for (size_t i = 0; i < n; i++) move_uninitialized(Out[i], In[i]);
  }
}

// wrapper to reduce copies and avoid modifying In when not inplace
// In and Tmp can be the same, but Out must be different
template <class SeqIn, class Slice, class GetKey>
void seq_radix_sort(SeqIn const &In, Slice Out, Slice Tmp, GetKey const &g,
                    size_t key_bits, bool inplace = true) {
  bool odd = ((key_bits - 1) / radix) & 1;
  size_t n = In.size();
  if (slice_eq(In.slice(), Tmp)) {  // inplace
    seq_radix_sort_(Tmp.slice(), Out, g, key_bits, inplace);
  } else {
    if (odd) {
      for (size_t i = 0; i < n; i++) move_uninitialized(Tmp[i], In[i]);
      seq_radix_sort_(Tmp, Out, g, key_bits, false);
    } else {
      for (size_t i = 0; i < n; i++) move_uninitialized(Out[i], In[i]);
      seq_radix_sort_(Out, Tmp, g, key_bits, true);
    }
  }
}

// a top down recursive radix sort
// g extracts the integer keys from In
// key_bits specifies how many bits there are left
// if inplace is true, then result will be in Tmp, otherwise in Out
// In and Out cannot be the same, but In and Tmp should be same if inplace
template <typename SeqIn, typename Slice, typename Get_Key>
sequence<size_t> integer_sort_r(SeqIn const &In, Slice Out, Slice Tmp,
                                Get_Key const &g, size_t key_bits,
                                size_t num_buckets, bool inplace,
                                float parallelism = 1.0) {
  using T = typename SeqIn::value_type;
  size_t n = In.size();
  timer t("integer sort", false);
  size_t cache_per_thread = 1000000;
  size_t base_bits = log2_up(2 * (size_t)sizeof(T) * n / cache_per_thread);
  // keep between 8 and 13
  base_bits = std::max<size_t>(8, std::min<size_t>(13, base_bits));
  sequence<size_t> offsets;
  bool one_bucket;
  bool return_offsets = (num_buckets > 0);

  if (key_bits == 0) {
    if (!inplace) parallel_for(0, In.size(), [&](size_t i) { Out[i] = In[i]; });
    return sequence<size_t>();

    // for small inputs or little parallelism use sequential radix sort
  } else if ((n < (1 << 17) || parallelism < .0001) && !return_offsets) {
    seq_radix_sort(In, Out, Tmp, g, key_bits, inplace);
    return sequence<size_t>();

    // few bits, just do a single parallel count sort
  } else if (key_bits <= base_bits) {
    size_t mask = (1 << key_bits) - 1;
    auto f = [&](size_t i) { return g(In[i]) & mask; };
    auto get_bits = delayed_seq<size_t>(n, f);
    size_t num_bkts = (num_buckets == 0) ? (1 << key_bits) : num_buckets;
    // only uses one bucket optimization (last argument) if inplace
    std::tie(offsets, one_bucket) =
        count_sort(In.slice(), Out, get_bits, num_bkts, parallelism, inplace);
    if (inplace && !one_bucket)
      parallel_for(0, n, [&](size_t i) { move_uninitialized(Tmp[i], Out[i]); });
    if (return_offsets)
      return offsets;
    else
      return sequence<size_t>();

    // recursive case
  } else {
    size_t bits = 8;
    size_t shift_bits = key_bits - bits;
    size_t num_outer_buckets = (1 << bits);
    size_t num_inner_buckets = return_offsets ? ((size_t)1 << shift_bits) : 0;
    size_t mask = num_outer_buckets - 1;
    auto f = [&](size_t i) { return (g(In[i]) >> shift_bits) & mask; };
    auto get_bits = delayed_seq<size_t>(n, f);

    // divide into buckets
    std::tie(offsets, one_bucket) =
        count_sort(In.slice(), Out, get_bits, num_outer_buckets, parallelism,
                   !return_offsets);

    // if all but one bucket are empty, try again on lower bits
    if (one_bucket) {
      return integer_sort_r(In, Out, Tmp, g, shift_bits, 0, inplace,
                            parallelism);
    }

    sequence<size_t> inner_offsets(return_offsets ? num_buckets + 1 : 0);
    if (return_offsets) inner_offsets[num_buckets] = n;

    // recursively sort each bucket
    parallel_for(
        0, num_outer_buckets,
        [&](size_t i) {
          size_t start = offsets[i];
          size_t end = offsets[i + 1];
          auto a = Out.slice(start, end);
          auto b = Tmp.slice(start, end);
          sequence<size_t> r =
              integer_sort_r(a, b, a, g, shift_bits, num_inner_buckets,
                             !inplace, (parallelism * (end - start)) / (n + 1));
          if (return_offsets) {
            size_t bstart = std::min(i * num_inner_buckets, num_buckets);
            size_t bend = std::min((i + 1) * num_inner_buckets, num_buckets);
            size_t m = (bend - bstart);
            for (size_t j = 0; j < m; j++)
              inner_offsets[bstart + j] = offsets[i] + r[j];
          }
        },
        1);
    return inner_offsets;
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
// If num_buckets is non-zero then the output sequence will contain
// the offsets of each bucket (num_bucket of them)
// num_bucket must be less than or equal to 2^bits
template <typename SeqIn, typename IterOut, typename Get_Key>
sequence<size_t> integer_sort_(SeqIn const &In, range<IterOut> Out,
                               range<IterOut> Tmp, Get_Key const &g,
                               size_t bits, size_t num_buckets, bool inplace) {
  if (slice_eq(In.slice(), Out)) {
    std::cout << "in integer_sort : input and output must be different locations" << std::endl;
    exit(-1);
  }
  if (bits == 0) {
    auto get_key = [&](size_t i) { return g(In[i]); };
    auto keys = delayed_seq<size_t>(In.size(), get_key);
    num_buckets = reduce(keys, maxm<size_t>()) + 1;
    bits = log2_up(num_buckets);
  }
  return integer_sort_r(In, Out, Tmp, g, bits, num_buckets, inplace);
}

template <typename T, typename Get_Key>
void integer_sort_inplace(range<T *> In, Get_Key const &g,
                          size_t num_bits = 0) {
  sequence<T> Tmp = sequence<T>::no_init(In.size());
  integer_sort_(In, Tmp.slice(), In, g, num_bits, 0, true);
}

template <typename Seq, typename Get_Key>
sequence<typename Seq::value_type> integer_sort(Seq const &In, Get_Key const &g,
                                                size_t num_bits = 0) {
  using T = typename Seq::value_type;
  sequence<T> Out = sequence<T>::no_init(In.size());
  sequence<T> Tmp = sequence<T>::no_init(In.size());
  integer_sort_(In, Out.slice(), Tmp.slice(), g, num_bits, 0, false);
  return Out;
}

// Given a sorted input sequence with integer keys in the range
// `[0,..,num_buckets)`, returns a sequence of length `num_buckets` in which the
// i-th element is the number of times key i appears in the input.
template <typename Tint = size_t, typename Seq, typename Get_Key>
sequence<Tint> get_counts(Seq const &In, Get_Key const &g, size_t num_buckets) {
  size_t n = In.size();
  if (n == 0) {
    return sequence<Tint>(num_buckets, (Tint)0);
  }
  sequence<Tint> starts(num_buckets, (Tint)0);
  sequence<Tint> ends(num_buckets, (Tint)0);
  parallel_for(0, n - 1, [&](size_t i) {
    if (g(In[i]) != g(In[i + 1])) {
      starts[g(In[i + 1])] = i + 1;
      ends[g(In[i])] = i + 1;
    };
  });
  ends[g(In[n - 1])] = n;
  return sequence<Tint>(num_buckets,
                        [&](size_t i) { return ends[i] - starts[i]; });
}

template <typename Tint = size_t, typename Seq, typename Get_Key>
std::pair<sequence<typename Seq::value_type>, sequence<Tint>>
integer_sort_with_counts(Seq const &In, Get_Key const &g, size_t num_buckets) {
  size_t bits = log2_up(num_buckets);
  auto R = integer_sort(In, g, bits);
  return std::make_pair(std::move(R), get_counts<Tint>(R, g, num_buckets));
}

}  // namespace pbbs
