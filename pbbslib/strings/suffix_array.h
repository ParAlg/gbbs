// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2019 Guy Blelloch, Julian Shun and the PBBS team
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

// A modified version of the Apostolico, Iliopoulos, Landau, Schieber, and Vishkin
// Suffix Array algorithm (actually originally designed as a suffix tree algorithm)
// It does O(n log n) work in the worst case, but for most inputs
//    it does O(n) work beyond a sort on fixed length keys.
// The depth is O(log^2 n) assuming the sort is within that bound.
// Each round doubles the substring lengths that have been sorted but drops any strings
// that are already sorted (i.e. have no equal strings for the current length).
// For most inputs, most strings (suffixes) drop out early.

// Supports the following interface returning a suffix array for c
//   indexT is the type of integer for the suffix indices
//   It needs to have at least log_2 n bits.  It can be unsigned.

//  template <typename indexT>
//  pbbs::sequence<indexT> pbbs::suffixArray(pbbs::sequence<unsigned char> const &s);
#pragma once

#include "sequence.h"
#include <math.h>
#include "parallel.h"
#include "get_time.h"
#include "sample_sort.h"

namespace pbbs {

  constexpr bool verbose = false;

  using namespace std;

  using uchar = unsigned char;
  using uint128 = unsigned __int128;

  template <typename indexT>
  using ipair = std::pair<indexT,indexT>;

  timer sa_timer("Suffix Array", false);

  template <typename indexT>
  struct seg {
    indexT start;
    indexT length;
    seg<indexT>(indexT s, indexT l) : start(s), length(l) {}
    seg<indexT>() {}
  };

  template <typename indexT>
  void split_segment(range<seg<indexT>*> segOut,
		    indexT start,
		    sequence<indexT> &ranks,
		    range<ipair<indexT>*> Cs) {
    indexT l = segOut.size();
    if (l < 5000) { // sequential version

      indexT name = 0;
      ranks[Cs[0].second] = name + start + 1;
      for (indexT i=1; i < l; i++) {
	if (Cs[i-1].first != Cs[i].first) name = i;
	ranks[Cs[i].second] = name + start + 1;
      }

      name = 0;
      for (indexT i=1; i < l; i++) {
	if (Cs[i-1].first != Cs[i].first) {
	  segOut[i-1] = seg<indexT>(name+start,i-name);
	  name = i;
	} else segOut[i-1] = seg<indexT>(0,0);
      }
      segOut[l-1] = seg<indexT>(name+start,l-name);

    } else { // parallel version
      sequence<indexT> names(l);

      // mark start of each segment with equal keys
      parallel_for (1, l, [&] (size_t i) {
	  names[i] = (Cs[i].first != Cs[i-1].first) ? i : 0;});
      names[0] = 0;

      // scan start i across each segment
      scan_inplace(names.slice(), maxm<indexT>(), fl_scan_inclusive);

      // write new rank into original location
      parallel_for (0, l, [&] (size_t i) {
	  ranks[Cs[i].second] = names[i]+start+1;});

      // get starts and lengths of new segments
      parallel_for (1, l, [&] (size_t i) {
	  if (names[i] == i)
	    segOut[i-1] = seg<indexT>(start+names[i-1],i-names[i-1]);
	  else segOut[i-1] = seg<indexT>(0,0);
	});
      segOut[l-1] = seg<indexT>(start+names[l-1],l-names[l-1]);
    }
  }

  template <class indexT>
  sequence<ipair<indexT>>
  split_segment_top(sequence<seg<indexT>> &segOut,
		  sequence<indexT> &ranks,
		  sequence<uint128> const &Cs) {
    size_t n = segOut.size();
    sequence<indexT> names(n);
    size_t mask = ((((size_t) 1) << 32) - 1);

    // mark start of each segment with equal keys
    parallel_for (1, n, [&] (size_t i) {
	names[i] = ((Cs[i] >> 32) != (Cs[i-1] >> 32)) ? i : 0;});
    names[0] = 0;
    sa_timer.next("names");

    // scan start i across each segment
    scan_inplace(names.slice(), maxm<indexT>(), fl_scan_inclusive);

    sequence<ipair<indexT>> C(n);
    // write new rank into original location
    parallel_for (0, n, [&] (size_t i) {
	ranks[Cs[i] & mask] = names[i]+1;
	C[i].second = Cs[i] & mask;
      });
    sa_timer.next("write rank and copy");

    // get starts and lengths of new segments
    parallel_for (1, n, [&] (size_t i) {
	if (names[i] == i)
	  segOut[i-1] = seg<indexT>(names[i-1],i-names[i-1]);
	else segOut[i-1] = seg<indexT>(0,0);
      });
    segOut[n-1] = seg<indexT>(names[n-1],n-names[n-1]);

    sa_timer.next("segments");
    return C;
  }

  template <class indexT>
  sequence<indexT> suffix_array(sequence<uchar> const &ss) {
    if (verbose) sa_timer.start();
    size_t n = ss.size();

    // renumber characters densely
    // start numbering at 1 leaving 0 to indicate end-of-string
    size_t pad = 48;
    sequence<indexT> flags(256, (indexT) 0);
    parallel_for (0, n, [&] (size_t i) {
	if (!flags[ss[i]]) flags[ss[i]] = 1;}, 1000);
    auto add = [&] (indexT a, indexT b) {return a + b;};
    indexT m;
    std::tie(flags, m) = scan(flags, make_monoid(add,(indexT) 1));

    // pad the end of string with 0s
    sequence<uchar> s(n + pad, [&] (size_t i) {
	return (i < n) ? flags[ss[i]] : 0;});

    if (verbose) std::cout << "distinct characters = " << m-1 << std::endl;

    // pack characters into 128-bit word, along with the location i
    // 96 bits for characters, and 32 for location
    double logm = log2((double) m);
    indexT nchars = floor(96.0/logm);

    sequence<uint128> Cl(n, [&] (size_t i) {
	uint128 r = s[i];
	for (indexT j=1; j < nchars; j++) r = r*m + s[i+j];
	return (r << 32) + i;
      });
    sa_timer.next("copy into 128bit int");

    // sort based on packed words
    sample_sort_inplace(Cl.slice(), std::less<uint128>());
    sa_timer.next("sort");

    // identify segments of equal values
    sequence<indexT> ranks(n);
    sequence<seg<indexT>> seg_outs(n);
    sequence<ipair<indexT>> C = split_segment_top(seg_outs, ranks, Cl);
    Cl.clear();
    sa_timer.next("split");

    indexT offset = nchars;
    uint round =0;
    indexT nKeys = n;

    // offset is how many characters for each suffix have already been sorted
    // each round doubles offset so there should be at most log n rounds
    // The segments keep regions that have not yet been fully sorted
    while (1) {
      if (round++ > 40) {
	cout << "Suffix Array:  Too many rounds" << std::endl;
	abort();
      }

      auto is_seg = [&] (seg<indexT> s) {return s.length > 1;};
      // only keep segments that are longer than 1 (otherwise already sorted)
      sequence<seg<indexT>> Segs = filter(seg_outs.slice(0,nKeys), is_seg);
      indexT nSegs = Segs.size();
      if (nSegs == 0) break;
      sa_timer.next("filter and scan");

      sequence<indexT> offsets(nSegs);
      parallel_for (0, nSegs, [&] (size_t i) {
	  indexT start = Segs[i].start;
	  indexT l = Segs[i].length;
	  auto Ci = C.slice(start, start + l);
	  offsets[i] = l;

	  // grab rank from offset locations ahead
	  parallel_for (0, l, [&] (size_t j) {
	      indexT o = Ci[j].second + offset;
	      Ci[j].first = (o >= n) ? 0 : ranks[o];
	    }, 100);

	  // sort within each segment based on ranks
	  auto less = [&] (ipair<indexT> A, ipair<indexT> B) {
	    return A.first < B.first;};
	  if (l >= n/10) sample_sort_inplace(Ci, less);
	  else quicksort(Ci, less);
	});
      sa_timer.next("sort");

      // starting offset for each segment
      nKeys = scan_inplace(offsets.slice(), addm<indexT>());

      // Split each segment into subsegments if neighbors differ.
      parallel_for (0, nSegs, [&] (size_t i) {
	  indexT start = Segs[i].start;
	  indexT l = Segs[i].length;
	  indexT o = offsets[i];
	  split_segment(seg_outs.slice(o, o + l),
			start,
			ranks,
			C.slice(start, start+l));
	}, 100);
      sa_timer.next("split");

      if (verbose)
	cout << "length: " << offset << " keys remaining: " << nKeys << std::endl;

      offset = 2 * offset;
    }
    parallel_for (0, n, [&] (size_t i) {
	ranks[i] = C[i].second;});
    return ranks;
  }

}  // namespace pbbs
