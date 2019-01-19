// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011-2016 Guy Blelloch and the PBBS team
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
#include "utilities.h"
#include "get_time.h"

namespace pbbs {

  constexpr const size_t TRANS_THRESHHOLD = PAR_GRANULARITY/4;

  size_t split(size_t n) {
    return n/2;
    //return ((((size_t) 1) << log2_up(n) != n) ? n/2 : (7*(n+1))/16);
  }

  template <class E>
  struct transpose {
    E *A, *B;
    transpose(E *AA, E *BB) : A(AA), B(BB) {}

      void transR(size_t rStart, size_t rCount, size_t rLength,
		  size_t cStart, size_t cCount, size_t cLength) {
	if (cCount*rCount < TRANS_THRESHHOLD) {
	  for (size_t i=rStart; i < rStart+ rCount; i++)
	    for (size_t j=cStart; j < cStart + cCount; j++) 
	      B[j*cLength + i] = A[i*rLength + j];
	} else if (cCount > rCount) {
	  //size_t l1 = cCount/2;
	  size_t l1 = split(cCount);
	  size_t l2 = cCount - l1;
	  auto left = [&] () {
	    transR(rStart,rCount,rLength,cStart,l1,cLength);};
	  auto right = [&] () {
	    transR(rStart,rCount,rLength,cStart + l1,l2,cLength);};
	  par_do(left, right);
	} else {
	  //size_t l1 = rCount/2;
	  size_t l1 = split(cCount);
	  size_t l2 = rCount - l1;
	  auto left = [&] () {
	    transR(rStart,l1,rLength,cStart,cCount,cLength);};
	  auto right = [&] () {
	    transR(rStart + l1,l2,rLength,cStart,cCount,cLength);};
	  par_do(left, right);
	}	
      }

      void trans(size_t rCount, size_t cCount) {
#if defined(OPENMP)
#pragma omp parallel
#pragma omp single
#endif
	transR(0,rCount,cCount,0,cCount,rCount);
      }
    };

  template <class E, class int_t>
    struct blockTrans {
      E *A, *B;
      int_t *OA, *OB, *L;

    blockTrans(E *AA, E *BB, int_t *OOA, int_t *OOB, int_t *LL) 
    : A(AA), B(BB), OA(OOA), OB(OOB), L(LL) {}

      void transR(size_t rStart, size_t rCount, size_t rLength,
		  size_t cStart, size_t cCount, size_t cLength) {
	if (cCount*rCount < TRANS_THRESHHOLD/2) {
	  par_for(rStart, rStart+rCount, [&] (size_t i) {
	    for (size_t j=cStart; j < cStart + cCount; j++) {
	      E* pa = A+OA[i*rLength + j];
	      E* pb = B+OB[j*cLength + i];
	      size_t l = L[i*rLength + j];
	      for (size_t k =0; k < l; k++)
		move_uninitialized(pb[k], pa[k]);
	      //const size_t bytes = l * sizeof(E);
	      //for (size_t k=0; k < bytes; k++)
	      //((char*) pb)[k] = ((char *) pa)[k];
	    }
          });
	} else if (cCount > rCount) {
	  //size_t l1 = cCount/2;
	  size_t l1 = split(cCount);
	  size_t l2 = cCount - l1;
	  auto left = [&] () {
	    transR(rStart,rCount,rLength,cStart,l1,cLength);};
	  auto right = [&] () {
	    transR(rStart,rCount,rLength,cStart + l1,l2,cLength);};
	  par_do(left, right);
	} else {
	  //size_t l1 = rCount/2;
	  size_t l1 = split(cCount);
	  size_t l2 = rCount - l1;
	  auto left = [&] () {
	    transR(rStart,l1,rLength,cStart,cCount,cLength);};
	  auto right = [&] () {
	    transR(rStart + l1,l2,rLength,cStart,cCount,cLength);};
	  par_do(left, right);
	}	
      }
 
      void trans(size_t rCount, size_t cCount) {
#if defined(OPENMP)
#pragma omp parallel
#pragma omp single
#endif
	transR(0,rCount,cCount,0,cCount,rCount);
      }
  
    } ;

  // Moves values from blocks to buckets
  // From is sorted by key within each block, in block major
  // counts is the # of keys in each bucket for each block, in block major
  // From and To are of lenght n
  // counts is of length num_blocks * num_buckets
  // Data is memcpy'd into To avoiding initializers and overloaded =
  template<typename E, typename s_size_t>
  size_t* transpose_buckets(E* From, E* To, s_size_t* counts, size_t n,
			    size_t block_size,
			    size_t num_blocks, size_t num_buckets) {
    timer t;
    t.start();
    size_t m = num_buckets * num_blocks;
    sequence<s_size_t> dest_offsets(m);
    //if (n > 100000) std::cout << "blocks: " << num_blocks << " buckets: " << num_buckets << std::endl;

    // for smaller input do non-cache oblivious version
    if (n < (1 << 22) || num_buckets <= 512 || num_blocks <= 512) {
      size_t block_bits = log2_up(num_blocks);      
      size_t block_mask = num_blocks-1;
      if ((size_t) 1 << block_bits != num_blocks) {
	cout << "in transpose_buckets: num_blocks must be a power or 2"
	     << endl;
	abort();
      }

      // determine the destination offsets
      auto get = [&] (size_t i) {
	return counts[(i>>block_bits) + num_buckets*(i&block_mask)];};
      s_size_t sum = scan_add(make_sequence<s_size_t>(m, get), dest_offsets);
      if (sum != n) abort();

      // send each key to correct location within its bucket
      //parallel_for (size_t i = 0; i < num_blocks; i++) {
      auto f = [&] (size_t i) {
	size_t s_offset = i * block_size;
	for (size_t j= 0; j < num_buckets; j++) {
	  size_t d_offset = dest_offsets[i+ num_blocks*j];
	  size_t len = counts[i*num_buckets+j];
	  for (size_t k =0; k < len; k++)
	    move_uninitialized(To[d_offset++], From[s_offset++]);
	  //memcpy((char*) (To+d_offset), (char*) (From+s_offset),len*sizeof(E));
	}
      };
      par_for(0, num_blocks, 1, f);
    } else { // for larger input do cache efficient transpose
      sequence<s_size_t> source_offsets(m);
      sequence<s_size_t> seq_counts(counts,m);
      //std::cout << "blocks: " << num_blocks << " buckets: " << num_buckets << std::endl;
      //t.next("trans head");
      
      scan_add(seq_counts, source_offsets);
      transpose<s_size_t>(counts, dest_offsets.as_array()).trans(num_blocks,
								 num_buckets);
      //t.next("trans");
      scan_add(dest_offsets, dest_offsets);
      blockTrans<E,s_size_t>(From, To, source_offsets.as_array(),
			     dest_offsets.as_array(), counts).trans(num_blocks,
								    num_buckets);
      //t.next("block trans");
    }

    size_t *bucket_offsets = new_array_no_init<size_t>(num_buckets+1);
    for (s_size_t i=0; i < num_buckets; i++)
      bucket_offsets[i] = dest_offsets[i*num_blocks];
    // last element is the total size n
    bucket_offsets[num_buckets] = n;
    return bucket_offsets;
  }
}
