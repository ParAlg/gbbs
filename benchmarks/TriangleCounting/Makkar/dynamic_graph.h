// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once

#include <algorithm>
#include <cmath>
#include "gbbs/gbbs.h"
#include "pbbslib/sample_sort.h"

namespace gbbs {

  struct DynamicGraph {

    static constexpr size_t initial_vertex_size = 16;

    size_t n; // num_vertices
    pbbs::sequence<pbbs::sequence<uintE>> A; // adjacency lists
    pbbs::sequence<uintE> initial_vertex_memory; // the memory used to store adjacency lists of small vertices
    pbbs::sequence<bool> allocated; // a flag per vertex indicating whether 

    pbbs::sequence<std::pair<uintE, uintE>> post_ins_offs;

    using edge = gbbs::gbbs_io::Edge<int>;

    // For now, num_vertices is a fixed upper-bound on the number of vertices.
    // This code should be able to be easily extended to support vertex
    // insertions.
    DynamicGraph(size_t num_vertices) : n(num_vertices) {
      A = pbbs::sequence<pbbs::sequence<uintE>>(num_vertices);
      initial_vertex_memory = pbbs::sequence<uintE>(num_vertices*initial_vertex_size);
      allocated = pbbs::sequence<bool>(n);
      parallel_for(0, n, [&] (size_t i) {
        A[i] = pbbs::sequence<uintE>();
        allocated[i] = false;
      });
      std::cout << "DONE!" << std::endl;
    }

    ~DynamicGraph() {
      parallel_for(0, n, [&] (size_t i) {
        if (allocated[i]) { // TODO
//          ~A[i];
        }
        A[i].to_array(); // clears ownership
      });
    }

    // B is some container wrapping gbbs::gbbs_io::Edge<int> (a set of mixed
    // updates)
    template <class B>
    void process_batch(B& batch) {
      size_t batch_size = batch.size();

      auto insertions = pbbs::filter(batch, [&] (const auto& e) {
        return e.weight == 1;
      });
      auto deletions = pbbs::filter(batch, [&] (const auto& e) {
        return e.weight == 0;
      });

      process_insertions(insertions);
      process_deletions(std::move(deletions));
    }

    // An unsorted batch of updates
    template <class B>
    void process_insertions(B& batch) {
      auto sort_f = [&] (const edge& l, const edge& r) {
        return (l.from < r.from) || ((l.from == r.from) && (l.to < r.to));
      };
      timer sort_t; sort_t.start();
      pbbs::sample_sort_inplace(batch.slice(), sort_f);
      sort_t.stop(); sort_t.reportTotal("insert: sort time");

      // (i) generate vertex offsets into the batch
      auto index_seq = pbbs::delayed_seq<std::pair<uintE, size_t>>(batch.size(), [&] (size_t i) {
        return std::make_pair(batch[i].from, i);
      });

      // (vertex id, indexof start in batch) pairs.
      auto starts = pbbs::filter(index_seq, [&] (const std::pair<uintE, size_t>& p) {
        size_t ind = p.second;
        return ind == 0 || (index_seq[ind].first != index_seq[ind-1].first);
      });

      std::cout << "batchsize = " << batch.size() << " starts.size = " << starts.size() << std::endl;

      // now merge all
    }

    // An unsorted batch of updates
    template <class B>
    void process_deletions(B&& batch) {
    }

  };


}  // namespace gbbs
