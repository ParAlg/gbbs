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

    // the number of vertices
    size_t n;

    // adjacency lists
    pbbs::sequence<uintE*> A;
    // degrees
    pbbs::sequence<uintE> D;
    // memory used to store adjacency lists of small vertices
    pbbs::sequence<uintE> initial_vertex_memory;
    // a flag per vertex indicating whether its adjacency list in A is allocated,
    // or from initial_vertex_memory
    pbbs::sequence<bool> allocated;
    // a dense array of offsets used to index into the batch
    // post_ins_offs[v] (if v is in the batch) indicates the index into "starts"
    // corresponding to v (see process_insertions).
    pbbs::sequence<std::pair<uintE, uintE>> post_ins_offs;

    // The edge type used. A batch is a sequence of edges.
    using edge = gbbs::gbbs_io::Edge<int>;

    // For now, num_vertices is a fixed upper-bound on the number of vertices.
    // This code should be able to be easily extended to support vertex
    // insertions.
    DynamicGraph(size_t num_vertices) : n(num_vertices) {
      A = pbbs::sequence<uintE*>(num_vertices);
      initial_vertex_memory = pbbs::sequence<uintE>(num_vertices*initial_vertex_size);
      allocated = pbbs::sequence<bool>(n);
      parallel_for(0, n, [&] (size_t i) {
        A[i] = &(initial_vertex_memory[i*initial_vertex_size]);
        allocated[i] = false;
        post_ins_offs[i] = std::make_pair(UINT_E_MAX, UINT_E_MAX); // TODO
      });
      std::cout << "Initialized!" << std::endl;
    }

    ~DynamicGraph() {
      // TODO
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

    template <class Container>
    size_t do_merge_ins(uintE* A, uintE nA, Container& B, uintE* out) {
      size_t a = 0, b = 0, nB = B.size(), k = 0;
      while (a < nA && b < nB) {
        if (A[a] == B[b].second) {
          out[k] = A[a];
          k++; a++; b++;
        } else if (A[a] < B[b].second) {
          out[k] = A[a];
          k++; a++;
        } else {
          out[k] = B[b];
          k++; b++;
        }
      }
      while (a < nA) {
        out[k] = A[a];
        k++; a++;
      }
      while (b < nB) {
        out[k] = B[b];
        k++; b++;
      }
      return k;
    }

    template <class C>
    void merge_small_insertions(uintE v, C& updates, uintE current_degree) {
      uintE[initial_vertex_size] tmp;
      uintE* ngh_v = A[v];
      size_t new_degree = do_merge_ins(ngh_v, current_degree, updates, &tmp);
      for (size_t i=0; i<new_degree; i++) {
        ngh_v[i] = tmp[i]; // update the actual edges
      }
      D[v] = new_degree; // update the degree in the table
    }

    template <class C>
    void merge_insertions(uintE v, C& updates) {
      uintE current_degree = D[v];
      uintE max_new_degree = current_degree + updates.size();
      if (max_new_degree <= initial_vertex_size) {
        merge_small_insertions(v, updates, current_degree);
        return;
      }
      // otherwise allocate
      uintE* ngh_v = A[v];
      uintE* new_array = pbbs::new_array_no_init<uintE>(max_new_degree);
      size_t new_degree = do_merge_ins(ngh_v, current_degree, updates, new_array);
      D[v] = new_degree; // update the degree in the table
    }

    // An unsorted batch of updates
    template <class B>
    void process_insertions(B& unsorted_batch) {
      // (i) sort the batch
      auto sort_f = [&] (const edge& l, const edge& r) {
        return (l.from < r.from) || ((l.from == r.from) && (l.to < r.to));
      };
      timer sort_t; sort_t.start();
      pbbs::sample_sort_inplace(unsorted_batch.slice(), sort_f);

      // (ii) define the unweighted version (with ins/del) dropped and filter to
      // remove duplicates in the input batch
      auto unsorted_batch_unweighted = pbbs::delayed_seq<std::pair<uintE, uintE>>(unsorted_batch.size(), [&] (const size_t i) {
        return std::make_pair(unsorted_batch[i].from, unsorted_batch[i].to);
      });
      auto batch = pbbs::filter_index(unsorted_batch_unweighted, [&] (const pair<uintE, uintE>& p, size_t ind) {
        return (ind == 0) || (p != unsorted_batch_unweighted[ind-1]);
      });

      // (iii) generate vertex offsets into the de-duplicated batch
      auto index_seq = pbbs::delayed_seq<std::pair<uintE, size_t>>(batch.size(), [&] (size_t i) {
        return std::make_pair(batch[i].first, i);
      });
      // starts contains (vertex id, indexof start in batch) pairs.
      auto starts = pbbs::filter(index_seq, [&] (const std::pair<uintE, size_t>& p) {
        size_t ind = p.second;
        return ind == 0 || (index_seq[ind].first != index_seq[ind-1].first)
      });
      std::cout << "batchsize = " << batch.size() << " starts.size = " << starts.size() << std::endl;

      parallel_for(0, starts.size(), [&] (size_t i) {
        const auto[v, index] = starts[i];
        post_ins_offs[v] = i; // TODO
      });

      // (iv) insert the updates into the old graph
      parallel_for(0, starts.size(), [&] (size_t i) {
        const auto [v, index] = starts[i];
        size_t v_deg = ((i == starts.size()) ? batch.size() : starts[i+].second) - index;
        if (v_deg == 0) abort();

        auto v_inserts = batch.slice(index, index + v_deg);
        merge_insertions(v, v_inserts);
      });

    }

    // An unsorted batch of updates
    template <class B>
    void process_deletions(B&& batch) {
    }

  };


}  // namespace gbbs
