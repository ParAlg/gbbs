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

    // the triangle count
    size_t T;

    // adjacency lists
    pbbs::sequence<uintE*> A;
    // degrees
    pbbs::sequence<uintE> D;
    // memory used to store adjacency lists of small vertices
    pbbs::sequence<uintE> initial_vertex_memory;
    // a flag per vertex indicating whether its adjacency list in A is allocated,
    // or from initial_vertex_memory
    pbbs::sequence<bool> allocated;
    // a dense array of offsets used to index into the batch. UINT_E_MAX if the
    // vertex is not updated in this batch.
    pbbs::sequence<uintE> starts_offsets;

    // The edge type used. A batch is a sequence of edges.
    using edge = gbbs::gbbs_io::Edge<int>;

    // For now, num_vertices is a fixed upper-bound on the number of vertices.
    // This code should be able to be easily extended to support vertex
    // insertions.
    DynamicGraph(size_t num_vertices) : n(num_vertices) {
      A = pbbs::sequence<uintE*>(num_vertices);
      D = pbbs::sequence<uintE>(n);
      T = 0;
      initial_vertex_memory = pbbs::sequence<uintE>(num_vertices*initial_vertex_size);
      allocated = pbbs::sequence<bool>(n);
      starts_offsets = pbbs::sequence<uintE>(n);
      parallel_for(0, n, [&] (size_t i) {
        A[i] = &(initial_vertex_memory[i*initial_vertex_size]);
        D[i] = 0;
        allocated[i] = false;
        starts_offsets[i] = UINT_E_MAX;
        // post_ins_offs[i] = std::make_pair(UINT_E_MAX, UINT_E_MAX); // TODO
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
      auto insertions = pbbs::filter(batch, [&] (const auto& e) {
        return e.weight == 1;
      });
      auto deletions = pbbs::filter(batch, [&] (const auto& e) {
//        return e.weight == 0;
        return e.weight == 1;
      });

      process_insertions(insertions);

      process_deletions(deletions);

      process_insertions(insertions);
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
          out[k] = B[b].second;
          k++; b++;
        }
      }
      while (a < nA) {
        out[k] = A[a];
        k++; a++;
      }
      while (b < nB) {
        out[k] = B[b].second;
        k++; b++;
      }
      return k;
    }

    template <class C>
    void merge_small_insertions(uintE v, C& updates, uintE current_degree) {
      uintE tmp[initial_vertex_size];
      uintE* ngh_v = A[v];
      size_t new_degree = do_merge_ins(ngh_v, current_degree, updates, tmp);
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
      if (allocated[v]) { pbbs::free_array(ngh_v); } // free old neighbors if nec.
      else { allocated[v] = true; } // update allocated[v] if nec.
      A[v] = new_array; // update to the new neighbors
    }

    size_t intersect_A(uintE a_id, uintE b_id) {
      size_t a = 0, b = 0, nA = D[a_id], nB = D[b_id], k = 0;
      uintE* A_arr = A[a_id]; uintE* B_arr = A[b_id];
      while (a < nA && b < nB) {
        if (A_arr[a] == B_arr[b]) {
          k++; a++; b++;
        } else if (A_arr[a] < B_arr[b]) {
          a++;
        } else {
          b++;
        }
      }
      return k;
    }

    template <class Container>
    size_t truncated_intersect(uintE a_id, Container& B) {
      size_t a = 0, b = 0, nA = D[a_id], nB = B.size(), k = 0;
      uintE* A_arr = A[a_id];
      while (a < nA && b < nB) {
        if (A_arr[a] <= a_id) { // filter edges to vertices < a_id
          a++;
          continue;
        }
        if (A_arr[a] == B[b].second) {
          k++; a++; b++;
        } else if (A_arr[a] < B[b].second) {
          a++;
        } else {
          b++;
        }
      }
      return k;
    }

    template <class ContainerA, class ContainerB>
    size_t intersect_new(ContainerA& A, ContainerB& B) {
      size_t a = 0, b = 0, nA = A.size(), nB = B.size(), k = 0;
      while (a < nA && b < nB) {
        if (A[a].second == B[b].second) {
          k++; a++; b++;
        } else if (A[a].second < B[b].second) {
          a++;
        } else {
          b++;
        }
      }
      return k;
    }

    // An unsorted batch of updates
    template <class B>
    void process_insertions(B& unsorted_batch) {
      auto duplicated_batch = pbbs::sequence<edge>(2*unsorted_batch.size());
      parallel_for(0, unsorted_batch.size(), [&] (size_t i) {
        duplicated_batch[2*i] = unsorted_batch[i];
        duplicated_batch[2*i+1].from = unsorted_batch[i].to;
        duplicated_batch[2*i+1].to = unsorted_batch[i].from;
        duplicated_batch[2*i+1].weight = 1;
      });

      // (i) sort the batch
      auto sort_f = [&] (const edge& l, const edge& r) {
        return (l.from < r.from) || ((l.from == r.from) && (l.to < r.to));
      };
      timer sort_t; sort_t.start();
      pbbs::sample_sort_inplace(duplicated_batch.slice(), sort_f);

      // (ii) define the unweighted version (with ins/del) dropped and filter to
      // remove duplicates in the batch
      auto duplicated_batch_unweighted = pbbs::delayed_seq<std::pair<uintE, uintE>>(duplicated_batch.size(), [&] (const size_t i) {
        return std::make_pair(duplicated_batch[i].from, duplicated_batch[i].to);
      });
      auto batch = pbbs::filter_index(duplicated_batch_unweighted, [&] (const pair<uintE, uintE>& p, size_t ind) {
        return (ind == 0) || (p != duplicated_batch_unweighted[ind-1]);
      });

      // (iii) generate vertex offsets into the de-duplicated batch
      auto index_seq = pbbs::delayed_seq<std::pair<uintE, size_t>>(batch.size(), [&] (size_t i) {
        return std::make_pair(batch[i].first, i);
      });
      // starts contains (vertex id, indexof start in batch) pairs.
      auto starts = pbbs::filter(index_seq, [&] (const std::pair<uintE, size_t>& p) {
        size_t ind = p.second;
        return ind == 0 || (index_seq[ind].first != index_seq[ind-1].first);
      });
      std::cout << "batchsize = " << batch.size() << " starts.size = " << starts.size() << std::endl;

      // (iv) update starts_offsets mapping for batch vertices (reset at the end of this batch)
      parallel_for(0, starts.size(), [&] (size_t i) {
        const auto[v, index] = starts[i];
        starts_offsets[v] = i;
      });

      // (v) insert the updates into the old graph
      parallel_for(0, starts.size(), [&] (size_t i) {
        const auto [v, index] = starts[i];
        size_t v_deg = ((i == starts.size()-1) ? batch.size() : starts[i+1].second) - index;
        if (v_deg == 0) abort();

        auto v_inserts = batch.slice(index, index + v_deg);
        merge_insertions(v, v_inserts);
      });

      /////////////////////////////////////////////////////////////////////////
      // the next phase of the algorithm intersects the batch + the updated
      // graph, and updates the triangle counts

      // G(u) intersect G(v)
      pbbs::sequence<size_t> counts_one(batch.size());
      parallel_for(0, batch.size(), [&] (size_t b) {
        auto [u, v] = batch[b];
        if (u < v) {
          counts_one[b] = intersect_A(u, v);
        } else {
          counts_one[b] = 0;
        }
      });

      auto get_new_edgelist = [&] (const uintE& u, size_t u_starts_offset) {
        auto [up, index] = starts[u_starts_offset];
        assert(up == u);
        size_t u_deg = ((u_starts_offset == starts.size()-1) ? batch.size() : starts[u_starts_offset+1].second) - index;
        if (u_deg == 0) abort();
        auto u_inserts = batch.slice(index, index + u_deg);
        return u_inserts;
      };

      // truncated G(u) intersect G'(v) and truncated G(v) intersect G'(u)
      pbbs::sequence<size_t> counts_two(batch.size());
      parallel_for(0, batch.size(), [&] (size_t b) {
        auto [u, v] = batch[b];
        size_t count = 0;

        if (u < v) {
          size_t v_starts_offset = starts_offsets[v];
          if (v_starts_offset != UINT_E_MAX) { // non-zero number of updates for v in this batch
            auto v_inserts = get_new_edgelist(v, v_starts_offset);
            count += truncated_intersect(u, v_inserts);
          }

          size_t u_starts_offset = starts_offsets[u];
          if (u_starts_offset != UINT_E_MAX) { // non-zero number of updates for v in this batch
            auto u_inserts = get_new_edgelist(u, u_starts_offset);
            count += truncated_intersect(v, u_inserts);
          }
        }
        counts_two[b] = count;
      });

      pbbs::sequence<size_t> counts_three(batch.size());
      parallel_for(0, batch.size(), [&] (size_t b) {
        auto [u, v] = batch[b];
        size_t count = 0;
        if (u < v) {
          size_t v_starts_offset = starts_offsets[v];
          size_t u_starts_offset = starts_offsets[u];
          if (u_starts_offset != UINT_E_MAX && v_starts_offset != UINT_E_MAX) {
            auto u_inserts = get_new_edgelist(u, u_starts_offset);
            auto v_inserts = get_new_edgelist(v, v_starts_offset);
            count += intersect_new(u_inserts, v_inserts);
          }
        }
        counts_three[b] = count;
      });

      size_t first_count = pbbslib::reduce_add(counts_one.slice());
      size_t second_count = pbbslib::reduce_add(counts_two.slice());
      size_t third_count = pbbslib::reduce_add(counts_three.slice());

      size_t new_triangles = first_count - second_count + (third_count/3);
      T += new_triangles;
      std::cout << "first_count = " << first_count << " second_count = " << second_count << " third_count = " << third_count << std::endl;
      std::cout << "Number of new triangles:" << new_triangles << " total count = " << T << std::endl;

      // (cleanup) update starts_offsets mapping for batch vertices (reset at the end of this batch)
      parallel_for(0, starts.size(), [&] (size_t i) {
        const auto[v, index] = starts[i];
        starts_offsets[v] = UINT_E_MAX;
      });
    }


    template <class Container>
    size_t merge_deletions_inplace(uintE* A, uintE nA, Container& B) {
      size_t a = 0, b = 0, nB = B.size(), k = 0;
      while (a < nA && b < nB) {
        if (A[a] == B[b].second) { // deleted
          a++; b++;
        } else if (A[a] < B[b].second) { // not deleted
          A[k] = A[a];
          k++; a++;
        } else { // deletion not found
          b++;
        }
      }
      while (a < nA) {
        A[k] = A[a];
        k++; a++;
      }
      return k;
    }

    template <class C>
    void merge_deletions(uintE v, C& updates) {
      uintE current_degree = D[v];
      uintE* ngh_v = A[v];
      size_t new_degree = merge_deletions_inplace(ngh_v, current_degree, updates);
      D[v] = new_degree;
    }

    // An unsorted batch of updates
    template <class B>
    void process_deletions(B& unsorted_batch) {
      auto duplicated_batch = pbbs::sequence<edge>(2*unsorted_batch.size());
      parallel_for(0, unsorted_batch.size(), [&] (size_t i) {
        duplicated_batch[2*i] = unsorted_batch[i];
        duplicated_batch[2*i+1].from = unsorted_batch[i].to;
        duplicated_batch[2*i+1].to = unsorted_batch[i].from;
        duplicated_batch[2*i+1].weight = 1;
      });

      // (i) sort the batch
      auto sort_f = [&] (const edge& l, const edge& r) {
        return (l.from < r.from) || ((l.from == r.from) && (l.to < r.to));
      };
      timer sort_t; sort_t.start();
      pbbs::sample_sort_inplace(duplicated_batch.slice(), sort_f);


      // (ii) define the unweighted version (with ins/del) dropped and filter to
      // remove duplicates in the batch
      auto duplicated_batch_unweighted = pbbs::delayed_seq<std::pair<uintE, uintE>>(duplicated_batch.size(), [&] (const size_t i) {
        return std::make_pair(duplicated_batch[i].from, duplicated_batch[i].to);
      });
      auto batch = pbbs::filter_index(duplicated_batch_unweighted, [&] (const pair<uintE, uintE>& p, size_t ind) {
        return (ind == 0) || (p != duplicated_batch_unweighted[ind-1]);
      });

      // (iii) generate vertex offsets into the de-duplicated batch
      auto index_seq = pbbs::delayed_seq<std::pair<uintE, size_t>>(batch.size(), [&] (size_t i) {
        return std::make_pair(batch[i].first, i);
      });
      // starts contains (vertex id, indexof start in batch) pairs.
      auto starts = pbbs::filter(index_seq, [&] (const std::pair<uintE, size_t>& p) {
        size_t ind = p.second;
        return ind == 0 || (index_seq[ind].first != index_seq[ind-1].first);
      });
      std::cout << "batchsize = " << batch.size() << " starts.size = " << starts.size() << std::endl;

      // (iv) update starts_offsets mapping for batch vertices (reset at the end of this batch)
      parallel_for(0, starts.size(), [&] (size_t i) {
        const auto[v, index] = starts[i];
        starts_offsets[v] = i;
      });

      // (v) insert the updates into the old graph
      parallel_for(0, starts.size(), [&] (size_t i) {
        const auto [v, index] = starts[i];
        size_t v_deg = ((i == starts.size()-1) ? batch.size() : starts[i+1].second) - index;
        if (v_deg == 0) abort();

        auto v_inserts = batch.slice(index, index + v_deg);
        merge_deletions(v, v_inserts);
      });

      /////////////////////////////////////////////////////////////////////////
      // the next phase of the algorithm intersects the batch + the updated
      // graph, and updates the triangle counts

      // G(u) intersect G(v)
      pbbs::sequence<size_t> counts_one(batch.size());
      parallel_for(0, batch.size(), [&] (size_t b) {
        auto [u, v] = batch[b];
        if (u < v) {
          counts_one[b] = intersect_A(u, v);
        } else {
          counts_one[b] = 0;
        }
      });

      auto get_new_edgelist = [&] (const uintE& u, size_t u_starts_offset) {
        auto [up, index] = starts[u_starts_offset];
        assert(up == u);
        size_t u_deg = ((u_starts_offset == starts.size()-1) ? batch.size() : starts[u_starts_offset+1].second) - index;
        if (u_deg == 0) abort();
        auto u_inserts = batch.slice(index, index + u_deg);
        return u_inserts;
      };

      // truncated G(u) intersect G'(v) and truncated G(v) intersect G'(u)
      pbbs::sequence<size_t> counts_two(batch.size());
      parallel_for(0, batch.size(), [&] (size_t b) {
        auto [u, v] = batch[b];
        size_t count = 0;

        if (u < v) {
          size_t v_starts_offset = starts_offsets[v];
          if (v_starts_offset != UINT_E_MAX) { // non-zero number of updates for v in this batch
            auto v_inserts = get_new_edgelist(v, v_starts_offset);
            count += truncated_intersect(u, v_inserts);
          }

          size_t u_starts_offset = starts_offsets[u];
          if (u_starts_offset != UINT_E_MAX) { // non-zero number of updates for v in this batch
            auto u_inserts = get_new_edgelist(u, u_starts_offset);
            count += truncated_intersect(v, u_inserts);
          }
        }
        counts_two[b] = count;
      });

      pbbs::sequence<size_t> counts_three(batch.size());
      parallel_for(0, batch.size(), [&] (size_t b) {
        auto [u, v] = batch[b];
        size_t count = 0;
        if (u < v) {
          size_t v_starts_offset = starts_offsets[v];
          size_t u_starts_offset = starts_offsets[u];
          if (u_starts_offset != UINT_E_MAX && v_starts_offset != UINT_E_MAX) {
            auto u_inserts = get_new_edgelist(u, u_starts_offset);
            auto v_inserts = get_new_edgelist(v, v_starts_offset);
            count += intersect_new(u_inserts, v_inserts);
          }
        }
        counts_three[b] = count;
      });

      size_t first_count = pbbslib::reduce_add(counts_one.slice());
      size_t second_count = pbbslib::reduce_add(counts_two.slice());
      size_t third_count = pbbslib::reduce_add(counts_three.slice());

      size_t triangles_deleted = first_count + second_count + third_count;
      T -= triangles_deleted;
      std::cout << "first_count = " << first_count << " second_count = " << second_count << " third_count = " << third_count << std::endl;
      std::cout << "Number of new triangles deleted:" << triangles_deleted << " total count = " << T << std::endl;
    }

    void report_stats() {
      // compute reduction
       auto seq_copy = pbbs::sequence<size_t>(n, [&] (size_t i) { return D[i]; });
       size_t sum_deg = pbbslib::reduce_add(seq_copy.slice());
       std::cout << "Deg = " << sum_deg << std::endl;
    }

  };


}  // namespace gbbs
