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
#include <exception>
#include <fstream>
#include <string>
#include <sys/mman.h>
#include <type_traits>
#include <vector>

#include "gbbs/gbbs.h"
#include "gbbs/graph_io.h"

namespace gbbs {

template <class weight_type>
struct DynamicEdge {
  uintE from;
  uintE to;
  weight_type weight;
  bool insert;
  uintE batch_id;

  DynamicEdge() {}
  DynamicEdge(uintE _from) : from(_from), to(0) {}
  DynamicEdge(uintE _from, uintE _to, char _insert, uintE _batch_id)
    : from(_from)
    , to(_to)
    , insert(_insert == '+')
    , batch_id(_batch_id) {}
  DynamicEdge(const uintE _from, const uintE _to, const weight_type _weight, char _insert, uintE _batch_id)
    : from(_from)
    , to(_to)
    , weight(_weight)
    , insert(_insert == '+')
    , batch_id(_batch_id) {}
};

template <class weight_type>
struct BatchDynamicEdges {
  sequence<sequence<DynamicEdge<weight_type>>> edges;
  uintE max_vertex = 0;
};

// Read weighted edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <+/-> <edge 1 first endpoint> <edge 1 second endpoint> <edge 1 weight> <batch id>
//     <+/-> <edge 2 first endpoint> <edge 2 second endpoint> <edge 2 weight> <batch id>
//     <+/-> <edge 3 first endpoint> <edge 3 second endpoint> <edge 3 weight> <batch id>
//     ...
//     <+/-> <edge m first endpoint> <edge m second endpoint> <edge m weight> <batch id>
template <class weight_type>
std::vector<DynamicEdge<weight_type>> read_dynamic_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  gbbs_io::internal::skip_ifstream_comments(&file);

  std::vector<DynamicEdge<weight_type>> edge_list;
  char insert;
  uintE from;
  uintE to;
  weight_type weight;
  uintE batch_id;
  while (file >> insert >> from >> to >> weight >> batch_id) {
    edge_list.emplace_back(insert, from, to, weight, batch_id);
  }
  return edge_list;
}

// Read edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <+/-> <edge 1 first endpoint> <edge 1 second endpoint> <batch id>
//     <+/-> <edge 2 first endpoint> <edge 2 second endpoint> <batch id>
//     <+/-> <edge 3 first endpoint> <edge 3 second endpoint> <batch id>
//     ...
//     <+/-> <edge m first endpoint> <edge m second endpoint> <batch id>
std::vector<DynamicEdge<pbbslib::empty>>
read_dynamic_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  gbbs_io::internal::skip_ifstream_comments(&file);

  std::vector<DynamicEdge<pbbslib::empty>> edge_list;
  char insert;
  uintE from;
  uintE to;
  uintE batch_id;
  while (file >> insert >> from >> to >> batch_id) {
    edge_list.emplace_back(insert, from, to, batch_id);
  }
  return edge_list;
}

template <class weight_type>
std::vector<DynamicEdge<weight_type>> read_ordered_dynamic_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  gbbs_io::internal::skip_ifstream_comments(&file);

  std::vector<DynamicEdge<weight_type>> edge_list;
  char insert;
  uintE from;
  uintE to;
  weight_type weight;
  uintE batch_id = 0;
  char prev_insert = '\0';
  while (file >> insert >> from >> to >> weight) {
    if (insert != prev_insert) {
      prev_insert = insert;
      batch_id++;
    }
    edge_list.emplace_back(insert, from, to, weight, batch_id);
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return edge_list;
}

// This should work with the Sun et al. format
std::vector<DynamicEdge<pbbslib::empty>>
read_ordered_dynamic_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  gbbs_io::internal::skip_ifstream_comments(&file);

  std::vector<DynamicEdge<pbbslib::empty>> edge_list;
  char insert;
  uintE from;
  uintE to;
  uintE batch_id = 0;
  char prev_insert = '\0';
  while (file >> insert >> from >> to) {
    if (insert != prev_insert) {
      prev_insert = insert;
      batch_id++;
    }
    edge_list.emplace_back(insert, from, to, batch_id);
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return edge_list;
}

template <class A>
parlay::sequence<A> GetBoundaryIndices(
    std::size_t num_keys,
    const std::function<bool(std::size_t, std::size_t)>& key_eq_func) {
  parlay::sequence<A> mark_keys(num_keys + 1);
  auto null_key = std::numeric_limits<A>::max();
  parallel_for(0, num_keys, [&](std::size_t i) {
    if (i != 0 && key_eq_func(i, i - 1))
      mark_keys[i] = null_key;
    else
      mark_keys[i] = i;
  });
  mark_keys[num_keys] = num_keys; 
  return filter(mark_keys, [&null_key](A x) -> bool { return x != null_key; });
}

template <class weight_type>
parlay::sequence<parlay::sequence<DynamicEdge<weight_type>>> batch_dynamic_edge_list(
  std::vector<DynamicEdge<weight_type>>& edge_list, bool is_sorted = false) {
  using D = DynamicEdge<weight_type>;
  using batch_type = parlay::sequence<DynamicEdge<weight_type>>;
  if (edge_list.size() == 0) return parlay::sequence<batch_type>();
  // Sort in increasing order of batch_id
  batch_type sorted_edges;
  if (!is_sorted) {
    auto edge_list_f = [&](std::size_t i) {return edge_list[i];};
    auto edge_list_seq = parlay::delayed_seq<D>(edge_list.size(), edge_list_f);
    sorted_edges = parlay::internal::sample_sort(parlay::make_slice(edge_list_seq), [](const DynamicEdge<weight_type>& a, const DynamicEdge<weight_type>& b){
        return (a.batch_id < b.batch_id); // || (a.batch_id == b.batch_id && a.insert > b.insert)
      }, true);
  } else {
    sorted_edges = batch_type::from_function(edge_list.size(), [&](std::size_t i){ return edge_list[i]; });
  }

  // Extract ranges of equal batch_ids
  auto filtered_mark_ids = GetBoundaryIndices<uintE>(
    sorted_edges.size(), [&](std::size_t i, std::size_t j) {
      return sorted_edges[i].batch_id == sorted_edges[j].batch_id; //&& sorted_edges[i].insert == sorted_edges[j].insert
    });
  auto num_filtered_mark_ids = filtered_mark_ids.size() - 1;

  //auto sorted_edges_arr = sorted_edges.to_array();

  // Now rearrange sorted_edges into appropriate nested sequences
  parlay::sequence<batch_type> batch_edge_list = parlay::sequence<batch_type>::from_function(num_filtered_mark_ids, [&](std::size_t i){
    auto start_id_index = filtered_mark_ids[i];
    auto end_id_index = filtered_mark_ids[i + 1];
    return batch_type::from_function(end_id_index - start_id_index, [&](std::size_t j){
      return sorted_edges[start_id_index + j];
    });
    //return batch_type::from_function(sorted_edges_arr + start_id_index, end_id_index - start_id_index);
  });
  return batch_edge_list;
} 

template <class weight_type>
BatchDynamicEdges<weight_type> read_batch_dynamic_edge_list(const char* filename, bool ordered) {
  using D = DynamicEdge<weight_type>;
  std::vector<D> edge_list = ordered ? read_ordered_dynamic_edge_list(filename) : read_dynamic_edge_list(filename);
  auto edge_list_f = [&](std::size_t i) -> D{return edge_list[i];};
  auto edge_list_seq = parlay::delayed_seq<D>(edge_list.size(), edge_list_f);
  auto max_f = [](const D& l, const D& r) __attribute__((always_inline)) {
    return D(std::max({l.from, l.to, r.from, r.to}));
  };
  auto max_monoid = parlay::make_monoid(max_f, D(0));
  uintE max_vertex = reduce(edge_list_seq, max_monoid).from + 1;
  return BatchDynamicEdges<weight_type>{batch_dynamic_edge_list(edge_list, ordered), max_vertex};
}

// something that will read dynamic edge list, and write to graph form (output)

}  // namespace gbbs
