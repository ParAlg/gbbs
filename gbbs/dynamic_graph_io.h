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

  DynamicEdge() {}
  DynamicEdge(uintE _from) : from(_from), to(0) {}
  DynamicEdge(char _insert, uintE _from, uintE _to)
    : from(_from)
    , to(_to)
    , insert(_insert == '+') {}
  DynamicEdge(char _insert, const uintE _from, const uintE _to, const weight_type _weight)
    : from(_from)
    , to(_to)
    , weight(_weight)
    , insert(_insert == '+') {}
};

template <class weight_type>
struct BatchDynamicEdges {
  std::vector<DynamicEdge<weight_type>> edges;
  uintE max_vertex = 0;
};

// Read weighted edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <+/-> <edge 1 first endpoint> <edge 1 second endpoint> <edge 1 weight>
//     <+/-> <edge 2 first endpoint> <edge 2 second endpoint> <edge 2 weight>
//     <+/-> <edge 3 first endpoint> <edge 3 second endpoint> <edge 3 weight>
//     ...
//     <+/-> <edge m first endpoint> <edge m second endpoint> <edge m weight>
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
  while (file >> insert >> from >> to >> weight) {
    edge_list.emplace_back(insert, from, to, weight);
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return edge_list;
}

// Read edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <+/-> <edge 1 first endpoint> <edge 1 second endpoint>
//     <+/-> <edge 2 first endpoint> <edge 2 second endpoint>
//     <+/-> <edge 3 first endpoint> <edge 3 second endpoint>
//     ...
//     <+/-> <edge m first endpoint> <edge m second endpoint>
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
  while (file >> insert >> from >> to) {
    edge_list.emplace_back(insert, from, to);
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return edge_list;
}

template <class weight_type>
BatchDynamicEdges<weight_type> read_batch_dynamic_edge_list(const char* filename) {
  using D = DynamicEdge<weight_type>;
  std::vector<D> edge_list = read_dynamic_edge_list(filename);

  // Find the max vertex
  auto edge_list_f = [&](size_t i) -> D {return edge_list[i];};
  auto edge_list_seq = parlay::delayed_seq<D>(edge_list.size(), edge_list_f);
  auto max_f = [](const D& l, const D& r) __attribute__((always_inline)) {
    return D(std::max({l.from, l.to, r.from, r.to}));
  };
  auto max_monoid = parlay::make_monoid(max_f, D(0));
  uintE max_vertex = parlay::reduce(edge_list_seq, max_monoid).from + 1;

  return BatchDynamicEdges<weight_type>{std::move(edge_list), max_vertex};
}

// something that will read dynamic edge list, and write to graph form (output)

}  // namespace gbbs
