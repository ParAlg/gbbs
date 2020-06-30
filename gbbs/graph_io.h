// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
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

#include <algorithm>
#include <exception>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

#include "gbbs/graph.h"
#include "gbbs/io.h"
#include "gbbs/macros.h"
#include "gbbs/vertex.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/stlalgs.h"

namespace gbbs {
namespace gbbs_io {

template <class weight_type>
struct Edge {
  uintE from;
  uintE to;
  weight_type weight;

  Edge() {}
  Edge(uintE _from, uintE _to)
    : from(_from)
    , to(_to)
    , weight(0) {}
  Edge(const uintE _from, const uintE _to, const weight_type _weight)
    : from(_from)
    , to(_to)
    , weight(_weight) {}
};

template <>
Edge<pbbslib::empty>::Edge(const uintE _from, const uintE _to);

namespace internal {  // Internal declarations

// Header string expected at the top of unweighted adjacency graph files.
const std::string kUnweightedAdjGraphHeader = "AdjacencyGraph";
// Header string expected at the top of weighted adjacency graph files.
const std::string kWeightedAdjGraphHeader = "WeightedAdjacencyGraph";

template <class weight_type>
size_t get_num_vertices_from_edges(const pbbs::sequence<Edge<weight_type>>&);

template <class weight_type>
pbbs::sequence<vertex_data> sorted_edges_to_vertex_data_array(
    size_t, const pbbs::sequence<Edge<weight_type>>&);

} // namespace internal

/* Returns a tuple containing (n, m, offsets, edges) --- the number of
 * vertices, edges, the vertex offsets, and the edge values, after
 * parsing the input (weighted) graph file */
std::tuple<size_t, size_t, uintT*, std::tuple<uintE, intE>*>
parse_weighted_graph(
  const char* fname,
  bool mmap,
  char* bytes = nullptr,
  size_t bytes_size = std::numeric_limits<size_t>::max());

symmetric_graph<symmetric_vertex, intE> read_weighted_symmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

asymmetric_graph<asymmetric_vertex, intE> read_weighted_asymmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

/* Returns a tuple containing (n, m, offsets, edges) --- the number of
 * vertices, edges, the vertex offsets, and the edge values, after
 * parsing the input graph file */
std::tuple<size_t, size_t, uintT*, uintE*> parse_unweighted_graph(
    const char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

symmetric_graph<symmetric_vertex, pbbslib::empty> read_unweighted_symmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

asymmetric_graph<asymmetric_vertex, pbbslib::empty> read_unweighted_asymmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

std::tuple<char*, size_t> parse_compressed_graph(
    const char* fname, bool mmap, bool mmapcopy);

template <class weight_type>
symmetric_graph<csv_bytepd_amortized, weight_type>
read_compressed_symmetric_graph(const char* fname, bool mmap, bool mmapcopy) {
  char* bytes;
  size_t bytes_size;
  std::tie(bytes, bytes_size) = parse_compressed_graph(fname, mmap, mmapcopy);

  long* sizes = (long*)bytes;
  uint64_t n = sizes[0], m = sizes[1];

  debug(uint64_t totalSpace = sizes[2];
  std::cout << "# n = " << n << " m = " << m << " totalSpace = " << totalSpace
            << "\n");

  uintT* offsets = (uintT*)(bytes + 3 * sizeof(long));
  uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees = (uintE*)(bytes + skip);
  skip += n * sizeof(intE);
  uchar* edges = (uchar*)(bytes + skip);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = Degrees[i];
  });

  std::function<void()> deletion_fn = [=] () { pbbslib::free_arrays(v_data, bytes); };
  if (mmap && !mmapcopy) {
    deletion_fn = [v_data, bytes, bytes_size] () {
      pbbslib::free_array(v_data);
      unmmap(bytes, bytes_size);
    };
  }
  symmetric_graph<csv_bytepd_amortized, weight_type> G(v_data, n, m, deletion_fn, edges);
  return G;
}

template <class weight_type>
asymmetric_graph<cav_bytepd_amortized, weight_type>
read_compressed_asymmetric_graph(const char* fname, bool mmap, bool mmapcopy) {
  char* bytes;
  size_t bytes_size;
  std::tie(bytes, bytes_size) = parse_compressed_graph(fname, mmap, mmapcopy);

  long* sizes = (long*)bytes;
  uint64_t n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  debug(std::cout << "# n = " << n << " m = " << m << " totalSpace = " << totalSpace
            << "\n");

  uintT* offsets = (uintT*)(bytes + 3 * sizeof(long));
  uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees = (uintE*)(bytes + skip);
  skip += n * sizeof(intE);
  uchar* edges = (uchar*)(bytes + skip);

  uintT* inOffsets;
  uchar* inEdges;
  uintE* inDegrees;

  skip += totalSpace;
  uchar* inData = (uchar*)(bytes + skip);
  sizes = (long*)inData;
  debug(size_t inTotalSpace = sizes[0];
  std::cout << "# inTotalSpace = " << inTotalSpace << "\n";);
  skip += sizeof(long);
  inOffsets = (uintT*)(bytes + skip);
  skip += (n + 1) * sizeof(uintT);
  inDegrees = (uintE*)(bytes + skip);
  skip += n * sizeof(uintE);
  inEdges = (uchar*)(bytes+ skip);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  auto v_in_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = Degrees[i];

    v_in_data[i].offset = inOffsets[i];
    v_in_data[i].degree = inDegrees[i];
  });

  std::function<void()> deletion_fn = [=] () { pbbslib::free_arrays(v_data, v_in_data, bytes); };
  if (mmap && !mmapcopy) {
    deletion_fn = [v_data, v_in_data, bytes, bytes_size] () {
      pbbslib::free_array(v_data);
      pbbslib::free_array(v_in_data);
      unmmap(bytes, bytes_size);
    };
  }

  asymmetric_graph<cav_bytepd_amortized, weight_type> G(v_data, v_in_data, n, m, deletion_fn, edges, inEdges);
  return G;
}

// Read weighted edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <edge 1 first endpoint> <edge 1 second endpoint> <edge 1 weight>
//     <edge 2 first endpoint> <edge 2 second endpoint> <edge 2 weight>
//     <edge 3 first endpoint> <edge 3 second endpoint> <edge 3 weight>
//     ...
//     <edge m first endpoint> <edge m second endpoint> <edge m weight>
std::vector<Edge<intT>> read_weighted_edge_list(const char* filename);

// Read edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <edge 1 first endpoint> <edge 1 second endpoint>
//     <edge 2 first endpoint> <edge 2 second endpoint>
//     <edge 3 first endpoint> <edge 3 second endpoint>
//     ...
//     <edge m first endpoint> <edge m second endpoint>
std::vector<Edge<pbbslib::empty>>
read_unweighted_edge_list(const char* filename);

// Converts edge list into an asymmetric graph.
//
// Adjacency lists of the output graph are sorted by increasing neighbor vertex
// ID.
//
// Duplicate edges are removed. If there are multiple edges between the same
// endpoints with different weights, an arbitrary one is kept.
template <class weight_type>
asymmetric_graph<asymmetric_vertex, weight_type>
edge_list_to_asymmetric_graph(const std::vector<Edge<weight_type>>& edge_list) {
  using edge_type = typename asymmetric_vertex<weight_type>::edge_type;

  if (edge_list.empty()) {
    return asymmetric_graph<asymmetric_vertex, weight_type>{};
  }

  constexpr auto compare_endpoints = [](
      const Edge<weight_type>& left,
      const Edge<weight_type>& right) {
    return std::tie(left.from, left.to) < std::tie(right.from, right.to);
  };
  const auto edge_list_seq = pbbs::delayed_seq<Edge<weight_type>>(
      edge_list.size(), [&](const size_t i) { return edge_list[i]; });
  pbbs::sequence<Edge<weight_type>> out_edges =
    pbbs::remove_duplicates_ordered(edge_list_seq, compare_endpoints);
  const size_t num_edges = out_edges.size();
  const size_t num_vertices = internal::get_num_vertices_from_edges(out_edges);
  pbbs::sequence<vertex_data> vertex_out_data =
    internal::sorted_edges_to_vertex_data_array(num_vertices, out_edges);

  pbbs::sequence<Edge<weight_type>> in_edges =
    pbbs::map<Edge<weight_type>>(out_edges, [&](const Edge<weight_type>& edge) {
      return Edge<weight_type>{edge.to, edge.from, edge.weight};
    });
  pbbs::sample_sort_inplace(in_edges.slice(), compare_endpoints);
  pbbs::sequence<vertex_data> vertex_in_data =
    internal::sorted_edges_to_vertex_data_array(num_vertices, in_edges);

  edge_type* out_edges_array = pbbs::new_array_no_init<edge_type>(num_edges);
  edge_type* in_edges_array = pbbs::new_array_no_init<edge_type>(num_edges);
  par_for(0, num_edges, [&](const size_t i) {
    const Edge<weight_type>& out_edge = out_edges[i];
    out_edges_array[i] = std::make_tuple(out_edge.to, out_edge.weight);
    const Edge<weight_type>& in_edge = in_edges[i];
    in_edges_array[i] = std::make_tuple(in_edge.to, in_edge.weight);
  });

  auto vertex_out_data_array = vertex_out_data.to_array();
  auto vertex_in_data_array = vertex_in_data.to_array();
  return asymmetric_graph<asymmetric_vertex, weight_type>{
    vertex_out_data_array,
    vertex_in_data_array,
    num_vertices,
    num_edges,
    [=] () { pbbslib::free_arrays(vertex_out_data_array, vertex_in_data_array, out_edges_array, in_edges_array); },
    out_edges_array,
    in_edges_array};
}

// Converts a list of undirected edges into a symmetric graph.
//
// Adjacency lists of the output graph are sorted by increasing neighbor vertex
// ID.
//
// Duplicate edges are removed. If there are multiple edges between the same
// endpoints with different weights, an arbitrary one is kept.
template <class weight_type>
symmetric_graph<symmetric_vertex, weight_type>
edge_list_to_symmetric_graph(const std::vector<Edge<weight_type>>& edge_list) {
  using edge_type = typename symmetric_vertex<weight_type>::edge_type;

  if (edge_list.empty()) {
    return symmetric_graph<symmetric_vertex, weight_type>{};
  }

  pbbs::sequence<Edge<weight_type>> edges_both_directions(2 * edge_list.size());
  par_for(0, edge_list.size(), [&](const size_t i) {
      const Edge<weight_type>& edge = edge_list[i];
      edges_both_directions[2 * i] = edge;
      edges_both_directions[2 * i + 1] =
        Edge<weight_type>{edge.to, edge.from, edge.weight};
  });
  constexpr auto compare_endpoints = [](
      const Edge<weight_type>& left,
      const Edge<weight_type>& right) {
    return std::tie(left.from, left.to) < std::tie(right.from, right.to);
  };
  pbbs::sequence<Edge<weight_type>> edges =
    pbbs::remove_duplicates_ordered(edges_both_directions, compare_endpoints);
  const size_t num_edges = edges.size();
  const size_t num_vertices = internal::get_num_vertices_from_edges(edges);
  pbbs::sequence<vertex_data> vertex_data =
    internal::sorted_edges_to_vertex_data_array(num_vertices, edges);

  edge_type* edges_array = pbbs::new_array_no_init<edge_type>(num_edges);
  par_for(0, num_edges, [&](const size_t i) {
    const Edge<weight_type>& edge = edges[i];
    edges_array[i] = std::make_tuple(edge.to, edge.weight);
  });

  auto vertex_data_array = vertex_data.to_array();
  return symmetric_graph<symmetric_vertex, weight_type>{
    vertex_data_array,
    num_vertices,
    num_edges,
    [=] () { pbbslib::free_arrays(vertex_data_array, edges_array); },
    edges_array};
}

// Write graph in adjacency graph format to file.
template <class Graph>
void write_graph_to_file(const char* filename, Graph& graph) {
  using weight_type = typename Graph::weight_type;

  std::ofstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }

  constexpr bool is_weighted_graph{
    !std::is_same<weight_type, pbbslib::empty>::value};
  const size_t num_vertices{graph.n};
  const size_t num_edges{graph.m};

  pbbs::sequence<size_t> offsets{num_vertices,
    [&](const size_t i) {
      return graph.get_vertex(i).getOutDegree();
    }};
  pbbslib::scan_add_inplace(offsets);

  file << (is_weighted_graph
      ? internal::kWeightedAdjGraphHeader
      : internal::kUnweightedAdjGraphHeader) << '\n';
  file << num_vertices << '\n';
  file << num_edges  << '\n';
  for (size_t i = 0; i < num_vertices; i++) {
    file << offsets[i] << '\n';
  }
  for (size_t i = 0; i < num_vertices; i++) {
    constexpr bool kParallel{false};
    const auto print_neighbor{[&](uintE, const uintE neighbor_id, weight_type) {
      file << neighbor_id << '\n';
    }};
    graph.get_vertex(i).mapOutNgh(i, print_neighbor, kParallel);
  }
  if constexpr (is_weighted_graph) {
    for (size_t i = 0; i < num_vertices; i++) {
      constexpr bool kParallel{false};
      const auto print_weight{[&](uintE, uintE, const weight_type& weight) {
        file << weight << '\n';
      }};
      graph.get_vertex(i).mapOutNgh(i, print_weight, kParallel);
    }
  }
}

namespace internal {  // Internal definitions

// Given a list of edges on a graph with vertex IDs {0, 1, 2, 3,..., n - 1},
// return the minimal valid value of n.
template <class weight_type>
size_t
get_num_vertices_from_edges(const pbbs::sequence<Edge<weight_type>>& edges) {
  const auto max_endpoints = pbbs::delayed_seq<size_t>(
    edges.size(),
    [&](const size_t i) {
      return std::max(edges[i].from, edges[i].to);
    });
    return pbbslib::reduce_max(max_endpoints) + 1;
}

// Given a list of edges sorted by their first endpoint, return a corresponding
// vertex_data array.
template <class weight_type>
pbbs::sequence<vertex_data> sorted_edges_to_vertex_data_array(
    const size_t num_vertices,
    const pbbs::sequence<Edge<weight_type>>& edges) {
  if (edges.empty()) {
    return {};
  }

  const size_t num_edges = edges.size();
  pbbs::sequence<vertex_data> data(num_vertices);
  par_for(0, edges[0].from + 1, [&](const size_t j) {
      data[j].offset = 0;
  });
  par_for(1, num_edges, [&](const size_t i) {
    const size_t from = edges[i].from;
    const size_t previous_from = edges[i - 1].from;
    if (from != previous_from) {
      par_for(previous_from + 1, from + 1, [&](const size_t j) {
        data[j].offset = i;
      });
    }
  });
  par_for(edges[num_edges - 1].from + 1, num_vertices, [&](const size_t j) {
    data[j].offset = num_edges;
  });

  par_for(0, num_vertices, [&](const size_t j) {
    const size_t next_offset =
      j == num_vertices - 1 ? num_edges : data[j + 1].offset;
    data[j].degree = next_offset - data[j].offset;
  });

  return data;
}

}  // namespace internal

}  // namespace gbbs_io
}  // namespace gbbs
