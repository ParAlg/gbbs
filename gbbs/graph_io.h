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

#include <sys/mman.h>
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

#include "parlay/io.h"

namespace gbbs {
namespace gbbs_io {

template <class weight_type>
struct Edge {
  uintE from;
  uintE to;
  weight_type weight;

  Edge() {}
  Edge(uintE _from, uintE _to) : from(_from), to(_to), weight(0) {}
  Edge(const uintE _from, const uintE _to, const weight_type _weight)
      : from(_from), to(_to), weight(_weight) {}
};

template <>
Edge<gbbs::empty>::Edge(uintE _from, uintE _to);

namespace internal {  // Internal declarations

// Header string expected at the top of unweighted adjacency graph files.
const std::string kUnweightedAdjGraphHeader = "AdjacencyGraph";
// Header string expected at the top of weighted adjacency graph files.
const std::string kWeightedAdjGraphHeader = "WeightedAdjacencyGraph";

void skip_ifstream_comments(std::ifstream* stream);

template <class weight_type>
size_t get_num_vertices_from_edges(const sequence<Edge<weight_type>>&);

template <class weight_type>
vertex_data* sorted_edges_to_vertex_data_array(
    size_t, const sequence<Edge<weight_type>>&);

template <class weight_type>
std::tuple<size_t, size_t, uintT*, std::tuple<uintE, weight_type>*>
parse_weighted_graph(const char* fname, bool mmap, bool binary,
                     char* bytes = nullptr,
                     size_t bytes_size = std::numeric_limits<size_t>::max());

// Output a list of sorted edges with no duplicates and no self-loop edges.  If
// there are multiple edges between the same endpoints with different weights,
// an arbitrary one is kept.
//
// The argument is passed by value, so use `std::move` where appropriate.
template <class weight_type>
sequence<Edge<weight_type>> sort_and_dedupe(sequence<Edge<weight_type>> edges);

}  // namespace internal

/* Returns a tuple containing (n, m, offsets, edges) --- the number of
 * vertices, edges, the vertex offsets, and the edge values, after
 * parsing the input graph file */
std::tuple<size_t, size_t, uintT*, uintE*> parse_unweighted_graph(
    const char* fname, bool mmap, bool binary, char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

symmetric_graph<symmetric_vertex, gbbs::empty> read_unweighted_symmetric_graph(
    const char* fname, bool mmap, bool binary, char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

asymmetric_graph<asymmetric_vertex, gbbs::empty>
read_unweighted_asymmetric_graph(
    const char* fname, bool mmap, bool binary, char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

std::tuple<char*, size_t> parse_compressed_graph(const char* fname, bool mmap);

template <class weight_type>
symmetric_graph<symmetric_vertex, weight_type> read_weighted_symmetric_graph(
    const char* fname, bool mmap, bool binary, char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max()) {
  size_t n, m;
  uintT* offsets;
  std::tuple<uintE, weight_type>* edges;
  std::tie(n, m, offsets, edges) = internal::parse_weighted_graph<weight_type>(
      fname, mmap, binary, bytes, bytes_size);

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i + 1] - v_data[i].offset;
  });
  if (!binary) {
    gbbs::free_array(offsets, n + 1);
  }

  return symmetric_graph<symmetric_vertex, weight_type>(
      v_data, n, m,
      [=]() {
        gbbs::free_array(v_data, n);
        if (!binary) {
          gbbs::free_array(edges, 2 * m);
        }
      },
      edges);
}

template <class weight_type>
asymmetric_graph<asymmetric_vertex, weight_type> read_weighted_asymmetric_graph(
    const char* fname, bool mmap, bool binary, char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max()) {
  using id_and_weight = std::tuple<uintE, weight_type>;
  using triple = std::pair<uintE, std::pair<uintE, weight_type>>;

  size_t n, m;
  uintT* offsets;
  std::tuple<uintE, weight_type>* edges;
  if (binary) {
    std::cout << "Todo: implement binary support for asymmetric graphs"
              << std::endl;
    exit(-1);
  }
  std::tie(n, m, offsets, edges) = internal::parse_weighted_graph<weight_type>(
      fname, mmap, binary, bytes, bytes_size);

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i + 1] - v_data[i].offset;
  });
  gbbs::free_array(offsets, n + 1);

  auto tOffsets = sequence<uintT>::uninitialized(n + 1);
  parallel_for(0, n, [&](size_t i) { tOffsets[i] = INT_T_MAX; });
  triple* temp = gbbs::new_array_no_init<triple>(m);
  parallel_for(0, n, [&](size_t i) {
    uintT o = v_data[i].offset;
    uintE deg = v_data[i].degree;
    for (uintT j = 0; j < deg; j++) {
      auto& cur_edge = (edges + o)[j];
      temp[o + j] =
          std::make_pair(std::get<0>(cur_edge),
                         std::make_pair((uintE)i, std::get<1>(cur_edge)));
    }
  });

  auto temp_seq = gbbs::make_slice(temp, m);
  parlay::integer_sort_inplace(temp_seq,
                               [&](const triple& p) { return p.first; });

  tOffsets[temp[0].first] = 0;
  id_and_weight* inEdges = gbbs::new_array_no_init<id_and_weight>(m);
  inEdges[0] = std::make_tuple(temp[0].second.first, temp[0].second.second);

  parallel_for(1, m, [&](size_t i) {
    inEdges[i] = std::make_tuple(temp[i].second.first, temp[i].second.second);
    if (temp[i].first != temp[i - 1].first) {
      tOffsets[temp[i].first] = i;
    }
  });

  gbbs::free_array(temp, m);

  // fill in offsets of degree 0 vertices by taking closest non-zero
  // offset to the right
  auto t_seq = parlay::make_slice(tOffsets.rbegin(), tOffsets.rend());
  auto M = parlay::minimum<uintT>();
  M.identity = m;
  parlay::scan_inclusive_inplace(t_seq, M);

  auto v_in_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    v_in_data[i].offset = tOffsets[i];
    v_in_data[i].degree = tOffsets[i + 1] - v_in_data[i].offset;
  });

  return asymmetric_graph<asymmetric_vertex, weight_type>(
      v_data, v_in_data, n, m,
      [=]() {
        gbbs::free_array(v_data, n);
        gbbs::free_array(v_in_data, n);
        gbbs::free_array(edges, m);
        gbbs::free_array(inEdges, m);
      },
      edges, inEdges);
}

template <class weight_type>
symmetric_graph<csv_bytepd_amortized, weight_type>
read_compressed_symmetric_graph(const char* fname, bool mmap) {
  char* bytes;
  size_t bytes_size;
  std::tie(bytes, bytes_size) = parse_compressed_graph(fname, mmap);

  long* sizes = (long*)bytes;
  uint64_t n = sizes[0], m = sizes[1];

  gbbs_debug(uint64_t totalSpace = sizes[2];
        std::cout << "# n = " << n << " m = " << m
                  << " totalSpace = " << totalSpace << "\n");

  uintT* offsets = (uintT*)(bytes + 3 * sizeof(long));
  uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees = (uintE*)(bytes + skip);
  skip += n * sizeof(intE);
  uchar* edges = (uchar*)(bytes + skip);

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = Degrees[i];
  });

  std::function<void()> deletion_fn = [=]() {
    gbbs::free_array(v_data, n);
    gbbs::free_array(bytes, bytes_size);
  };

  if (mmap) {
    deletion_fn = [=]() {
      gbbs::free_array(v_data, n);
      unmmap(bytes, bytes_size);
    };
  }
  symmetric_graph<csv_bytepd_amortized, weight_type> G(v_data, n, m,
                                                       std::move(deletion_fn), edges);
  return G;
}

template <class weight_type>
asymmetric_graph<cav_bytepd_amortized, weight_type>
read_compressed_asymmetric_graph(const char* fname, bool mmap) {
  char* bytes;
  size_t bytes_size;
  std::tie(bytes, bytes_size) = parse_compressed_graph(fname, mmap);

  long* sizes = (long*)bytes;
  uint64_t n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  gbbs_debug(std::cout << "# n = " << n << " m = " << m
                  << " totalSpace = " << totalSpace << "\n");

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
  gbbs_debug(size_t inTotalSpace = sizes[0];
        std::cout << "# inTotalSpace = " << inTotalSpace << "\n";);
  skip += sizeof(long);
  inOffsets = (uintT*)(bytes + skip);
  skip += (n + 1) * sizeof(uintT);
  inDegrees = (uintE*)(bytes + skip);
  skip += n * sizeof(uintE);
  inEdges = (uchar*)(bytes + skip);

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  auto v_in_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = Degrees[i];

    v_in_data[i].offset = inOffsets[i];
    v_in_data[i].degree = inDegrees[i];
  });

  std::function<void()> deletion_fn = [=]() {
    gbbs::free_array(v_data, n);
    gbbs::free_array(v_in_data, n);
    gbbs::free_array(bytes, bytes_size);
  };
  if (mmap) {
    deletion_fn = [=]() {
      gbbs::free_array(v_data, n);
      gbbs::free_array(v_in_data, n);
      unmmap(bytes, bytes_size);
    };
  }

  asymmetric_graph<cav_bytepd_amortized, weight_type> G(
      v_data, v_in_data, n, m, deletion_fn, edges, inEdges);
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
template <class weight_type>
std::vector<Edge<weight_type>> read_weighted_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  internal::skip_ifstream_comments(&file);

  std::vector<Edge<weight_type>> edge_list;
  uintE from;
  uintE to;
  weight_type weight;
  while (file >> from >> to >> weight) {
    edge_list.emplace_back(from, to, weight);
  }
  return edge_list;
}

// Read edges from a file that has the following format:
//     # There can be comments at the top of the file as long as each line of
//     # the comment starts with '#'.
//     <edge 1 first endpoint> <edge 1 second endpoint>
//     <edge 2 first endpoint> <edge 2 second endpoint>
//     <edge 3 first endpoint> <edge 3 second endpoint>
//     ...
//     <edge m first endpoint> <edge m second endpoint>
std::vector<Edge<gbbs::empty>> read_unweighted_edge_list(const char* filename);

// Converts edge list into an asymmetric graph.
//
// Adjacency lists of the output graph are sorted by increasing neighbor vertex
// ID.
//
// Duplicate edges and self-loop edges are removed. If there are multiple edges
// between the same endpoints with different weights, an arbitrary one is kept.
template <class weight_type>
asymmetric_graph<asymmetric_vertex, weight_type> edge_list_to_asymmetric_graph(
    const std::vector<Edge<weight_type>>& edge_list) {
  using edge_type = typename asymmetric_vertex<weight_type>::edge_type;

  if (edge_list.empty()) {
    return asymmetric_graph<asymmetric_vertex, weight_type>{};
  }

  const sequence<Edge<weight_type>> out_edges =
      internal::sort_and_dedupe(sequence<Edge<weight_type>>::from_function(
          edge_list.size(), [&](const size_t i) { return edge_list[i]; }));
  const size_t num_edges = out_edges.size();
  const size_t num_vertices = internal::get_num_vertices_from_edges(out_edges);
  vertex_data* vertex_out_data =
      internal::sorted_edges_to_vertex_data_array(num_vertices, out_edges);

  sequence<Edge<weight_type>> in_edges = parlay::map<Edge<weight_type>>(
      out_edges, [&](const Edge<weight_type>& edge) {
        return Edge<weight_type>{edge.to, edge.from, edge.weight};
      });
  constexpr auto compare_endpoints = [](const Edge<weight_type>& left,
                                        const Edge<weight_type>& right) {
    return std::tie(left.from, left.to) < std::tie(right.from, right.to);
  };
  parlay::sample_sort_inplace(make_slice(in_edges), compare_endpoints);
  vertex_data* vertex_in_data =
      internal::sorted_edges_to_vertex_data_array(num_vertices, in_edges);

  edge_type* out_edges_array = gbbs::new_array_no_init<edge_type>(num_edges);
  edge_type* in_edges_array = gbbs::new_array_no_init<edge_type>(num_edges);
  parallel_for(0, num_edges, [&](const size_t i) {
    const Edge<weight_type>& out_edge = out_edges[i];
    out_edges_array[i] = std::make_tuple(out_edge.to, out_edge.weight);
    const Edge<weight_type>& in_edge = in_edges[i];
    in_edges_array[i] = std::make_tuple(in_edge.to, in_edge.weight);
  });

  return asymmetric_graph<asymmetric_vertex, weight_type>{
      vertex_out_data, vertex_in_data, num_vertices, num_edges,
      [=]() {
        gbbs::free_array(vertex_out_data, num_vertices);
        gbbs::free_array(vertex_in_data, num_vertices);
        gbbs::free_array(out_edges_array, num_edges);
        gbbs::free_array(in_edges_array, num_edges);
      },
      out_edges_array, in_edges_array};
}

// Converts a list of undirected edges into a symmetric graph.
//
// Adjacency lists of the output graph are sorted by increasing neighbor vertex
// ID.
//
// Duplicate edges and self-loop edges are removed. If there are multiple edges
// between the same endpoints with different weights, an arbitrary one is kept.
template <class weight_type>
symmetric_graph<symmetric_vertex, weight_type> edge_list_to_symmetric_graph(
    const std::vector<Edge<weight_type>>& edge_list) {
  using edge_type = typename symmetric_vertex<weight_type>::edge_type;

  if (edge_list.empty()) {
    return symmetric_graph<symmetric_vertex, weight_type>{};
  }

  sequence<Edge<weight_type>> edges_both_directions(2 * edge_list.size());
  parallel_for(0, edge_list.size(), [&](const size_t i) {
    const Edge<weight_type>& edge = edge_list[i];
    edges_both_directions[2 * i] = edge;
    edges_both_directions[2 * i + 1] =
        Edge<weight_type>{edge.to, edge.from, edge.weight};
  });
  const sequence<Edge<weight_type>> edges =
      internal::sort_and_dedupe(std::move(edges_both_directions));
  const size_t num_edges = edges.size();
  const size_t num_vertices = internal::get_num_vertices_from_edges(edges);
  vertex_data* vertex_data =
      internal::sorted_edges_to_vertex_data_array(num_vertices, edges);

  edge_type* edges_array = gbbs::new_array_no_init<edge_type>(num_edges);
  parallel_for(0, num_edges, [&](const size_t i) {
    const Edge<weight_type>& edge = edges[i];
    edges_array[i] = std::make_tuple(edge.to, edge.weight);
  });

  return symmetric_graph<symmetric_vertex, weight_type>{
      vertex_data, num_vertices, num_edges,
      [=]() {
        gbbs::free_array(vertex_data, num_vertices);
        gbbs::free_array(edges_array, num_edges);
      },
      edges_array};
}

template <class weight_type>
symmetric_graph<symmetric_vertex, weight_type> edge_list_to_symmetric_graph(
    const edge_array<weight_type>& edge_list) {
  using edge_type = typename symmetric_vertex<weight_type>::edge_type;
  using Edge = gbbs_io::Edge<weight_type>;

  if (edge_list.E.size() == 0) {
    return symmetric_graph<symmetric_vertex, weight_type>{};
  }

  sequence<Edge> edges_both_directions(2 * edge_list.E.size());
  parallel_for(0, edge_list.E.size(), [&](const size_t i) {
    const auto& orig_edge = edge_list.E[i];
    Edge edge(std::get<0>(orig_edge), std::get<1>(orig_edge),
              std::get<2>(orig_edge));
    edges_both_directions[2 * i] = edge;
    edges_both_directions[2 * i + 1] = Edge{edge.to, edge.from, edge.weight};
  });
  const sequence<Edge> edges =
      internal::sort_and_dedupe(std::move(edges_both_directions));
  const size_t num_edges = edges.size();
  const size_t num_vertices = internal::get_num_vertices_from_edges(edges);
  vertex_data* vertex_data =
      internal::sorted_edges_to_vertex_data_array(num_vertices, edges);

  edge_type* edges_array = gbbs::new_array_no_init<edge_type>(num_edges);
  parallel_for(0, num_edges, [&](const size_t i) {
    const Edge& edge = edges[i];
    edges_array[i] = std::make_tuple(edge.to, edge.weight);
  });

  return symmetric_graph<symmetric_vertex, weight_type>{
      vertex_data, num_vertices, num_edges,
      [=]() {
        gbbs::free_array(vertex_data, num_vertices);
        gbbs::free_array(edges_array, num_edges);
      },
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
      !std::is_same<weight_type, gbbs::empty>::value};
  const size_t num_vertices{graph.n};
  const size_t num_edges{graph.m};

  sequence<size_t> offsets = sequence<size_t>::from_function(
      num_vertices,
      [&](const size_t i) { return graph.get_vertex(i).out_degree(); });
  parlay::scan_inplace(make_slice(offsets));

  file << (is_weighted_graph ? internal::kWeightedAdjGraphHeader
                             : internal::kUnweightedAdjGraphHeader)
       << '\n';
  file << num_vertices << '\n';
  file << num_edges << '\n';
  for (size_t i = 0; i < num_vertices; i++) {
    file << offsets[i] << '\n';
  }
  for (size_t i = 0; i < num_vertices; i++) {
    constexpr bool kParallel{false};
    const auto print_neighbor{[&](uintE, const uintE neighbor_id, weight_type) {
      file << neighbor_id << '\n';
    }};
    graph.get_vertex(i).out_neighbors().map(print_neighbor, kParallel);
  }
  if
    constexpr(is_weighted_graph) {
      for (size_t i = 0; i < num_vertices; i++) {
        constexpr bool kParallel{false};
        const auto print_weight{[&](uintE, uintE, const weight_type& weight) {
          file << weight << '\n';
        }};
        graph.get_vertex(i).out_neighbors().map(print_weight, kParallel);
      }
    }
}

namespace internal {  // Internal definitions

// For use in `static_assert(false)`. See
// https://stackoverflow.com/a/53945549/4865149 .
template <class...>
constexpr std::false_type always_false{};

// Given a list of edges on a graph with vertex IDs {0, 1, 2, 3,..., n - 1},
// return the minimal valid value of n.
template <class weight_type>
size_t get_num_vertices_from_edges(const sequence<Edge<weight_type>>& edges) {
  const auto max_endpoints = parlay::delayed_seq<size_t>(
      edges.size(),
      [&](const size_t i) { return std::max(edges[i].from, edges[i].to); });
  return parlay::reduce_max(max_endpoints) + 1;
}

// Given a list of edges sorted by their first endpoint, return a corresponding
// vertex_data array.
template <class weight_type>
vertex_data* sorted_edges_to_vertex_data_array(
    const size_t num_vertices, const sequence<Edge<weight_type>>& edges) {
  if (edges.empty()) {
    return {};
  }

  const size_t num_edges = edges.size();
  vertex_data* data = gbbs::new_array_no_init<vertex_data>(num_vertices);
  parallel_for(0, edges[0].from + 1,
               [&](const size_t j) { data[j].offset = 0; });
  parallel_for(1, num_edges, [&](const size_t i) {
    const size_t from = edges[i].from;
    const size_t previous_from = edges[i - 1].from;
    if (from != previous_from) {
      parallel_for(previous_from + 1, from + 1,
                   [&](const size_t j) { data[j].offset = i; });
    }
  });
  parallel_for(edges[num_edges - 1].from + 1, num_vertices,
               [&](const size_t j) { data[j].offset = num_edges; });

  parallel_for(0, num_vertices, [&](const size_t j) {
    const size_t next_offset =
        j == num_vertices - 1 ? num_edges : data[j + 1].offset;
    data[j].degree = next_offset - data[j].offset;
  });

  return data;
}

template <class weight_type>
weight_type string_to_weight(const char* str) {
  if
    constexpr(std::is_integral<weight_type>::value) { return atol(str); }
  else if
    constexpr(std::is_floating_point<weight_type>::value) { return atof(str); }
  else {
    static_assert(always_false<weight_type>,
                  "Converting a char* to weight_type is not supported");
  }
}

/* Returns a tuple containing (n, m, offsets, edges) --- the number of
 * vertices, edges, the vertex offsets, and the edge values, after
 * parsing the input (weighted) graph file. */
template <class weight_type>
std::tuple<size_t, size_t, uintT*, std::tuple<uintE, weight_type>*>
parse_weighted_graph(const char* fname, bool mmap, bool binary, char* bytes,
                     size_t bytes_size) {
  using id_and_weight = std::tuple<uintE, weight_type>;

  uintT* offsets;
  id_and_weight* edges;
  uint64_t n, m;

  sequence<char> S;

  if (!binary) {
    if (bytes == nullptr) {
      if (mmap) {
        std::pair<char*, size_t> MM = mmapStringFromFile(fname);
        S = sequence<char>(MM.second);
        // Cannot mutate the graph unless we copy.
        parallel_for(0, S.size(), [&](size_t i) { S[i] = MM.first[i]; });
        if (munmap(MM.first, MM.second) == -1) {
          perror("munmap");
          exit(-1);
        }
      } else {
        S = readStringFromFile(fname);
      }
    }
    auto tokens = parlay::map_tokens(
        parlay::make_slice(S), [](auto x) { return parlay::make_slice(x); });
    gbbs_debug(std::string header = std::string(tokens[0].begin(), tokens[0].size());
          assert(header == internal::kWeightedAdjGraphHeader););

    uint64_t len = tokens.size() - 1;

    n = parlay::internal::chars_to_int_t<unsigned long>(tokens[1]);
    m = parlay::internal::chars_to_int_t<unsigned long>(tokens[2]);

    if (len != (n + 2 * m + 2)) {
      std::cout << "len = " << len << "\n";
      std::cout << "n = " << n << " m = " << m << "\n";
      std::cout << "should be : " << (n + 2 * m + 2) << "\n";
      assert(false);  // invalid format
    }

    offsets = gbbs::new_array_no_init<uintT>(n + 1);
    edges = gbbs::new_array_no_init<id_and_weight>(2 * m);

    parallel_for(0, n, [&](size_t i) {
      offsets[i] =
          parlay::internal::chars_to_int_t<unsigned long>(tokens[i + 3]);
    });
    offsets[n] = m; /* make sure to set the last offset */
    parallel_for(0, m, [&](size_t i) {
      auto& wgh = tokens[i + n + m + 3];
      auto wgh_str = std::string(wgh.begin(), wgh.size());
      edges[i] = std::make_tuple(
          parlay::internal::chars_to_int_t<uintE>(tokens[i + n + 3]),
          std::stof(wgh_str));
    });
    S.clear();
    tokens.clear();
  } else {
    std::pair<char*, size_t> MM = mmapStringFromFile(fname);
    auto mmap_file = MM.first;

    long* sizes = (long*)mmap_file;
    n = sizes[0], m = sizes[1];

    offsets = (uintT*)(mmap_file + 3 * sizeof(long));
    uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
    edges = (id_and_weight*)(mmap_file + skip);
  }
  return std::make_tuple(n, m, offsets, edges);
}

template <class weight_type>
sequence<Edge<weight_type>> sort_and_dedupe(sequence<Edge<weight_type>> edges) {
  constexpr auto compare_endpoints{
      [](const Edge<weight_type>& left, const Edge<weight_type>& right) {
        return std::tie(left.from, left.to) < std::tie(right.from, right.to);
      }};
  constexpr auto unequal_endpoints{
      [](const Edge<weight_type>& left, const Edge<weight_type>& right) {
        return std::tie(left.from, left.to) != std::tie(right.from, right.to);
      }};
  parlay::sample_sort_inplace(make_slice(edges), compare_endpoints);
  return parlay::pack(
      edges, parlay::delayed_seq<bool>(edges.size(), [&](const size_t i) {
        const auto edge{edges[i]};
        const bool is_self_loop{edge.from == edge.to};
        return !is_self_loop &&
               (i == 0 || unequal_endpoints(edges[i - 1], edge));
      }));
}

}  // namespace internal

}  // namespace gbbs_io
}  // namespace gbbs
