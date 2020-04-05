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

#include "ligra/graph.h"
#include "ligra/io.h"
#include "pbbslib/seq.h"

namespace gbbs_io {

template <class weight_type>
struct edge {
  uintE from;
  uintE to;
  weight_type weight;
};

/* Returns a tuple containing (n, m, offsets, edges) --- the number of
 * vertices, edges, the vertex offsets, and the edge values, after
 * parsing the input (weighted) graph file */
std::tuple<size_t, size_t, uintT*, std::tuple<uintE, intE>*>
parse_weighted_graph(
  char* fname,
  bool mmap,
  char* bytes = nullptr,
  size_t bytes_size = std::numeric_limits<size_t>::max());

symmetric_graph<symmetric_vertex, intE> read_weighted_symmetric_graph(
    char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

asymmetric_graph<asymmetric_vertex, intE> read_weighted_asymmetric_graph(
    char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

/* Returns a tuple containing (n, m, offsets, edges) --- the number of
 * vertices, edges, the vertex offsets, and the edge values, after
 * parsing the input graph file */
std::tuple<size_t, size_t, uintT*, uintE*> parse_unweighted_graph(
    char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

symmetric_graph<symmetric_vertex, pbbslib::empty> read_unweighted_symmetric_graph(
    char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

asymmetric_graph<asymmetric_vertex, pbbslib::empty> read_unweighted_asymmetric_graph(
    char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max());

std::tuple<char*, size_t> parse_compressed_graph(
    char* fname, bool mmap, bool mmapcopy);

template <class weight_type>
symmetric_graph<csv_bytepd_amortized, weight_type>
read_compressed_symmetric_graph(char* fname, bool mmap, bool mmapcopy) {
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

  auto deletion_fn = get_deletion_fn(v_data, bytes);
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
read_compressed_asymmetric_graph(char* fname, bool mmap, bool mmapcopy) {
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

  auto deletion_fn = get_deletion_fn(v_data, v_in_data, bytes);
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
//   <edge 1 first endpoint> <edge 1 second endpoint> <edge 1 weight>
//   <edge 2 first endpoint> <edge 2 second endpoint> <edge 2 weight>
//   <edge 3 first endpoint> <edge 3 second endpoint> <edge 3 weight>
//   ...
//   <edge m first endpoint> <edge m second endpoint> <edge m weight>
pbbs::sequence<edge<intT>> read_weighted_edge_list(const char* filename);

// Read edges from a file that has the following format:
//   <edge 1 first endpoint> <edge 1 second endpoint>
//   <edge 2 first endpoint> <edge 2 second endpoint>
//   <edge 3 first endpoint> <edge 3 second endpoint>
//   ...
//   <edge m first endpoint> <edge m second endpoint>
pbbs::sequence<edge<pbbslib::empty>>
read_unweighted_edge_list(const char* filename);

// Converts edge list into an asymmetric graph.
//
// Adjacency lists of the output graph are sorted by increasing neighbor vertex
// ID. Duplicate edges are removed.
template <class weight_type>
asymmetric_graph<symmetric_vertex, weight_type> edge_list_to_asymmetric_graph(
    const pbbs::sequence<edge<weight_type>>& edge_list);

// Converts a list of undirected edges into a symmetric graph.
//
// Adjacency lists of the output graph are sorted by increasing neighbor vertex
// ID. Duplicate edges are removed.
template <class weight_type>
symmetric_graph<symmetric_vertex, weight_type> edge_list_to_symmetric_graph(
    const pbbs::sequence<edge<weight_type>>& edge_list);

// Write graph in adjacency graph format to file.
template <class Graph>
void write_graph_to_file(const char* filename, Graph* graph);

} // namespace gbbs_io
