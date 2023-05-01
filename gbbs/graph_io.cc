#include "graph_io.h"

#include <sys/mman.h>

namespace gbbs {
namespace gbbs_io {

typedef std::pair<uintE, uintE> intPair;

namespace internal {

// Starting from the current position, skips all consecutive lines in the stream
// that start with '#' or are empty.
//
// The intent here is that lines starting with '#' are interpreted to be
// comments that should be ignored.
void skip_ifstream_comments(std::ifstream* stream) {
  std::string line;
  while (*stream) {
    std::streampos current_position = stream->tellg();
    std::getline(*stream, line);
    if (!(line.empty() || line[0] == '#')) {
      stream->seekg(current_position);
      return;
    }
  }
}

}  // namespace internal

template <>
Edge<gbbs::empty>::Edge(uintE _from, uintE _to)
    : from(_from), to(_to) {}

std::tuple<size_t, size_t, uintT*, uintE*> parse_unweighted_graph(
    const char* fname, bool mmap, bool binary, char* bytes, size_t bytes_size) {
  uintT* offsets;
  uintE* edges;
  uint64_t n, m;

  if (!binary) {
    sequence<char> S;

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
    sequence<slice<char>> tokens = parlay::map_tokens(
        parlay::make_slice(S), [](auto x) { return parlay::make_slice(x); });

    gbbs_debug(std::string header = std::string(tokens[0].begin(), tokens[0].size());
          assert(header == internal::kUnweightedAdjGraphHeader););

    n = parlay::internal::chars_to_int_t<unsigned long>(tokens[1]);
    m = parlay::internal::chars_to_int_t<unsigned long>(tokens[2]);

    gbbs_debug(std::cout << "# n = " << n << " m = " << m
                    << " len = " << (tokens.size() - 1) << "\n";
          uint64_t len = tokens.size() - 1; assert(len == n + m + 2););

    offsets = gbbs::new_array_no_init<uintT>(n + 1);
    edges = gbbs::new_array_no_init<uintE>(m);

    parallel_for(0, n, [&](size_t i) {
      offsets[i] =
          parlay::internal::chars_to_int_t<unsigned long>(tokens[i + 3]);
    });
    offsets[n] = m; /* make sure to set the last offset */
    parallel_for(0, m, [&](size_t i) {
      edges[i] =
          parlay::internal::chars_to_int_t<unsigned long>(tokens[i + n + 3]);
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
    edges = (uintE*)(mmap_file + skip);
  }

  return std::make_tuple(n, m, offsets, edges);
}

symmetric_graph<symmetric_vertex, gbbs::empty> read_unweighted_symmetric_graph(
    const char* fname, bool mmap, bool binary, char* bytes, size_t bytes_size) {
  size_t n, m;
  uintT* offsets;
  uintE* edges;
  std::tie(n, m, offsets, edges) =
      parse_unweighted_graph(fname, mmap, binary, bytes, bytes_size);

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i + 1] - v_data[i].offset;
  });
  if (!binary) {
    gbbs::free_array(offsets, n + 1);
  }

  return symmetric_graph<symmetric_vertex, gbbs::empty>(
      v_data, n, m,
      [=]() {
        gbbs::free_array(v_data, n);
        if (!binary) gbbs::free_array(edges, m);
      },
      (std::tuple<uintE, gbbs::empty>*)edges);
}

asymmetric_graph<asymmetric_vertex, gbbs::empty>
read_unweighted_asymmetric_graph(const char* fname, bool mmap, char* bytes,
                                 size_t bytes_size) {
  size_t n, m;
  uintT* offsets;
  uintE* edges;
  // Reads in just the out-edges and computes the in-edges.
  std::tie(n, m, offsets, edges) =
      parse_unweighted_graph(fname, mmap, false, bytes, bytes_size);

  auto v_data = gbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&](size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i + 1] - v_data[i].offset;
  });
  gbbs::free_array(offsets, n + 1);

  /* construct transpose of the graph */
  sequence<uintT> tOffsets = sequence<uintT>::uninitialized(n + 1);
  parallel_for(0, n + 1, [&](size_t i) { tOffsets[i] = INT_T_MAX; });
  intPair* temp = gbbs::new_array_no_init<intPair>(m);
  parallel_for(0, n, [&](size_t i) {
    uintT o = v_data[i].offset;
    uintT deg = v_data[i].degree;
    for (uintT j = 0; j < deg; j++) {
      temp[o + j] = std::make_pair((edges + o)[j], i);
    }
  });

  auto temp_seq = gbbs::make_slice(temp, m);
  parlay::integer_sort_inplace(temp_seq,
                               [&](const intPair& p) { return p.first; });

  tOffsets[temp[0].first] = 0;
  uintE* inEdges = gbbs::new_array_no_init<uintE>(m);
  inEdges[0] = temp[0].second;
  parallel_for(1, m, [&](size_t i) {
    inEdges[i] = temp[i].second;
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

  return asymmetric_graph<asymmetric_vertex, gbbs::empty>(
      v_data, v_in_data, n, m,
      [=]() {
        gbbs::free_array(v_data, n);
        gbbs::free_array(v_in_data, n);
        gbbs::free_array(inEdges, m);
        gbbs::free_array(edges, m);
      },
      (std::tuple<uintE, gbbs::empty>*)edges,
      (std::tuple<uintE, gbbs::empty>*)inEdges);
}

asymmetric_graph<asymmetric_vertex, gbbs::empty>
read_unweighted_asymmetric_graph(const char* fname, bool mmap, bool binary,
                                 char* bytes, size_t bytes_size) {
  if (!binary) {
    return read_unweighted_asymmetric_graph(fname, mmap, bytes, bytes_size);
  } else {  // binary input
    size_t n, m;
    std::pair<char*, size_t> MM = mmapStringFromFile(fname);
    auto mmap_file = MM.first;
    long* sizes = (long*)mmap_file;
    n = sizes[0], m = sizes[1];
    uintT* out_offsets = (uintT*)(mmap_file + 3 * sizeof(long));
    uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(uintT);
    uintE* out_edges = (uintE*)(mmap_file + skip);
    skip += m * sizeof(uintE);

    uintT* in_offsets = (uintT*)(mmap_file + skip + 3 * sizeof(long));
    skip += 3 * sizeof(long) + (n + 1) * sizeof(uintT);
    uintE* in_edges = (uintE*)(mmap_file + skip);

    auto v_out_data = gbbs::new_array_no_init<vertex_data>(n);
    parallel_for(0, n, [&](size_t i) {
      v_out_data[i].offset = out_offsets[i];
      v_out_data[i].degree = out_offsets[i + 1] - v_out_data[i].offset;
    });

    auto v_in_data = gbbs::new_array_no_init<vertex_data>(n);
    parallel_for(0, n, [&](size_t i) {
      v_in_data[i].offset = in_offsets[i];
      v_in_data[i].degree = in_offsets[i + 1] - v_in_data[i].offset;
    });

    return asymmetric_graph<asymmetric_vertex, gbbs::empty>(
        v_out_data, v_in_data, n, m,
        [=]() {
          gbbs::free_array(v_out_data, n);
          gbbs::free_array(v_in_data, n);
        },
        (std::tuple<uintE, gbbs::empty>*)out_edges,
        (std::tuple<uintE, gbbs::empty>*)in_edges);
  }
}

std::tuple<char*, size_t> parse_compressed_graph(const char* fname, bool mmap) {
  char* bytes;
  size_t bytes_size;

  if (mmap) {
    std::tie(bytes, bytes_size) = mmapStringFromFile(fname);
  } else {
    std::tie(bytes, bytes_size) = read_o_direct(fname);
  }
  return std::make_tuple(bytes, bytes_size);
}

std::vector<Edge<gbbs::empty>> read_unweighted_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  internal::skip_ifstream_comments(&file);

  std::vector<Edge<gbbs::empty>> edge_list;
  uintE from;
  uintE to;
  while (file >> from >> to) {
    edge_list.emplace_back(from, to);
  }
  return edge_list;
}

}  // namespace gbbs_io
}  // namespace gbbs
