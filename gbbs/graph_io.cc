#include "graph_io.h"

#include <sys/mman.h>

#include "pbbslib/strings/string_basics.h"

namespace gbbs {
namespace gbbs_io {

typedef std::pair<uintE, uintE> intPair;
typedef std::pair<uintE, std::pair<uintE, intE>> intTriple;

namespace {

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

}  // namespace

template <>
Edge<pbbslib::empty>::Edge(const uintE _from, const uintE _to)
  : from(_from)
  , to(_to) {}

std::tuple<size_t, size_t, uintT*, std::tuple<uintE, intE>*>
parse_weighted_graph(const char* fname, bool mmap, char* bytes, size_t bytes_size) {
  sequence<char*> tokens;
  sequence<char> S;
  if (bytes == nullptr) {
    if (mmap) {
      std::pair<char*, size_t> MM = mmapStringFromFile(fname);
      S = sequence<char>(MM.second);
      // Cannot mutate the graph unless we copy.
      par_for(0, S.size(), pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { S[i] = MM.first[i]; });
      if (munmap(MM.first, MM.second) == -1) {
        perror("munmap");
        exit(-1);
      }
    } else {
      S = readStringFromFile(fname);
    }
  }
  tokens = pbbs::tokenize(S, [] (const char c) { return pbbs::is_space(c); });
  assert(tokens[0] == internal::kWeightedAdjGraphHeader);

  uint64_t len = tokens.size() - 1;
  uint64_t n = atol(tokens[1]);
  uint64_t m = atol(tokens[2]);
  if (len != (n + 2 * m + 2)) {
    std::cout << tokens[0] << "\n";
    std::cout << "len = " << len << "\n";
    std::cout << "n = " << n << " m = " << m << "\n";
    std::cout << "should be : " << (n + 2 * m + 2) << "\n";
    assert(false);  // invalid format
  }

  uintT* offsets = pbbslib::new_array_no_init<uintT>(n+1);
  using id_and_weight = std::tuple<uintE, intE>;
  std::tuple<uintE, intE>* edges = pbbslib::new_array_no_init<id_and_weight>(2 * m);

  parallel_for(0, n, [&] (size_t i) { offsets[i] = atol(tokens[i + 3]); });
  offsets[n] = m; /* make sure to set the last offset */
  parallel_for(0, m, [&] (size_t i) {
    edges[i] = std::make_tuple(atol(tokens[i + n + 3]),
                               atol(tokens[i + n + m + 3]));
  });
  S.clear();
  tokens.clear();

  return std::make_tuple(n, m, offsets, edges);
}

symmetric_graph<symmetric_vertex, intE> read_weighted_symmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes,
    size_t bytes_size) {
  size_t n, m;
  uintT* offsets;
  std::tuple<uintE, intE>* edges;
  std::tie(n, m, offsets, edges) = parse_weighted_graph(fname, mmap, bytes, bytes_size);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i+1]-v_data[i].offset;
  });
  pbbs::free_array(offsets);

  return symmetric_graph<symmetric_vertex, intE>(v_data, n, m, [=]() { pbbslib::free_arrays(v_data, edges); }, edges);
}

asymmetric_graph<asymmetric_vertex, intE> read_weighted_asymmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes,
    size_t bytes_size) {
  using id_and_weight = std::tuple<uintE, intE>;

  size_t n, m;
  uintT* offsets;
  std::tuple<uintE, intE>* edges;
  std::tie(n, m, offsets, edges) = parse_weighted_graph(fname, mmap, bytes, bytes_size);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i+1]-v_data[i].offset;
  });
  pbbs::free_array(offsets);

  uintT* tOffsets = pbbslib::new_array_no_init<uintT>(n+1);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { tOffsets[i] = INT_T_MAX; });
  intTriple* temp = pbbslib::new_array_no_init<intTriple>(m);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintT o = v_data[i].offset;
    uintE deg = v_data[i].degree;
    for (uintT j = 0; j < deg; j++) {
      auto& cur_edge = (edges + o)[j];
      temp[o + j] = std::make_pair(std::get<0>(cur_edge),
                                   std::make_pair((uintE)i, std::get<1>(cur_edge)));
    }
  });

  auto temp_seq = pbbslib::make_sequence(temp, m);
  pbbslib::integer_sort_inplace(temp_seq.slice(), [&] (const intTriple& p) { return p.first; }, pbbs::log2_up(n));

  tOffsets[temp[0].first] = 0;
  id_and_weight* inEdges = pbbslib::new_array_no_init<id_and_weight>(m);
  inEdges[0] = std::make_tuple(temp[0].second.first, temp[0].second.second);

  par_for(1, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    inEdges[i] = std::make_tuple(temp[i].second.first, temp[i].second.second);
    if (temp[i].first != temp[i - 1].first) {
      tOffsets[temp[i].first] = i;
    }
  });

  pbbslib::free_array(temp);

  // fill in offsets of degree 0 vertices by taking closest non-zero
  // offset to the right
  debug(cout << "# scan I back " << std::endl;);
  auto t_seq = pbbslib::make_sequence(tOffsets, n+1).rslice();
  auto M = pbbslib::minm<uintT>();
  M.identity = m;
  pbbslib::scan_inplace(t_seq, M, pbbslib::fl_scan_inclusive);

  auto v_in_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_in_data[i].offset = tOffsets[i];
    v_in_data[i].degree = tOffsets[i+1]-v_in_data[i].offset;
  });
  pbbs::free_array(tOffsets);

  return asymmetric_graph<asymmetric_vertex, intE>(v_data, v_in_data, n, m, [=]() { pbbslib::free_arrays(v_data, v_in_data, edges, inEdges); }, edges, inEdges);
}

std::tuple<size_t, size_t, uintT*, uintE*> parse_unweighted_graph(
    const char* fname,
    bool mmap,
    char* bytes,
    size_t bytes_size) {
  sequence<char*> tokens;
  sequence<char> S;

  if (bytes == nullptr) {
    if (mmap) {
      std::pair<char*, size_t> MM = mmapStringFromFile(fname);
      S = sequence<char>(MM.second);
      // Cannot mutate the graph unless we copy.
      par_for(0, S.size(), pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { S[i] = MM.first[i]; });
      if (munmap(MM.first, MM.second) == -1) {
        perror("munmap");
        exit(-1);
      }
    } else {
      S = readStringFromFile(fname);
    }
  }
  tokens = pbbs::tokenize(S, [] (const char c) { return pbbs::is_space(c); });

  assert(tokens[0] == internal::kUnweightedAdjGraphHeader);

  uint64_t n = atol(tokens[1]);
  uint64_t m = atol(tokens[2]);

  debug(std::cout << "# n = " << n << " m = " << m << " len = " << (tokens.size() - 1) << "\n";
  uint64_t len = tokens.size() - 1;
  assert(len == n + m + 2););

  uintT* offsets = pbbslib::new_array_no_init<uintT>(n+1);
  uintE* edges = pbbslib::new_array_no_init<uintE>(m);

  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { offsets[i] = atol(tokens[i + 3]); });
  offsets[n] = m; /* make sure to set the last offset */
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { edges[i] = atol(tokens[i + n + 3]); });
  S.clear();
  tokens.clear();

  return std::make_tuple(n, m, offsets, edges);
}

symmetric_graph<symmetric_vertex, pbbslib::empty> read_unweighted_symmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes,
    size_t bytes_size) {
  size_t n, m;
  uintT* offsets;
  uintE* edges;
  std::tie(n, m, offsets, edges) = parse_unweighted_graph(fname, mmap, bytes, bytes_size);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i+1]-v_data[i].offset;
  });
  pbbs::free_array(offsets);

  return symmetric_graph<symmetric_vertex, pbbs::empty>(
      v_data, n, m, [=](){ pbbslib::free_arrays(v_data, edges); }, (std::tuple<uintE, pbbs::empty>*)edges);
}

asymmetric_graph<asymmetric_vertex, pbbslib::empty> read_unweighted_asymmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes,
    size_t bytes_size) {
  size_t n, m;
  uintT* offsets;
  uintE* edges;
  std::tie(n, m, offsets, edges) = parse_unweighted_graph(fname, mmap, bytes, bytes_size);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = offsets[i+1]-v_data[i].offset;
  });
  pbbs::free_array(offsets);

  /* construct transpose of the graph */
  uintT* tOffsets = pbbslib::new_array_no_init<uintT>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { tOffsets[i] = INT_T_MAX; });
  intPair* temp = pbbslib::new_array_no_init<intPair>(m);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintT o = v_data[i].offset;
    uintT deg = v_data[i].degree;
    for (uintT j = 0; j < deg; j++) {
      temp[o + j ] = std::make_pair((edges + o)[j], i);
    }
  });

  auto temp_seq = pbbslib::make_sequence(temp, m);
  pbbslib::integer_sort_inplace(temp_seq.slice(), [&] (const intPair& p) { return p.first; }, pbbs::log2_up(n));

  tOffsets[temp[0].first] = 0;
  uintE* inEdges = pbbslib::new_array_no_init<uintE>(m);
  inEdges[0] = temp[0].second;
  par_for(1, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    inEdges[i] = temp[i].second;
    if (temp[i].first != temp[i - 1].first) {
      tOffsets[temp[i].first] = i;
    }
  });

  pbbslib::free_array(temp);

  // fill in offsets of degree 0 vertices by taking closest non-zero
  // offset to the right
  auto t_seq = pbbslib::make_sequence(tOffsets, n).rslice();
  auto M = pbbslib::minm<uintT>();
  M.identity = m;
  pbbslib::scan_inplace(t_seq, M, pbbslib::fl_scan_inclusive);

  auto v_in_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_in_data[i].offset = tOffsets[i];
    v_in_data[i].degree = tOffsets[i+1]-v_in_data[i].offset;
  });
  pbbs::free_array(tOffsets);

  return asymmetric_graph<asymmetric_vertex, pbbs::empty>(
      v_data, v_in_data, n, m, [=]() {pbbslib::free_arrays(v_data, v_in_data, inEdges, edges);},(std::tuple<uintE, pbbs::empty>*)edges, (std::tuple<uintE, pbbs::empty>*)inEdges);
}

std::tuple<char*, size_t> parse_compressed_graph(
    const char* fname, bool mmap, bool mmapcopy) {
  char* bytes;
  size_t bytes_size;

  if (mmap) {
    std::tie(bytes, bytes_size) = mmapStringFromFile(fname);
    if (mmapcopy) {
      debug(std::cout << "# Copying compressed graph due to mmapcopy being set."
                << "\n";);
      char* next_bytes = pbbslib::new_array_no_init<char>(bytes_size);
      par_for(0, bytes_size, pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { next_bytes[i] = bytes[i]; });
      if (munmap(bytes, bytes_size) == -1) {
        perror("munmap");
        exit(-1);
      }
      bytes = next_bytes;
    }
  } else {
    std::tie(bytes, bytes_size) = read_o_direct(fname);
  }
  return std::make_tuple(bytes, bytes_size);
}

std::vector<Edge<intT>> read_weighted_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  skip_ifstream_comments(&file);

  std::vector<Edge<intT>> edge_list;
  uintE from;
  uintE to;
  intT weight;
  while (file >> from >> to >> weight) {
    edge_list.emplace_back(from, to, weight);
  }
  return edge_list;
}

std::vector<Edge<pbbslib::empty>>
read_unweighted_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  skip_ifstream_comments(&file);

  std::vector<Edge<pbbslib::empty>> edge_list;
  uintE from;
  uintE to;
  while (file >> from >> to) {
    edge_list.emplace_back(from, to);
  }
  return edge_list;
}

}  // namespace gbbs_io
}  // namespace gbbs
