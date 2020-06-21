// This code is part of the project
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

#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/io.h"
#include "gbbs/graph_io.h"
#include "gbbs/macros.h"
#include "gbbs/vertex.h"
#include "pbbslib/sample_sort.h"
#include "pbbslib/stlalgs.h"

//#include <libpmem.h>
#include <fcntl.h>
#if defined(__APPLE__)
#else
#include <malloc.h>
#endif
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>


namespace gbbs {
namespace sage_io {

std::pair<char*, size_t> mmap_pmem(const char* filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG(sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  void* ret_arr;
  //size_t mapped_len;
  //int is_pmem;
  /* create a pmem file and memory map it */
//  if ((ret_arr = pmem_map_file(filename, sb.st_size, PMEM_FILE_CREATE, 0666,
//                               &mapped_len, &is_pmem)) == NULL) {
  if ((ret_arr = (mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0))) == NULL) {
    std::cout << "Error on pmem_map_file:\n";
    perror("pmem_map_file");
    exit(1);
  }

  char* p = static_cast<char*>(ret_arr);

  //return std::make_pair(p, mapped_len);
  return std::make_pair(p, sb.st_size);
}

// Currently specialized to 2 sockets, but note that the code below should be
// easy to generalize to more.
template <class weight_type>
symmetric_graph<csv_bytepd_amortized, weight_type>
read_compressed_symmetric_graph(const char* f1, const char* f2) {
  auto [s0, s0_size] = mmap_pmem(f1);
  auto [s1, s1_size] = mmap_pmem(f2);
  if (s0_size != s1_size) {
    std::cout << f1 << " and " << f2 << " have different file lengths, aborting" << std::endl;
    exit(-1);
  }

  long* sizes = (long*)s0;
  uint64_t n = sizes[0], m = sizes[1];

  debug(uint64_t totalSpace = sizes[2];
  std::cout << "# n = " << n << " m = " << m << " totalSpace = " << totalSpace
            << "\n");

  uintT* offsets = (uintT*)(s0 + 3 * sizeof(long));
  uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees = (uintE*)(s0 + skip);
  skip += n * sizeof(intE);

  uchar* edges0 = (uchar*)(s0 + skip);
  uchar* edges1 = (uchar*)(s1 + skip);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = Degrees[i];
  });

  std::function<void()> deletion_fn = [v_data, s0, s1, s0_size] () {
    pbbslib::free_array(v_data);
    gbbs_io::unmmap(s0, s0_size);
    gbbs_io::unmmap(s1, s0_size);
  };
  symmetric_graph<csv_bytepd_amortized, weight_type> G(v_data, n, m, deletion_fn, edges0, edges1);
  return G;
}

template <class weight_type>
asymmetric_graph<cav_bytepd_amortized, weight_type>
read_compressed_asymmetric_graph(const char* f1, const char* f2) {
  auto [s0, s0_size] = mmap_pmem(f1);
  auto [s1, s1_size] = mmap_pmem(f2);
  if (s0_size != s1_size) {
    std::cout << f1 << " and " << f2 << " have different file lengths, aborting" << std::endl;
    exit(-1);
  }

  long* sizes = (long*)s0;
  uint64_t n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  debug(std::cout << "# n = " << n << " m = " << m << " totalSpace = " << totalSpace
            << "\n");

  uintT* offsets = (uintT*)(s0 + 3 * sizeof(long));
  uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees = (uintE*)(s0 + skip);
  skip += n * sizeof(intE);

  uchar* out_edges0 = (uchar*)(s0 + skip);
  uchar* out_edges1 = (uchar*)(s1 + skip);

  skip += totalSpace;
  uchar* inData = (uchar*)(s0 + skip);
  sizes = (long*)inData;
  skip += sizeof(long);
  uintT* inOffsets = (uintT*)(s0 + skip);
  skip += (n + 1) * sizeof(uintT);
  uintE* inDegrees = (uintE*)(s0 + skip);
  skip += n * sizeof(uintE);
  uchar* in_edges0 = (uchar*)(s0 + skip);
  uchar* in_edges1 = (uchar*)(s1 + skip);

  auto v_data = pbbs::new_array_no_init<vertex_data>(n);
  auto v_in_data = pbbs::new_array_no_init<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    v_data[i].offset = offsets[i];
    v_data[i].degree = Degrees[i];

    v_in_data[i].offset = inOffsets[i];
    v_in_data[i].degree = inDegrees[i];
  });

  std::function<void()> deletion_fn = [v_data, v_in_data, s0, s1, s0_size] () {
    pbbslib::free_array(v_data);
    pbbslib::free_array(v_in_data);
    gbbs_io::unmmap(s0, s0_size);
    gbbs_io::unmmap(s1, s0_size);
  };
  asymmetric_graph<cav_bytepd_amortized, weight_type> G(v_data, v_in_data, n, m, deletion_fn, out_edges0, in_edges0, out_edges1, in_edges1);
  return G;
}

template <class weight_type>
symmetric_graph<symmetric_vertex, weight_type>
read_symmetric_binary_graph(const char* f1, const char* f2) {
  auto [s0, s0_size] = mmap_pmem(f1);
  auto [s1, s1_size] = mmap_pmem(f2);
  if (s0_size != s1_size) {
    std::cout << f1 << " and " << f2 << " have different file lengths, aborting" << std::endl;
    exit(-1);
  }

  long* sizes = (long*)s0;
  uint64_t n = sizes[0], m = sizes[1];

  std::cout << "n = " << n << " m = " << m << std::endl;
  using edge_type = std::tuple<uintE, weight_type>;

  uintT* offsets = (uintT*)(s0 + 3 * sizeof(long));
  uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);

  edge_type* edges0 = (edge_type*)(s0 + skip);
  edge_type* edges1 = (edge_type*)(s1 + skip);

  vertex_data* V = pbbslib::new_array_no_init<vertex_data>(n);

  par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
    uint64_t o = offsets[i];
    uint64_t next_o = offsets[i+1];
    uintT d = next_o - o;
    V[i].degree = d;
    V[i].offset = o;
  });

  auto deletion_fn = [V, s0, s1, s0_size]() {
    pbbslib::free_array(V);
    gbbs_io::unmmap(s0, s0_size);
    gbbs_io::unmmap(s1, s0_size);
  };
  return symmetric_graph<symmetric_vertex, weight_type>(V, n, m, deletion_fn, edges0, edges1);
}

template <class weight_type>
asymmetric_graph<asymmetric_vertex, weight_type>
read_asymmetric_binary_graph(const char* f1, const char* f2) {
  auto [s0, s0_size] = mmap_pmem(f1);
  auto [s1, s1_size] = mmap_pmem(f2);
  if (s0_size != s1_size) {
    std::cout << f1 << " and " << f2 << " have different file lengths, aborting" << std::endl;
    exit(-1);
  }

  long* sizes = (long*)s0;
  uint64_t n = sizes[0], m = sizes[1];

  std::cout << "n = " << n << " m = " << m << std::endl;
  using edge_type = std::tuple<uintE, weight_type>;

  uintT* offsets = (uintT*)(s0 + 3 * sizeof(long));
  uint64_t skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);

  edge_type* edges0 = (edge_type*)(s0 + skip);
  edge_type* edges1 = (edge_type*)(s1 + skip);

  skip = 3 * sizeof(long) + (n + 1) * sizeof(intT) + m * sizeof(edge_type);

  uintT* inOffsets = (uintT*)(s0 + skip);
  skip += (n+1) * sizeof(intT);
  edge_type* in_edges0 = (edge_type*)(s0 + skip);
  edge_type* in_edges1 = (edge_type*)(s1 + skip);

  auto v_data = pbbslib::new_array_no_init<vertex_data>(n);
  auto v_in_data = pbbs::new_array_no_init<vertex_data>(n);

  par_for(0, n, pbbslib::kSequentialForThreshold, [&](size_t i) {
    uint64_t o = offsets[i];
    uint64_t next_o = offsets[i+1];
    uintT d = next_o - o;
    v_data[i].degree = d;
    v_data[i].offset = o;

    o = inOffsets[i];
    next_o = inOffsets[i+1];
    uintT in_d = next_o - o;
    v_in_data[i].degree = in_d;
    v_in_data[i].offset = o;
  });

  auto deletion_fn = [v_data, v_in_data, s0, s1, s0_size]() {
    pbbslib::free_array(v_data);
    pbbslib::free_array(v_in_data);
    gbbs_io::unmmap(s0, s0_size);
    gbbs_io::unmmap(s1, s0_size);
  };
  return asymmetric_graph<asymmetric_vertex, weight_type>(v_data, v_in_data, n, m, deletion_fn, edges0, in_edges0, edges1, in_edges1);
}

}
}  // namespace gbbs
