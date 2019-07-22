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

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "bridge.h"
#include "pbbs_strings.h"
#include "graph.h"

#ifdef NVM
#include <libpmem.h>
#endif

typedef std::pair<uintE, uintE> intPair;
typedef std::pair<uintE, std::pair<uintE, intE>> intTriple;

char* compressed_mmap_bytes = nullptr;
size_t compressed_mmap_bytes_size = 0L;

template <class E>
struct pairFirstCmp {
  bool operator()(std::pair<uintE, E> a, std::pair<uintE, E> b) {
    return a.first < b.first;
  }
};

template <class E>
struct getFirst {
  uintE operator()(std::pair<uintE, E> a) { return a.first; }
};

template <class W,
          typename std::enable_if<!std::is_same<W, intE>::value, int>::type = 0>
inline std::string print_wgh(W wgh) {
  return "";
}

template <class W,
          typename std::enable_if<std::is_same<W, intE>::value, int>::type = 0>
inline std::string print_wgh(W wgh) {
  return std::to_string(wgh);
}


// returns a pointer and a length
// Note that the caller is responsible for unmmap'ing the file
inline std::pair<char*, size_t> mmapStringFromFile(const char* filename) {
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
  char* p =
      static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
  return std::make_pair(p, n);
}

#ifdef NVM
// used by the nvm graph algorithms (see pmdk).
inline std::pair<char*, size_t> pmem_from_file(const char* filename) {
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
	size_t mapped_len;
	int is_pmem;
  /* create a pmem file and memory map it */
  if ((ret_arr = pmem_map_file(filename, sb.st_size, PMEM_FILE_CREATE, 0666, &mapped_len, &is_pmem)) == NULL) {
		std::cout << "Error on pmem_map_file:\n";
		perror("pmem_map_file");
		exit(1);
  }

	cout << "is_pmem = " << is_pmem << endl;
	cout << "mapped_len = " << mapped_len << endl;

	char* p = static_cast<char*>(ret_arr);

  return std::make_pair(p, mapped_len);
}
#endif

void unmmap_if_needed() {
  if (compressed_mmap_bytes) {
    if (munmap(compressed_mmap_bytes, compressed_mmap_bytes_size) == -1) {
      perror("munmap");
      exit(-1);
    }
  }
}

inline sequence<char> readStringFromFile(char* fileName) {
  std::ifstream file(fileName, std::ios::in | std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    debug(std::cout << "Unable to open file: " << fileName << "\n";);
    abort();
  }
  uint64_t end = file.tellg();
  file.seekg(0, std::ios::beg);
  uint64_t n = end - file.tellg();
  auto bytes = sequence<char>(n); // n+1?
  file.read(bytes.begin(), n);
  file.close();
  return bytes;
}

template <template <typename W> class vertex>
inline graph<vertex, intE> readWeightedGraph(
    char* fname, bool isSymmetric, bool mmap, char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max()) {
  using W = intE;
  using wvtx = vertex<W>;
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
      sequence<char> S = readStringFromFile(fname);
    }
  }
  tokens = pbbslib::tokenize(S, [] (const char c) { return pbbs::is_space(c); });
  assert(tokens[0] == (std::string) "WeightedAdjacencyGraph");

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

  uintT* offsets = pbbslib::new_array_no_init<uintT>(n);
  using VW = std::tuple<uintE, intE>;
  std::tuple<uintE, intE>* edges = pbbslib::new_array_no_init<VW>(2 * m);

  {
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { offsets[i] = atol(tokens[i + 3]); });
  }
  {
    par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      edges[i] = std::make_tuple(atol(tokens[i + n + 3]),
                                 atol(tokens[i + n + m + 3]));
    });
  }
  S.clear();
  tokens.clear();
  // W.clear(); // to deal with performance bug in malloc

  wvtx* v = pbbslib::new_array_no_init<wvtx>(n);

  {
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      uintT o = offsets[i];
      uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
      v[i].setOutDegree(l);
      v[i].setOutNeighbors(edges + o);
    });
  }

  if (!isSymmetric) {
    uintT* tOffsets = pbbslib::new_array_no_init<uintT>(n);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { tOffsets[i] = INT_T_MAX; });
    intTriple* temp = pbbslib::new_array_no_init<intTriple>(m);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      uintT o = offsets[i];
      for (uintT j = 0; j < v[i].getOutDegree(); j++) {
        temp[o + j] = std::make_pair(v[i].getOutNeighbor(j),
                                     std::make_pair(i, v[i].getOutWeight(j)));
      }
    });
    pbbslib::free_array(offsets);

    auto temp_seq = pbbslib::make_sequence(temp, m);
    pbbslib::integer_sort_inplace(temp_seq.slice(), [&] (const intTriple& p) { return p.first; }, pbbs::log2_up(n));

    tOffsets[temp[0].first] = 0;
    VW* inEdges = pbbslib::new_array_no_init<VW>(m);
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

    debug(cout << "scan I back " << endl;);
    auto t_seq = pbbslib::make_sequence(tOffsets, n).rslice();
    auto M = pbbslib::minm<uintT>();
    M.identity = m;
    pbbslib::scan_inplace(t_seq, M, pbbslib::fl_scan_inclusive);

    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      uintT o = tOffsets[i];
      uintT l = ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
      v[i].setInDegree(l);
      v[i].setInNeighbors(inEdges + o);
    });

    pbbslib::free_array(tOffsets);
    return graph<vertex, W>(v, n, m, get_deletion_fn(v, inEdges, edges));
  } else {
    pbbslib::free_array(offsets);
    return graph<vertex, W>(v, n, m, get_deletion_fn(v, edges));
  }
}

template <template <typename W> class vertex>
inline graph<vertex, pbbslib::empty> readUnweightedGraph(
    char* fname, bool isSymmetric, bool mmap, char* bytes = nullptr,
    size_t bytes_size = std::numeric_limits<size_t>::max()) {
  using W = pbbslib::empty;
  using wvtx = vertex<W>;
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
  tokens = pbbslib::tokenize(S, [] (const char c) { return pbbs::is_space(c); });

  assert(tokens[0] == (std::string) "AdjacencyGraph");
  // TODO(laxmand): ensure that S is properly freed here

  // uint64_t len = tokens.size() - 1;
  uint64_t n = atol(tokens[1]);
  uint64_t m = atol(tokens[2]);

  debug(std::cout << "n = " << n << " m = " << m << " len = " << (tokens.size() - 1) << "\n";);
  // assert(len == n + m + 2);

  uintT* offsets = pbbslib::new_array_no_init<uintT>(n);
  uintE* edges = pbbslib::new_array_no_init<uintE>(m);

  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { offsets[i] = atol(tokens[i + 3]); });
  par_for(0, m, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { edges[i] = atol(tokens[i + n + 3]); });
  S.clear();
  tokens.clear();
  // W.clear();  // to deal with performance bug in malloc

  wvtx* v = pbbslib::new_array_no_init<wvtx>(n);

  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uintT o = offsets[i];
    uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
    v[i].setOutDegree(l);
    v[i].setOutNeighbors(((std::tuple<uintE, pbbslib::empty>*)(edges + o)));
  });

  if (!isSymmetric) {
    uintT* tOffsets = pbbslib::new_array_no_init<uintT>(n);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { tOffsets[i] = INT_T_MAX; });
    intPair* temp = pbbslib::new_array_no_init<intPair>(m);
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      uintT o = offsets[i];
      for (uintT j = 0; j < v[i].getOutDegree(); j++) {
        temp[o + j] = std::make_pair(v[i].getOutNeighbor(j), i);
      }
    });
    pbbslib::free_array(offsets);

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

    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
      uintT o = tOffsets[i];
      uintT l = ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
      v[i].setInDegree(l);
      v[i].setInNeighbors((std::tuple<uintE, pbbslib::empty>*)(inEdges + o));
    });

    pbbslib::free_array(tOffsets);

    return graph<vertex, W>(
        v, n, m, get_deletion_fn(v, inEdges, edges));
  } else {
    pbbslib::free_array(offsets);
    return graph<vertex, W>(
        v, n, m, get_deletion_fn(v, edges));
  }
}


// Handles both unweighted and weighted graphs.
template <template <typename W> class vertex, class W>
inline graph<vertex, W> readCompressedGraph(
    char* fname, bool isSymmetric, bool mmap, bool mmapcopy,
    char* bytes = nullptr, size_t bytes_size = std::numeric_limits<size_t>::max()) {
  using w_vertex = vertex<W>;

#ifndef NVM
  char* s;
  if (bytes == nullptr) {
    if (mmap) {
      std::pair<char*, size_t> S = mmapStringFromFile(fname);
      s = S.first;
      if (mmapcopy) {
        debug(std::cout << "Copying compressed graph"
                  << "\n";);
        // Cannot mutate graph unless we copy.
        char* bytes = pbbslib::new_array_no_init<char>(S.second);
        par_for(0, S.second, pbbslib::kSequentialForThreshold, [&] (size_t i)
                        { bytes[i] = S.first[i]; });
        if (munmap(S.first, S.second) == -1) {
          perror("munmap");
          exit(-1);
        }
        s = bytes;
      } else {
        compressed_mmap_bytes = S.first;
        compressed_mmap_bytes_size = S.second;
      }
    } else { // read file using O_DIRECT
      int fd;
      if ((fd = open(fname, O_RDONLY | O_DIRECT)) != -1) {
        debug(std::cout << "input opened!"
                  << "\n";);
      } else {
        std::cout << "can't open input file!";
      }
      //    posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);

      size_t fsize = lseek(fd, 0, SEEK_END);
      lseek(fd, 0, 0);
      s = (char*)memalign(4096 * 2, fsize + 4096);

      debug(std::cout << "fsize = " << fsize << "\n";);

      size_t sz = 0;

      size_t pgsize = getpagesize();
      debug(std::cout << "pgsize = " << pgsize << "\n";);

      size_t read_size = 1024 * 1024 * 1024;
      if (sz + read_size > fsize) {
        size_t k = std::ceil((fsize - sz) / pgsize);
        read_size = std::max(k * pgsize, pgsize);
        debug(std::cout << "set read size to: " << read_size << " " << (fsize - sz)
                  << " bytes left"
                  << "\n";);
      }

      while (sz + read_size < fsize) {
        void* buf = s + sz;
        debug(std::cout << "reading: " << read_size << "\n";);
        sz += read(fd, buf, read_size);
        debug(std::cout << "read: " << sz << " bytes"
                  << "\n";);
        if (sz + read_size > fsize) {
          size_t k = std::ceil((fsize - sz) / pgsize);
          read_size = std::max(k * pgsize, pgsize);
          debug(std::cout << "set read size to: " << read_size << " " << (fsize - sz)
                    << " bytes left"
                    << "\n";);
        }
      }
      if (sz < fsize) {
        debug(std::cout << "last read: rem = " << (fsize - sz) << "\n";);
        void* buf = s + sz;
        sz += read(fd, buf, pgsize);
        debug(std::cout << "read " << sz << " bytes "
                  << "\n";);
      }
      close(fd);
    }
  } else {
    s = bytes;
  }
  char* s0 = s;
#else
  // TODO(laxmand): there has to be a cleaner way to do this.
  std::pair<char*, size_t> S0 = pmem_from_file("/mnt/pmem12/hyperlink2012_sym.bytepda");
  char* s0 = S0.first;

  std::pair<char*, size_t> S1 = pmem_from_file("/mnt/pmem13/hyperlink2012_sym.bytepda");
  char* s1 = S1.first;
#endif

  long* sizes = (long*)s0;
  uint64_t n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  debug(std::cout << "n = " << n << " m = " << m << " totalSpace = " << totalSpace
            << "\n";
  std::cout << "reading file..."
            << "\n";);

  uintT* offsets0 = (uintT*)(s0 + 3 * sizeof(long));
  uint64_t skip0 = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees0 = (uintE*)(s0 + skip0);
  skip0 += n * sizeof(intE);
  uchar* edges0 = (uchar*)(s0 + skip0);
  w_vertex* V0 = pbbslib::new_array_no_init<w_vertex>(n);

#ifdef NVM
  uintT* offsets1 = (uintT*)(s1 + 3 * sizeof(long));
  uint64_t skip1 = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees1 = (uintE*)(s1 + skip1);
  skip1 += n * sizeof(intE);
  uchar* edges1 = (uchar*)(s1 + skip1);
  w_vertex* V1 = pbbslib::new_array_no_init<w_vertex>(n);
#endif

  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    uint64_t o = offsets0[i];
    uintT d = Degrees0[i];
    V0[i].setOutDegree(d);
    V0[i].setOutNeighbors(edges0 + o);

#ifdef NVM
    V1[i].setOutDegree(d);
    V1[i].setOutNeighbors(edges1 + o);
#endif
  });
  auto deletion_fn = get_deletion_fn(V0, s0);
  if (mmap && !mmapcopy) {
#ifndef NVM
    deletion_fn = [V0] () {
      pbbslib::free_array(V0);
      unmmap_if_needed();
    };
#else
    deletion_fn = [V0, V1] () {
      pbbslib::free_array(V0);
      pbbslib::free_array(V1);
      unmmap_if_needed();
    };
#endif
  }

  graph<vertex, W> G(V0, n, m, deletion_fn);
#ifdef NVM
  G.V0 = V0;
  G.V1 = V1;
#endif
  return G;
}
