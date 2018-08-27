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

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "graph.h"
#include "oldlib/block_radix_sort.h"
#include "oldlib/utils.h"

using namespace std;

typedef pair<uintE,uintE> intPair;
typedef pair<uintE, pair<uintE,intE> > intTriple;

template <class E>
struct pairFirstCmp {
  bool operator() (pair<uintE,E> a, pair<uintE,E> b) {
    return a.first < b.first; }
};

template <class E>
struct getFirst {uintE operator() (pair<uintE,E> a) {return a.first;} };

template <class IntType>
struct pairBothCmp {
  bool operator() (pair<uintE,IntType> a, pair<uintE,IntType> b) {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  }
};

// A structure that keeps a seq of strings all allocated from
// the same block of memory
struct words {
  long n; // total number of characters
  char* Chars;  // array storing all strings
  long m; // number of substrings
  char** Strings; // pointers to strings (all should be null terminated)
  words() {}
words(char* C, long nn, char** S, long mm)
: Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {free(Chars); free(Strings);}
};

inline bool isSpace(char c) {
  switch (c)  {
  case '\r':
  case '\t':
  case '\n':
  case 0:
  case ' ' : return true;
  default : return false;
  }
}

ligra_utils::_seq<char> mmapStringFromFile(const char *filename) {
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
  if (!S_ISREG (sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char *p = static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
//  char *bytes = newA(char, n);
//  parallel_for(size_t i=0; i<n; i++) {
//    bytes[i] = p[i];
//  }
//  if (munmap(p, sb.st_size) == -1) {
//    perror("munmap");
//    exit(-1);
//  }
//  cout << "mmapped" << endl;
//  free(bytes);
//  exit(0);
  return ligra_utils::_seq<char>(p, n);
}

ligra_utils::_seq<char> readStringFromFile(char *fileName) {
  ifstream file (fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long end = file.tellg();
  file.seekg (0, ios::beg);
  long n = end - file.tellg();
  char* bytes = newA(char,n+1);
  file.read (bytes,n);
  file.close();
  return ligra_utils::_seq<char>(bytes,n);
}

// parallel code for converting a string to words
words stringToWords(char *Str, long n) {
  {parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      if (isSpace(Str[i])) Str[i] = 0; }); }

  // mark start of words
  bool *FL = newA(bool,n);
  FL[0] = Str[0];
  {
    parallel_for_bc(i, 1, n, (n > pbbs::kSequentialForThreshold), { FL[i] = Str[i] && !Str[i-1];});
  }

//  cout << "n (strlen) = " << n << endl;
//  auto im = make_in_imap<size_t>(n, [&] (size_t i) { return FL[i]; });
//  cout << " sum is : " << pbbs::reduce_add(im) << endl;

  // offset for each start of word
  ligra_utils::_seq<long> Off = ligra_utils::seq::packIndex<long>(FL, n);
//  cout << "pack returned " << Off.n << endl;
  long m = Off.n;
  long *offsets = Off.A;

  // pointer to each start of word
  char **SA = newA(char*, m);
  {
    parallel_for_bc(j, 0, m, (m > pbbs::kSequentialForThreshold), {SA[j] = Str+offsets[j];});
  }

  free(offsets); free(FL);
  return words(Str,n,SA,m);
}





template <template <typename W> class vertex>
auto readWeightedGraph(char* fname, bool isSymmetric, bool mmap) {
  using wvtx = vertex<intE>;
  words W;
  if (mmap) {
    ligra_utils::_seq<char> S = mmapStringFromFile(fname);
    char *bytes = newA(char, S.n);
    // Cannot mutate the graph unless we copy.
    parallel_for_bc(i, 0, S.n, (S.n > pbbs::kSequentialForThreshold), {
      bytes[i] = S.A[i];
    });
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  } else {
    ligra_utils::_seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n);
  }
  assert(W.Strings[0] == (string) "WeightedAdjacencyGraph");

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);
  assert(len == n + 2*m + 2);

  uintT* offsets = newA(uintT,n);
  using VW = tuple<uintE, intE>;
  tuple<uintE, intE>* edges = newA(VW,2*m);

  {parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), { offsets[i] = atol(W.Strings[i + 3]);}); }
  {parallel_for_bc(i, 0, m, (m > pbbs::kSequentialForThreshold),{
    edges[i] = make_tuple(atol(W.Strings[i+n+3]), atol(W.Strings[i+n+m+3])); });
  }
  //W.del(); // to deal with performance bug in malloc

  wvtx* v = newA(wvtx,n);

  {
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    uintT o = offsets[i];
    uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].setOutDegree(l);
    v[i].setOutNeighbors(edges+o); });
  }

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), { tOffsets[i] = INT_T_MAX;});
    intTriple* temp = newA(intTriple,m);
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      uintT o = offsets[i];
      for(uintT j=0;j<v[i].getOutDegree();j++) {
        temp[o+j] = make_pair(v[i].getOutNeighbor(j),make_pair(i,v[i].getOutWeight(j)));
      }
    });
    free(offsets);

    intSort::iSort(temp,m,n+1,getFirst<intPair>());

    tOffsets[temp[0].first] = 0;
    VW* inEdges = newA(VW,m);
    inEdges[0] = make_tuple(temp[0].second.first, temp[0].second.second);

    parallel_for_bc(i, 1, m, (m > pbbs::kSequentialForThreshold), {
      inEdges[i] = make_tuple(temp[i].second.first, temp[i].second.second);
      if(temp[i].first != temp[i-1].first) {
      	tOffsets[temp[i].first] = i;
      }
    });

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    ligra_utils::seq::scanIBack(tOffsets,tOffsets,n,ligra_utils::minF<uintT>(),(uintT)m);

    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
      v[i].setInNeighbors(inEdges+o);
    });

    free(tOffsets);
    return graph<wvtx>(v,n,m,get_deletion_fn(v, inEdges, edges), get_copy_fn(v, inEdges, edges, n, m, m, m));
  }
  else {
    free(offsets);
    return graph<wvtx>(v,n,m,get_deletion_fn(v, edges), get_copy_fn(v, edges, n, m, m));
  }
}

template <template <typename W> class vertex>
auto readUnweightedGraph(char* fname, bool isSymmetric, bool mmap) {
  using wvtx = vertex<pbbs::empty>;
  words W;
  if (mmap) {
    ligra_utils::_seq<char> S = mmapStringFromFile(fname);
    char *bytes = newA(char, S.n);
    // Cannot mutate the graph unless we copy.
    parallel_for_bc(i, 0, S.n, (S.n > pbbs::kSequentialForThreshold), {
      bytes[i] = S.A[i];
    });
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  } else {
    ligra_utils::_seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n);
  }
  assert(W.Strings[0] == (string) "AdjacencyGraph");

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);

  cout << "n = " << n << " m = " << m << " len = " << len << endl;
  assert(len == n + m + 2);

  uintT* offsets = newA(uintT,n);
  uintE* edges = newA(uintE,m);

  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    offsets[i] = atol(W.Strings[i + 3]);});
  parallel_for_bc(i, 0, m, (m > pbbs::kSequentialForThreshold), {
    edges[i] = atol(W.Strings[i+n+3]);
  });
  //W.del(); // to deal with performance bug in malloc

  wvtx* v = newA(wvtx,n);

  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    uintT o = offsets[i];
    uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].setOutDegree(l);
    v[i].setOutNeighbors(((tuple<uintE, pbbs::empty>*)(edges+o)));
  });

  if(!isSymmetric) {
    uintT* tOffsets = newA(uintT,n);
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      tOffsets[i] = INT_T_MAX;});
    intPair* temp = newA(intPair,m);
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      uintT o = offsets[i];
      for(uintT j=0;j<v[i].getOutDegree();j++){
        temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
      }
    });
    free(offsets);

    intSort::iSort(temp,m,n+1,getFirst<uintE>());

    tOffsets[temp[0].first] = 0;
    uintE* inEdges = newA(uintE,m);
    inEdges[0] = temp[0].second;
    parallel_for_bc(i, 1, m, (m > pbbs::kSequentialForThreshold), {
      inEdges[i] = temp[i].second;
      if(temp[i].first != temp[i-1].first) {
        tOffsets[temp[i].first] = i;
      }
    });

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    ligra_utils::seq::scanIBack(tOffsets,tOffsets,n,ligra_utils::minF<uintT>(),(uintT)m);

    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      uintT o = tOffsets[i];
      uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
      v[i].setInDegree(l);
      v[i].setInNeighbors((tuple<uintE, pbbs::empty>*)(inEdges+o));
    });

    free(tOffsets);
    return graph<wvtx>(v,n,m,get_deletion_fn(v, inEdges, edges),
        get_copy_fn(v, (tuple<uintE, pbbs::empty>*)inEdges, (tuple<uintE, pbbs::empty>*)edges, n, m, m, m));
  }
  else {
    free(offsets);
    return graph<wvtx>(v,n,m,get_deletion_fn(v, edges), get_copy_fn(v, (tuple<uintE, pbbs::empty>*)edges, n, m, m));
  }
}

template <class W, typename std::enable_if<!std::is_same<W, intE>::value, int>::type=0>
string print_wgh(W wgh) {
  return "";
}

template <class W, typename std::enable_if<std::is_same<W, intE>::value, int>::type=0>
string print_wgh(W wgh) {
  return std::to_string(wgh);
}

template <template <typename W> class vertex, class W>
graph<vertex<W> > readCompressedGraph(char* fname, bool isSymmetric, bool mmap, bool mmapcopy) {
  char* s;
  using w_vertex = vertex<W>;
  if (mmap) {
    ligra_utils::_seq<char> S = mmapStringFromFile(fname);
    s = S.A;
    if (mmapcopy) {
      cout << "Copying compressed graph" << endl;
      // Cannot mutate graph unless we copy.
      char *bytes = newA(char, S.n);
      parallel_for_bc(i, 0, S.n, (S.n > pbbs::kSequentialForThreshold), {
        bytes[i] = S.A[i];
      });
      if (munmap(S.A, S.n) == -1) {
        perror("munmap");
        exit(-1);
      }
      s = bytes;
    }
  } else {

    int fd;
    if ( (fd = open(fname, O_RDONLY | O_DIRECT) ) != -1) {
      cout << "input opened!" << endl;
    } else {
      cout << "can't open input file!";
    }
//    posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);

    size_t fsize = lseek(fd, 0, SEEK_END);
    lseek(fd, 0, 0);
    s = (char*)memalign(4096 * 2, fsize + 4096);

    cout << "fsize = " << fsize << endl;

    size_t sz = 0;

    size_t pgsize = getpagesize();
    cout << "pgsize = " << pgsize << endl;

    size_t read_size = 1024*1024*1024;
    if (sz + read_size > fsize) {
      size_t k = std::ceil((fsize - sz) / pgsize);
      read_size = std::max(k*pgsize, pgsize);
      cout << "set read size to: " << read_size << " " << (fsize - sz) << " bytes left" << endl;
    }

    while (sz + read_size < fsize) {
      void* buf = s + sz;
      cout << "reading: " << read_size << endl;
      sz += read(fd, buf, read_size);
      cout << "read: " << sz << " bytes" << endl;
      if (sz + read_size > fsize) {
        size_t k = std::ceil((fsize - sz) / pgsize);
        read_size = std::max(k*pgsize, pgsize);
        cout << "set read size to: " << read_size << " " << (fsize - sz) << " bytes left" << endl;
      }
    }
    if (sz < fsize) {
      cout << "last read: rem = " << (fsize - sz) << endl;
      void* buf = s + sz;
      sz += read(fd, buf, pgsize);
      cout << "read " << sz << " bytes " << endl;
    }

//    while (sz < fsize) {
//      size_t rem = fsize - sz;
//      size_t read_size = pgsize;
////      size_t read_size = std::min(pgsize, rem);
//      void* buf = s + sz;
//      cout << "reading: " << read_size << endl;
//      sz += read(fd, buf, read_size);
//      cout << "read: " << sz << " bytes" << endl;
//    }
//    cout << "Finished: read " << sz << " out of fsize = " << fsize << endl;
    close(fd);

//    ifstream in(fname,ifstream::in |ios::binary);
//    in.seekg(0,ios::end);
//    long size = in.tellg();
//    in.seekg(0);
//    cout << "size = " << size << endl;
//    s = (char*) malloc(size);
//    in.read(s,size);
//    cout << "Finished read" << endl;
//    in.close();
  }

  long* sizes = (long*) s;
  long n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  cout << "n = "<<n<<" m = "<<m<<" totalSpace = "<<totalSpace<<endl;
  cout << "reading file..."<<endl;

  uintT* offsets = (uintT*) (s+3*sizeof(long));
  long skip = 3*sizeof(long) + (n+1)*sizeof(intT);
  uintE* Degrees = (uintE*) (s+skip);
  skip+= n*sizeof(intE);
  uchar* edges = (uchar*)(s+skip);

  uintT* inOffsets;
  uchar* inEdges;
  uintE* inDegrees;
  long inTotalSpace = 0;
  if(!isSymmetric){
    skip += totalSpace;
    uchar* inData = (uchar*)(s + skip);
    sizes = (long*) inData;
    inTotalSpace = sizes[0];
    cout << "inTotalSpace = "<<inTotalSpace<<endl;
    skip += sizeof(long);
    inOffsets = (uintT*) (s + skip);
    skip += (n+1)*sizeof(uintT);
    inDegrees = (uintE*)(s+skip);
    skip += n*sizeof(uintE);
    inEdges = (uchar*)(s + skip);
  } else {
    inOffsets = offsets;
    inEdges = edges;
    inDegrees = Degrees;
  }

  w_vertex *V = newA(w_vertex,n);
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    long o = offsets[i];
    uintT d = Degrees[i];
    V[i].setOutDegree(d);
    V[i].setOutNeighbors(edges+o);
  });

  if(!isSymmetric) {
    parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
      long o = inOffsets[i];
      uintT d = inDegrees[i];
      V[i].setInDegree(d);
      V[i].setInNeighbors(inEdges+o);
    });
    graph<w_vertex> G(V,n,m,get_deletion_fn(V, s), get_copy_fn(V, inEdges, edges, n, m, totalSpace, inTotalSpace));
    return G;
  } else {
    graph<w_vertex> G(V,n,m,get_deletion_fn(V, s), get_copy_fn(V, edges, n, m, totalSpace));
    return G;
  }
}
