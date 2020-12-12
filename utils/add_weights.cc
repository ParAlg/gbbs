// Usage: ./add_weights -enc <encoding> -o <output_file> input_graph
// Flags:
//   required:
//     -enc <oneof {adj, bytepd-amortized}>
//     -o <output file>
//   optional:
//     -s <if symmetric>
//     -c <if input file is compressed>
//     -m <if input graph should be mmaped>
//
// ex:
// > numactl -i all ./add_weights -enc bytepd-amortized -s -c -m -o
// hyperlink2014_sym_wgh.bytepda hyperlink2014_sym.bytepda

#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>

#include "gbbs/gbbs.h"

#include "pbbslib/random.h"
using namespace std;

namespace gbbs {
int max_weight = 32;
int* Choices;

template <class I>
struct wgh_it {
  uintE deg;
  I& it;
  uintE source;
  wgh_it(I& _it, uintE _src) : it(_it), source(_src), deg(0) {}
  tuple<uintE, int> cur() {
    auto cr = it.cur();
    uintE ind = pbbs::hash64(source ^ get<0>(cr));
    int wgh = Choices[ind % (2 * max_weight)];
    return make_tuple(get<0>(cr), wgh);
  }
  tuple<uintE, int> next() {
    deg++;
    auto nxt = it.next();
    uintE ind = pbbs::hash64(source ^ get<0>(nxt));
    int wgh = Choices[ind % (2 * max_weight)];
    return make_tuple(get<0>(nxt), wgh);
  }
  bool has_next() { return it.has_next(); }
};

template <class I, class GA>
struct degree_wgh_it {
  uintE deg;
  I& it;
  uintE source;
  uintE src_degree;
  GA& G;
  degree_wgh_it(I& _it, uintE _src, GA& G) : it(_it), source(_src), deg(0), G(G) {
    src_degree = G.V[source].getOutDegree();
  }
  tuple<uintE, int> cur() {
    auto cr = it.cur();
    uintE ngh = get<0>(cr);
    uintE ngh_degree = G.V[ngh].getOutDegree();
    int wgh = ngh_degree + src_degree; //(int)(((double)1)/((double)(ngh_degree + src_degree)));
    //uintE ind = pbbs::hash64(source ^ get<0>(cr));
    //int wgh = Choices[ind % (2 * max_weight)];
    return make_tuple(get<0>(cr), wgh);
  }
  tuple<uintE, int> next() {
    deg++;
    auto nxt = it.next();
    uintE ngh = get<0>(nxt);
    uintE ngh_degree = G.V[ngh].getOutDegree();
    int wgh = ngh_degree + src_degree; //(int)(((double)1)/((double)(ngh_degree + src_degree)));
//    uintE ind = pbbs::hash64(source ^ get<0>(nxt));
//    int wgh = Choices[ind % (2 * max_weight)];
    return make_tuple(get<0>(nxt), wgh);
  }
  bool has_next() { return it.has_next(); }
};

//template <class I, class GA>
//auto make_wgh_it(I& it, uintE source, GA& G) {
//  return wgh_it<I>(it, source);
//}

template <class I, class GA>
auto make_wgh_it(I& it, uintE source, GA& G) {
  return degree_wgh_it<I, GA>(it, source, G);
}

template <template <class W> class vertex, class W>
void writeWeightedAdj(graph<vertex<W>>& GA, string& outfile) {
  size_t n = GA.n;
  size_t m = GA.m;
  auto r = pbbs::random();

  auto degs = pbbs::sequence<uintT>(n + 1);
  par_for(0, n, [&] (size_t i) { degs[i] = GA.V[i].getOutDegree(); });
  degs[n] = 0;
  size_t total_offs = pbbs::scan_inplace(degs.slice(), pbbs::addm<uintT>());
  std::cout << "total offs = " << total_offs << " m = " << m << std::endl;

  auto edges = pbbs::sequence<uintT>(2 * m);
  uintT* wghs = edges.start() + m;

  parallel_for(0, n, [&] (size_t i) {
    size_t off = degs[i];
    size_t k = 0;
    auto itt = GA.V[i].getOutIter(i);
    auto it = make_wgh_it(itt, i, GA);
    if (itt.degree > 0) {
      while (true) {
        auto nghw = (k == 0) ? it.cur() : it.next();
        uintE ngh = get<0>(nghw);
        int weight = get<1>(nghw);
        edges[off + k] = ngh;
        wghs[off + k] = weight;
        k++;
        if (!it.has_next()) {
          break;
        }
      }
    }

    if (k != GA.V[i].getOutDegree()) {
      std::cout << "k = " << k << " deg = " << GA.V[i].getOutDegree() << std::endl;
    }
    assert(k == GA.V[i].getOutDegree());
  }, 1);

  ofstream file(outfile, ios::out | ios::binary);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << outfile << std::endl;
    exit(0);
  }
  file << "WeightedAdjacencyGraph" << std::endl;
  file << n << std::endl;
  file << m << std::endl;
  benchIO::writeArrayToStream(file, degs.start(), n);
  benchIO::writeArrayToStream(file, edges.start(), 2 * m);
  file.close();
  std::cout << "Wrote file." << std::endl;
  free(Choices);
}

namespace encodings {
namespace bytepd_amortized {
// Only writes the out-edges. Useful for comparing to other frameworks that
// run SSSP on a directed version of the graph.
template <template <class W> class vertex, class W>
void writeWeightedBytePDADirected(graph<vertex<W>>& GA, string& outfile) {
  size_t n = GA.n;

  ofstream out(outfile.c_str(), ofstream::out | ios::binary);
  // out-edges
  {
    auto degrees = pbbs::sequence<uintE>(n);
    auto byte_offsets = pbbs::sequence<uintT>(n + 1);
    auto r = pbbs::random();
    std::cout << "Calculating total size" << std::endl;
    pbbs::timer calc_t;
    calc_t.start();
    parallel_for(0, n, [&] (size_t i) {
      size_t total_bytes = 0;
      uintE last_ngh = 0;
      size_t deg = 0;
      uchar tmp[16];
      auto f = [&](const uintE& u, const uintE& v, const W& w) {
        long bytes = 0;
        uintE ind = pbbs::hash64(u ^ v);
        int weight = Choices[ind % (2 * max_weight)];
        if ((deg % PARALLEL_DEGREE) == 0) {
          bytes = compressFirstEdge(tmp, bytes, u, v);
          bytes = compressWeight<int32_t>(tmp, bytes, weight);
        } else {
          uintE diff = v - last_ngh;
          bytes = compressEdge(tmp, bytes, diff);
          bytes = compressWeight<int32_t>(tmp, bytes, weight);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
        return false;
      };

      if (GA.V[i].getOutDegree() > 0) {
        GA.V[i].mapOutNgh(i, f, false);
        assert(deg == GA.V[i].getOutDegree());

        size_t n_chunks = 1 + (deg - 1) / PARALLEL_DEGREE;
        // To account for the byte offsets
        total_bytes += (n_chunks - 1) * sizeof(uintE);
        // To account for the per-block counters
        total_bytes += (n_chunks) * sizeof(uintE);
        // To account for the virtual degree
        total_bytes += sizeof(uintE);
      }

      degrees[i] = deg;
      byte_offsets[i] = total_bytes;
    }, 1);
    calc_t.stop();
    calc_t.reportTotal("total size time");

    byte_offsets[n] = 0;
    size_t total_space = pbbs::scan_add(byte_offsets, byte_offsets);
    std::cout << "total space is: " << total_space << std::endl;
    auto edges = pbbs::sequence<uchar>(total_space);
    std::cout << "Allocated, compressing!" << std::endl;

    parallel_for(0, n, [&] (size_t i) {
      uintE deg = degrees[i];
      if (deg > 0) {
        auto iter = GA.V[i].getOutIter(i);
        auto it = make_wgh_it(iter, i, GA);
        long nbytes = sequentialCompressEdgeSet<int32_t>(
            edges.start() + byte_offsets[i], 0, deg, (uintE)i, it);
        if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
          std::cout << "nbytes = " << nbytes
               << " but offs = " << (byte_offsets[i + 1] - byte_offsets[i])
               << " deg = " << deg << " i = " << i << std::endl;
        }
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    }, 1);
    std::cout << "Compressed" << std::endl;

    long* sizes = pbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = total_space;
    out.write((char*)sizes, sizeof(long) * 3);  // write n, m and space used
    out.write((char*)byte_offsets.start(),
              sizeof(uintT) * (n + 1));  // write offsets
    out.write((char*)degrees.start(), sizeof(uintE) * n);
    out.write((char*)edges.start(), total_space);  // write edges
    free(sizes);
  }

  out.close();
}

template <template <class W> class vertex, class W>
void writeWeightedBytePDA(graph<vertex<W>>& GA, string& outfile,
                          bool symmetric) {
  if (!symmetric) {
    writeWeightedBytePDADirected(GA, outfile);
    return;
  }
  size_t n = GA.n;

  auto degrees = pbbs::sequence<uintE>(n);
  auto byte_offsets = pbbs::sequence<uintT>(n + 1);
  auto r = pbbs::random();
  std::cout << "Calculating total size" << std::endl;
  pbbs::timer calc_t;
  calc_t.start();
  parallel_for(0, n, [&] (size_t i) {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    uchar tmp[16];
    auto f = [&](const uintE& u, const uintE& v, const W& w) {
      long bytes = 0;
      uintE ind = pbbs::hash64(u ^ v);
      int weight = Choices[ind % (2 * max_weight)];
      if ((deg % PARALLEL_DEGREE) == 0) {
        bytes = compressFirstEdge(tmp, bytes, u, v);
        bytes = compressWeight<int32_t>(tmp, bytes, weight);
      } else {
        uintE diff = v - last_ngh;
        bytes = compressEdge(tmp, bytes, diff);
        bytes = compressWeight<int32_t>(tmp, bytes, weight);
      }
      last_ngh = v;
      total_bytes += bytes;
      deg++;
      return false;
    };

    if (GA.V[i].getOutDegree() > 0) {
      GA.V[i].mapOutNgh(i, f, false);
      assert(deg == GA.V[i].getOutDegree());

      size_t n_chunks = 1 + (deg - 1) / PARALLEL_DEGREE;
      // To account for the byte offsets
      total_bytes += (n_chunks - 1) * sizeof(uintE);
      // To account for the per-block counters
      total_bytes += (n_chunks) * sizeof(uintE);
      // To account for the virtual degree
      total_bytes += sizeof(uintE);
    }

    degrees[i] = deg;
    byte_offsets[i] = total_bytes;
  }, 1);
  calc_t.stop();
  calc_t.reportTotal("total size time");

  byte_offsets[n] = 0;
  size_t total_space = pbbs::scan_add(byte_offsets, byte_offsets);
  std::cout << "total space is: " << total_space << std::endl;
  auto edges = pbbs::sequence<uchar>(total_space);

  parallel_for(0, n, [&] (size_t i) {
    uintE deg = degrees[i];
    if (deg > 0) {
      auto iter = GA.V[i].getOutIter(i);
      auto it = make_wgh_it(iter, i, GA);
      long nbytes = sequentialCompressEdgeSet<int32_t>(
          edges.start() + byte_offsets[i], 0, deg, (uintE)i, it);
      if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
        std::cout << "nbytes = " << nbytes
             << " but offs = " << (byte_offsets[i + 1] - byte_offsets[i])
             << " deg = " << deg << " i = " << i << std::endl;
      }
      assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
    }
  }, 1);
  std::cout << "Compressed" << std::endl;

  long* sizes = pbbs::new_array_no_init<long>(3);
  sizes[0] = GA.n;
  sizes[1] = GA.m;
  sizes[2] = total_space;
  ofstream out(outfile.c_str(), ofstream::out | ios::binary);
  out.write((char*)sizes, sizeof(long) * 3);  // write n, m and space used
  out.write((char*)byte_offsets.start(),
            sizeof(uintT) * (n + 1));  // write offsets
  out.write((char*)degrees.start(), sizeof(uintE) * n);
  out.write((char*)edges.start(), total_space);  // write edges
  out.close();
}
};  // namespace bytepd_amortized
};  // namespace encodings

template <class vertex>
void Reencoder(graph<vertex>& GA, commandLine P) {
  auto outfile =
      P.getOptionValue("-o", "/ssd0/graphs/bench_experiments/out2.adj");
  auto encoding = P.getOptionValue("-enc", "adj");
  bool symmetric = P.getOptionValue("-s");
  bool unit_weights = P.getOptionValue("-unit");

  if (!unit_weights) {
    max_weight = log2(GA.n);
    int* Choices = pbbs::new_array_no_init<int>(2*max_weight);
    for (int i = 0; i < max_weight; i++) {
      Choices[2 * i] = i + 1;
      Choices[2 * i + 1] = i + 1;
    }
    std::cout << "mw = " << max_weight << std::endl;
  } else {
    max_weight = 1;
    int* Choices = pbbs::new_array_no_init<int>(2*max_weight);
    Choices[0] = 1; Choices[1] = 1;
  }

  if (encoding == "adj") {
    writeWeightedAdj(GA, outfile);
  } else if (encoding == "bytepd-amortized") {
    encodings::bytepd_amortized::writeWeightedBytePDA(GA, outfile, symmetric);
  }
  std::cout << "wrote output graph to: " << outfile << std::endl;
  // prevent running multiple times if -rounds 1 is not specified
  exit(0);
}
}  // namespace gbbs

generate_main(gbbs::Reencoder, false);
