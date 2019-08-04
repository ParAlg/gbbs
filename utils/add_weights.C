// Usage: ./add_weights -enc <encoding> -o <output_file> input_symmetric_graph
// Flags:
//   required:
//     -enc <oneof {adj, bytepd-amortized}>
//     -o <output file>
//   optional:
//     -s <if symmetric>
//     -c <if input file is compressed>
//     -m <if input symmetric_graph should be mmaped>
//
// ex:
// > numactl -i all ./add_weights -enc bytepd-amortized -s -c -m -o
// hyperlink2014_sym_wgh.bytepda hyperlink2014_sym.bytepda

#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>

#include "ligra.h"

#include "pbbslib/random.h"
#include "encodings/byte_pd_amortized.h"
using namespace std;

int max_weight = 32;
int* Choices;
size_t PAR_DEGREE_TWO = 256;

template <class I>
struct wgh_it {
  uintE deg;
  I& it;
  uintE source;
  size_t src_shift;
  wgh_it(I& _it, uintE _src) : it(_it), source(_src), deg(0) {
    src_shift = static_cast<size_t>(source) << 32L;
  }
  tuple<uintE, int> cur() {
    auto cr = it.cur();
    uintE ind = pbbs::hash64(src_shift ^ get<0>(cr));
    int wgh = Choices[ind % (2 * max_weight)];
    return make_tuple(get<0>(cr), wgh);
  }
  tuple<uintE, int> next() {
    deg++;
    auto nxt = it.next();
    uintE ind = pbbs::hash64(src_shift ^ get<0>(nxt));
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
    src_degree = G.get_vertex(source).getOutDegree();
  }
  tuple<uintE, int> cur() {
    auto cr = it.cur();
    uintE ngh = get<0>(cr);
    uintE ngh_degree = G.get_vertex(ngh).getOutDegree();
    int wgh = ngh_degree + src_degree; //(int)(((double)1)/((double)(ngh_degree + src_degree)));
    //uintE ind = pbbs::hash64(source ^ get<0>(cr));
    //int wgh = Choices[ind % (2 * max_weight)];
    return make_tuple(get<0>(cr), wgh);
  }
  tuple<uintE, int> next() {
    deg++;
    auto nxt = it.next();
    uintE ngh = get<0>(nxt);
    uintE ngh_degree = G.get_vertex(ngh).getOutDegree();
    int wgh = ngh_degree + src_degree; //(int)(((double)1)/((double)(ngh_degree + src_degree)));
//    uintE ind = pbbs::hash64(source ^ get<0>(nxt));
//    int wgh = Choices[ind % (2 * max_weight)];
    return make_tuple(get<0>(nxt), wgh);
  }
  bool has_next() { return it.has_next(); }
};

template <class I, class GA>
auto make_wgh_it(I& it, uintE source, GA& G) {
  return wgh_it<I>(it, source);
}

//template <template <class W> class vertex, class W>
//void writeWeightedAdj(symmetric_graph<vertex, W>& GA, string& outfile) {
//  size_t n = GA.n;
//  size_t m = GA.m;
//  auto r = pbbs::random();
//
//  auto degs = pbbs::sequence<uintT>(n + 1);
//  par_for(0, n, [&] (size_t i) { degs[i] = GA.V[i].getOutDegree(); });
//  degs[n] = 0;
//  size_t total_offs = pbbs::scan_inplace(degs.slice(), pbbs::addm<uintT>());
//  cout << "total offs = " << total_offs << " m = " << m << endl;
//
//  auto edges = pbbs::sequence<uintT>(2 * m);
//  uintT* wghs = edges.begin() + m;
//
//  parallel_for(0, n, [&] (size_t i) {
//    size_t off = degs[i];
//    size_t k = 0;
//    auto itt = GA.V[i].getOutIter(i);
//    auto it = make_wgh_it(itt, i, GA);
//    if (itt.degree > 0) {
//      while (true) {
//        auto nghw = (k == 0) ? it.cur() : it.next();
//        uintE ngh = get<0>(nghw);
//        int weight = get<1>(nghw);
//        edges[off + k] = ngh;
//        wghs[off + k] = weight;
//        k++;
//        if (!it.has_next()) {
//          break;
//        }
//      }
//    }
//
//    if (k != GA.V[i].getOutDegree()) {
//      cout << "k = " << k << " deg = " << GA.V[i].getOutDegree() << endl;
//    }
//    assert(k == GA.V[i].getOutDegree());
//  }, 1);
//
//  ofstream file(outfile, ios::out | ios::binary);
//  if (!file.is_open()) {
//    std::cout << "Unable to open file: " << outfile << std::endl;
//    exit(0);
//  }
//  file << "WeightedAdjacencyGraph" << endl;
//  file << n << endl;
//  file << m << endl;
//  benchIO::writeArrayToStream(file, degs.begin(), n);
//  benchIO::writeArrayToStream(file, edges.begin(), 2 * m);
//  file.close();
//  cout << "Wrote file." << endl;
//  free(Choices);
//}

namespace bytepd_amortized {
// Only writes the out-edges. Useful for comparing to other frameworks that
// run SSSP on a directed version of the symmetric_graph.
template <template <class W> class vertex, class W>
void writeWeightedBytePDADirected(symmetric_graph<vertex, W>& GA, string& outfile) {
  size_t n = GA.n;

  ofstream out(outfile.c_str(), ofstream::out | ios::binary);
  // out-edges
  {
    auto degrees = pbbs::sequence<uintE>(n);
    auto byte_offsets = pbbs::sequence<uintT>(n + 1);
    auto r = pbbs::random();
    cout << "Calculating total size" << endl;
    timer calc_t;
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
        if ((deg % PAR_DEGREE_TWO) == 0) {
          bytes = bytepd_amortized::compressFirstEdge(tmp, bytes, u, v);
          bytes = bytepd_amortized::compressWeight<int32_t>(tmp, bytes, weight);
        } else {
          uintE diff = v - last_ngh;
          bytes = bytepd_amortized::compressEdge(tmp, bytes, diff);
          bytes = bytepd_amortized::compressWeight<int32_t>(tmp, bytes, weight);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
        return false;
      };

      if (GA.get_vertex(i).getOutDegree() > 0) {
        GA.get_vertex(i).mapOutNgh(i, f, false);
        assert(deg == GA.get_vertex(i).getOutDegree());

        size_t n_chunks = 1 + (deg - 1) / PAR_DEGREE_TWO;
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
    size_t total_space = pbbslib::scan_add_inplace(byte_offsets.slice());
    cout << "total space is: " << total_space << endl;
    auto edges = pbbs::sequence<uchar>(total_space);
    cout << "Allocated, compressing!" << endl;

    parallel_for(0, n, [&] (size_t i) {
      uintE deg = degrees[i];
      if (deg > 0) {
        auto iter = GA.get_vertex(i).getOutIter(i);
        auto it = make_wgh_it(iter, i, GA);
        long nbytes = sequentialCompressEdgeSet<int32_t>(
            edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
        if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
          cout << "nbytes = " << nbytes
               << " but offs = " << (byte_offsets[i + 1] - byte_offsets[i])
               << " deg = " << deg << " i = " << i << endl;
        }
        size_t bytes_next = byte_offsets[i+1];
        size_t bytes_us = byte_offsets[i];
        size_t bytes_expec = bytes_next - bytes_us;
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    }, 1);
    cout << "Compressed" << endl;

    long* sizes = pbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = total_space;
    out.write((char*)sizes, sizeof(long) * 3);  // write n, m and space used
    out.write((char*)byte_offsets.begin(),
              sizeof(uintT) * (n + 1));  // write offsets
    out.write((char*)degrees.begin(), sizeof(uintE) * n);
    out.write((char*)edges.begin(), total_space);  // write edges
    free(sizes);
  }

  out.close();
}

template <template <class W> class vertex, class W>
void writeWeightedBytePDA(symmetric_graph<vertex, W>& GA, string& outfile,
                          bool symmetric, size_t n_batches=5) {
  assert(symmetric);
  size_t n = GA.n;

  auto degrees = pbbs::sequence<uintE>(n);
  auto byte_offsets = pbbs::sequence<uintT>(n + 1);
  auto r = pbbs::random();
  cout << "Calculating total size" << endl;
  timer calc_t;
  calc_t.start();
  parallel_for(0, n, [&] (size_t i) {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    size_t src_shift = static_cast<size_t>(i) << 32L;
    uchar tmp[16];
    auto f = [&](const uintE& u, const uintE& v, const W& w) {
      long bytes = 0;
      uintE ind = pbbs::hash64(src_shift | static_cast<size_t>(v));
      int weight = Choices[ind % (2 * max_weight)];
      if ((deg % PAR_DEGREE_TWO) == 0) {
        bytes = bytepd_amortized::compressFirstEdge(tmp, bytes, u, v);
        bytes = bytepd_amortized::compressWeight<int32_t>(tmp, bytes, weight);
      } else {
        uintE diff = v - last_ngh;
        bytes = bytepd_amortized::compressEdge(tmp, bytes, diff);
        bytes = bytepd_amortized::compressWeight<int32_t>(tmp, bytes, weight);
      }
      last_ngh = v;
      total_bytes += bytes;
      deg++;
      return false;
    };

    if (GA.get_vertex(i).getOutDegree() > 0) {
      GA.get_vertex(i).mapOutNgh(i, f, false);
      assert(deg == GA.get_vertex(i).getOutDegree());

      size_t n_chunks = 1 + (deg - 1) / PAR_DEGREE_TWO;
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
  size_t total_space = pbbslib::scan_add_inplace(byte_offsets.slice());

  long* sizes = pbbs::new_array_no_init<long>(3);
  sizes[0] = GA.n;
  sizes[1] = GA.m;
  sizes[2] = total_space;
  ofstream out(outfile.c_str(), ofstream::out | ios::binary);
  out.write((char*)sizes, sizeof(long) * 3);  // write n, m and space used
  out.write((char*)byte_offsets.begin(),
            sizeof(uintT) * (n + 1));  // write offsets
  out.write((char*)degrees.begin(), sizeof(uintE) * n);


  size_t block_size = pbbs::num_blocks(n, n_batches);
  size_t edges_written = 0;
  for (size_t b=0; b<n_batches; b++) {
    size_t start = b*block_size;
    size_t end = std::min(start + block_size, n);
    if (start >= end) break;
    cout << "writing vertices " << start << " to " << end << endl;

    // create slab of graph and write out
    size_t start_offset = byte_offsets[start];
    size_t end_offset = byte_offsets[end];
    size_t n_edge_bytes = end_offset - start_offset;
    auto edges = pbbs::new_array_no_init<uchar>(n_edge_bytes);

    parallel_for(start, end, [&] (size_t i) {
      uintE deg = degrees[i];
      if (deg > 0) {
        auto iter = GA.get_vertex(i).getOutIter(i);
        auto it = make_wgh_it(iter, i, GA);
        size_t our_offset = byte_offsets[i] - start_offset;
        long nbytes = sequentialCompressEdgeSet<int32_t>(
            edges + our_offset, 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
        if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
          cout << "nbytes = " << nbytes
               << " but offs = " << (byte_offsets[i + 1] - byte_offsets[i])
               << " deg = " << deg << " i = " << i << endl;
        }
       size_t bytes_next = byte_offsets[i+1];
       size_t bytes_us = byte_offsets[i];
       size_t bytes_expec = bytes_next - bytes_us;
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    }, 1);
    cout << "Compressed" << endl;
    out.write((char*)edges, n_edge_bytes);  // write edges
    cout << "Wrote" << endl;
    pbbs::free_array(edges);
  }
  out.close();
  cout << "Finished writing file" << endl;
}
};  // namespace bytepd_amortized

namespace binary_format {

  template <class G>
  void add_weights_binary_format(G& GA, string& outfile, bool symmetric, size_t n_batches=4) {
    size_t n = GA.n; size_t m = GA.m;
    using W = intE;
    using edge_type = std::tuple<uintE, W>;

    ofstream out(outfile.c_str(), ofstream::out | ios::binary);

    // 1. Calculate total size
    auto offsets = pbbs::sequence<uintT>(n+1);
    cout << "calculating size" << endl;
    parallel_for(0, n, [&] (size_t i) {
      offsets[i] = GA.get_vertex(i).getOutDegree();
    });
    offsets[n] = 0;
    size_t offset_scan = pbbslib::scan_add_inplace(offsets.slice());
    cout << "offset_scan = " << offset_scan << " m = " << m << endl;
    assert(offset_scan == m);

    long* sizes = pbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = sizeof(long) * 3 + sizeof(uintT)*(n+1) + sizeof(edge_type)*m;
    out.write((char*)sizes,sizeof(long)*3); //write n, m and space used
    out.write((char*)offsets.begin(),sizeof(uintT)*(n+1)); //write offsets

    // 2. Create compressed format in-memory (batched)
    size_t block_size = pbbs::num_blocks(n, n_batches);
    size_t edges_written = 0;
    for (size_t b=0; b<n_batches; b++) {
      size_t start = b*block_size;
      size_t end = std::min(start + block_size, n);
      if (start >= end) break;
      cout << "writing vertices " << start << " to " << end << endl;

      // create slab of graph and write out
      size_t start_offset = offsets[start];
      size_t end_offset = offsets[end];
      size_t n_edges = end_offset - start_offset;
      edge_type* edges = pbbs::new_array_no_init<edge_type>(n_edges);

      parallel_for(start, end, [&] (size_t i) {
        size_t our_offset = offsets[i] - start_offset;
        size_t src_shift = static_cast<size_t>(i) << 32L;
        auto map_f = [&] (const uintE& u, const uintE& v, const pbbs::empty& wgh) {
          return wgh;
        };
        auto write_f = [&] (const uintE& ngh, const uintT& offset, const pbbs::empty& val) {
          uintE ind = pbbs::hash64(src_shift | static_cast<size_t>(ngh));
          int weight = Choices[ind % (2 * max_weight)];
          edges[offset] = std::make_tuple(ngh, weight);
        };
        GA.get_vertex(i).copyOutNgh(i, our_offset, map_f, write_f);
      });

      size_t edge_space = sizeof(edge_type)*n_edges;
      out.write((char*)edges,edge_space); //write edges

      cout << "finished writing vertices " << start << " to " << end << endl;
      edges_written += n_edges;

      pbbs::free_array(edges);
    }
    cout << "Wrote " << edges_written << " edges in total, m = " << m << endl;
    assert(edges_written == m);
    out.close();
  }


}; // namespace binary_format

template <class G>
double Reencoder(G& GA, commandLine P) {
  auto outfile =
      P.getOptionValue("-o", "/ssd0/symmetric_graphs/bench_experiments/out2.adj");
  auto encoding = P.getOptionValue("-enc", "adj");
  bool symmetric = P.getOptionValue("-s");
  bool unit_weights = P.getOptionValue("-unit");
  PAR_DEGREE_TWO = P.getOptionLongValue("-bs", 256);

  if (!unit_weights) {
    max_weight = log2(GA.n);
    Choices = pbbs::new_array_no_init<int>(2*max_weight);
    for (int i = 0; i < max_weight; i++) {
      Choices[2 * i] = i + 1;
      Choices[2 * i + 1] = i + 1;
    }
    cout << "mw = " << max_weight << endl;
  } else {
    max_weight = 1;
    int* Choices = pbbs::new_array_no_init<int>(2*max_weight);
    Choices[0] = 1; Choices[1] = 1;
    for (int i = 0; i < max_weight; i++) {
      Choices[2 * i] = 1;
      Choices[2 * i + 1] = 1;
    }
  }

  if (encoding == "adj") {
//    writeWeightedAdj(GA, outfile);
  } else if (encoding == "bytepd-amortized") {
    bytepd_amortized::writeWeightedBytePDA(GA, outfile, symmetric);
  } else if (encoding == "binary") {
    binary_format::add_weights_binary_format(GA, outfile, symmetric);
  }
  cout << "wrote output symmetric_graph to: " << outfile << endl;
  // prevent running multiple times if -rounds 1 is not specified
  exit(0);
  return 1.0;
}

generate_main(Reencoder, false);
