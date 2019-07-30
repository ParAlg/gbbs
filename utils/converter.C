#include "ligra.h"
#include "IO.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

using namespace std;

// Provides utilities for converting between different compressed
// representations.
constexpr int max_weight = 32;
#define PAR_DEGREE_TWO 256

namespace byte {
  template <template <class W> class vertex, class W>
  void write_graph_byte_format(symmetric_graph<vertex, W>& GA, ofstream& out, bool symmetric) {
    size_t n = GA.n; size_t m = GA.m;

    // 1. Calculate total size
    auto degrees = pbbs::sequence<uintE>(n);
    auto byte_offsets = pbbs::sequence<uintT>(n+1);
    cout << "calculating size" << endl;
    parallel_for(0, n, [&] (size_t i) {
      size_t total_bytes = 0;
      uintE last_ngh = 0;
      size_t deg = 0;
      uchar tmp[16];
      auto f = [&] (uintE u, uintE v, W w) {
        long bytes = 0;
        if (deg == 0) {
          bytes = compressFirstEdge(tmp, bytes, u, v);
          bytes = compressWeight<W>(tmp, bytes, w);
        } else {
          bytes = compressEdge(tmp, bytes, v - last_ngh);
          bytes = compressWeight<W>(tmp, bytes, w);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
        return false;
      };
      GA.get_vertex(i).mapOutNgh(i, f, false);

      degrees[i] = deg;
      byte_offsets[i] = total_bytes;
    }, 1);
    byte_offsets[n] = 0;
    size_t total_space = pbbslib::scan_add_inplace(byte_offsets.slice());
    cout << "total_space = " << total_space << endl;

    // 2. Create compressed format in-memory
    auto edges = pbbs::sequence<uchar>(total_space);
    parallel_for(0, n, [&] (size_t i) {
      uintE deg = degrees[i];
      assert(deg == GA.get_vertex(i).getOutDegree());
      if (deg > 0) {
        auto it = GA.get_vertex(i).getOutIter(i);
        long nbytes = byte::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it);
        if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
          cout << "nbytes = " << nbytes << " but offs = " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << endl;
          exit(0);
        }
        assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
      }
    }, 1);

    long* sizes = pbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = total_space;
    out.write((char*)sizes,sizeof(long)*3); //write n, m and space used
    out.write((char*)byte_offsets.begin(),sizeof(uintT)*(n+1)); //write offsets
    out.write((char*)degrees.begin(),sizeof(uintE)*n);
    out.write((char*)edges.begin(),total_space); //write edges
    out.close();
  }
} // namespace byte

namespace bytepd_amortized {

  template <template <class W> class vertex, class W>
  void write_graph_bytepd_amortized_directed(symmetric_graph<vertex, W>& GA, ofstream& out) {
    size_t n = GA.n; size_t m = GA.m;

    // out-edges
    // 1. Calculate total size
    {
      auto degrees = pbbs::sequence<uintE>(n);
      auto byte_offsets = pbbs::sequence<uintT>(n+1);
      parallel_for(0, n, [&] (size_t i) {
        size_t total_bytes = 0;
        uintE last_ngh = 0;
        size_t deg = 0;
        uchar tmp[16];
        auto f = [&] (uintE u, uintE v, W w) {
          long bytes = 0;
          if ((deg % PAR_DEGREE_TWO) == 0) {
            bytes = compressFirstEdge(tmp, bytes, u, v);
            bytes = compressWeight<W>(tmp, bytes, w);
          } else {
            bytes = compressEdge(tmp, bytes, v - last_ngh);
            bytes = compressWeight<W>(tmp, bytes, w);
          }
          last_ngh = v;
          total_bytes += bytes;
          deg++;
          return false;
        };
        GA.get_vertex(i).mapOutNgh(i, f, false);

        if (deg > 0) {
          size_t n_chunks = 1+(deg-1)/PAR_DEGREE_TWO;
          // To account for the byte offsets
          total_bytes += (n_chunks-1)*sizeof(uintE);
          // To account for the per-block counters
          total_bytes += (n_chunks)*sizeof(uintE);
          // To account for the virtual degree
          total_bytes += sizeof(uintE);
        }

        degrees[i] = deg;
        byte_offsets[i] = total_bytes;
      }, 1);
      byte_offsets[n] = 0;
      size_t total_space = pbbslib::scan_add_inplace(byte_offsets.slice());
      cout << "total in-space = " << total_space << endl;

      // 2. Create compressed format in-memory
      auto edges = pbbs::sequence<uchar>(total_space);
      parallel_for(0, n, [&] (size_t i) {
        uintE deg = degrees[i];
        if (deg > 0) {
          auto it = GA.get_vertex(i).getOutIter(i);
          long nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
          if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
            cout << "nbytes = " << nbytes << ". Should be: " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << endl;
            exit(0);
          }
          assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
        }
      }, 1);
      cout << "Compressed" << endl;

      long* sizes = pbbs::new_array_no_init<long>(3);
      sizes[0] = GA.n;
      sizes[1] = GA.m;
      sizes[2] = total_space;
      out.write((char*)sizes,sizeof(long)*3); //write n, m and space used
      out.write((char*)byte_offsets.begin(),sizeof(uintT)*(n+1)); //write offsets
      out.write((char*)degrees.begin(),sizeof(uintE)*n);
      out.write((char*)edges.begin(),total_space); //write edges
    }

    {
      // in-edges
      // 1. Calculate total size
      auto degrees = pbbs::sequence<uintE>(n);
      auto byte_offsets = pbbs::sequence<uintT>(n+1);
      parallel_for(0, n, [&] (size_t i) {
        size_t total_bytes = 0;
        uintE last_ngh = 0;
        size_t deg = 0;
        uchar tmp[16];
        auto f = [&] (uintE u, uintE v, W w) {
          long bytes = 0;
          if ((deg % PAR_DEGREE_TWO) == 0) {
            bytes = compressFirstEdge(tmp, bytes, u, v);
            bytes = compressWeight<W>(tmp, bytes, w);
          } else {
            bytes = compressEdge(tmp, bytes, v - last_ngh);
            bytes = compressWeight<W>(tmp, bytes, w);
          }
          last_ngh = v;
          total_bytes += bytes;
          deg++;
          return false;
        };
        GA.get_vertex(i).mapInNgh(i, f, false);

        if (deg > 0) {
          size_t n_chunks = 1+(deg-1)/PAR_DEGREE_TWO;
          // To account for the byte offsets
          total_bytes += (n_chunks-1)*sizeof(uintE);
          // To account for the per-block counters
          total_bytes += (n_chunks)*sizeof(uintE);
          // To account for the virtual degree
          total_bytes += sizeof(uintE);
        }

        degrees[i] = deg;
        byte_offsets[i] = total_bytes;
      }, 1);
      byte_offsets[n] = 0;
      size_t total_space = pbbslib::scan_add_inplace(byte_offsets.slice());
      cout << "total in-space = " << total_space << endl;

      // 2. Create compressed format in-memory
      auto edges = pbbs::sequence<uchar>(total_space);
      parallel_for (0, n, [&] (size_t i) {
        uintE deg = degrees[i];
        if (deg > 0) {
          auto it = GA.get_vertex(i).getInIter(i);
          long nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
          if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
            cout << "nbytes = " << nbytes << ". Should be: " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << endl;
            exit(0);
          }
          assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
        }
      }, 1);
      cout << "Compressed" << endl;
      long inTotalSpace[1];
      inTotalSpace[0] = total_space;
      out.write((char*)inTotalSpace, sizeof(long)); // in-edges total space
      out.write((char*)byte_offsets.begin(),sizeof(uintT)*(n+1)); //write offsets
      out.write((char*)degrees.begin(),sizeof(uintE)*n);
      out.write((char*)edges.begin(),total_space); //write edges
    }
  }

  template <template <class W> class vertex, class W>
  void write_graph_bytepd_amortized_format(symmetric_graph<vertex, W>& GA, ofstream& out, bool symmetric) {
    if (!symmetric) {
      write_graph_bytepd_amortized_directed(GA, out);
      return;
    }
    size_t n = GA.n; size_t m = GA.m;

//    auto xors = pbbs::sequence<size_t>(n);
//    parallel_for(size_t i=0; i<n; i++) {
//      size_t xr = 0;
//      auto map_f = wrap_f<W>([&] (uintE src, uintE ngh) {
//        xr ^= (src ^ ngh);
//      });
//      GA.V[i].mapOutNgh(i, map_f, false);
//      xors[i] = xr;
//    }
//    cout << "input graph: output red = " << pbbs::reduce_xor(xors) << endl;
//
//  auto hash_or_lt = [&] (const uintE& src, const uintE& ngh) {
//    uint32_t src_h = pbbs::hash32(src);
//    uint32_t ngh_h = pbbs::hash32(ngh);
//    return (src_h < ngh_h) || ((src_h == ngh_h) && src < ngh);
//  };


//    auto self_arr = pbbs::sequence<size_t>(n);
//    parallel_for(size_t i=0; i<n; i++) {
//      uintE our_deg = pbbs::log2_up(GA.V[i].getOutDegree());
//      bool selfl = false;
//      size_t pri = 0;
//      auto map_f = wrap_f<W>([&] (uintE src, uintE ngh) {
//        uintE ngh_deg = pbbs::log2_up(GA.V[ngh].getOutDegree());
//        if (src == ngh) {
//          selfl = true;
//        }
//        if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) && hash_or_lt(src, ngh))) {
//          pri++;
//        }
//      });
//      GA.V[i].mapOutNgh(i, map_f, false);
//      self_arr[i] = selfl;
//      xors[i] = pri;
//    }
//    cout << "input graph: priorities = " << pbbslib::reduce_add(xors) << endl;
//    cout << "input graph: self-loops = " << pbbslib::reduce_add(self_arr) << endl;

//    parallel_for(size_t i=0; i<n; i++) {
//      uintE our_deg = pbbs::log2_up(GA.V[i].getOutDegree());
//      bool selfl = false;
//      size_t pri = 0;
//      auto it = GA.V[i].getOutIter(i);
//      size_t d = 0;
//      size_t degree = GA.V[i].getOutDegree();
//      uintE src = i;
//      if (degree > 0) {
//        uintE ngh = get<0>(it.cur());
//        if (src == ngh) { selfl = true; }
//        uintE ngh_deg = pbbs::log2_up(GA.V[ngh].getOutDegree());
//        if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) && hash_or_lt(src, ngh))) {
//          pri++;
//        }
//        for (size_t i=1; i<degree; i++) {
//          ngh = get<0>(it.next());
//          ngh_deg = pbbs::log2_up(GA.V[ngh].getOutDegree());
//          if (src == ngh) { selfl = true; }
//          if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) && hash_or_lt(src, ngh))) {
//            pri++;
//          }
//        }
//      }
//      self_arr[i] = selfl;
//      xors[i] = pri;
//    }
//    cout << "input graph: priorities = " << pbbslib::reduce_add(xors) << endl;
//    cout << "input graph: self-loops = " << pbbslib::reduce_add(self_arr) << endl;


    // 1. Calculate total size
    auto degrees = pbbs::sequence<uintE>(n);
    auto byte_offsets = pbbs::sequence<uintT>(n+1);
    parallel_for(0, n, [&] (size_t i) {
      size_t total_bytes = 0;
      uintE last_ngh = 0;
      size_t deg = 0;
      uchar tmp[16];
      auto f = [&] (uintE u, uintE v, W w) {
//        if (u == v) {
//          return;
//        }
        long bytes = 0;
        if ((deg % PAR_DEGREE_TWO) == 0) {
          bytes = compressFirstEdge(tmp, bytes, u, v);
          bytes = compressWeight<W>(tmp, bytes, w);
        } else {
          bytes = compressEdge(tmp, bytes, v - last_ngh);
          bytes = compressWeight<W>(tmp, bytes, w);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
      };
      GA.get_vertex(i).mapOutNgh(i, f, false);

      if (deg > 0) {
        size_t n_chunks = 1+(deg-1)/PAR_DEGREE_TWO;
        // To account for the byte offsets
        total_bytes += (n_chunks-1)*sizeof(uintE);
        // To account for the per-block counters
        total_bytes += (n_chunks)*sizeof(uintE);
        // To account for the virtual degree
        total_bytes += sizeof(uintE);
      }

      degrees[i] = deg;
      byte_offsets[i] = total_bytes;
    }, 1);
    byte_offsets[n] = 0;
    size_t total_space = pbbslib::scan_add_inplace(byte_offsets.slice());
    cout << "total space = " << total_space << endl;
    auto deg_im = pbbs::delayed_seq<size_t>(n, [&] (size_t i) { return degrees[i]; });
    cout << "sum degs = " << pbbslib::reduce_add(deg_im) << endl;

    // 2. Create compressed format in-memory
    auto edges = pbbs::sequence<uchar>(total_space);
    parallel_for (0, n, [&] (size_t i) {
      uintE deg = degrees[i];
      if (deg > 0) {
        auto it = GA.get_vertex(i).getOutIter(i);
        long nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);

//        uchar* edgeArray = edges.begin() + byte_offsets[i];
//        size_t degree = deg;
//
//        size_t current_offset = 0;
//        size_t num_blocks = 1+(degree-1)/PAR_DEGREE_TWO;
//        uintE* vertex_ctr = (uintE*)edgeArray;
//        *vertex_ctr = degree;
//
//        uintE* block_offsets = (uintE*)(edgeArray + sizeof(uintE));
//        current_offset += sizeof(uintE) + (num_blocks-1)*sizeof(uintE); // virtual deg + block_offs
//
//        size_t proc = 0;
//        size_t cur_block = 0;
//        uintE* prev_block_deg = nullptr;
//        uintE last_ngh;
//        auto map_f = [&] (uintE u, uintE v, W w) {
//          if ((proc % PAR_DEGREE_TWO) == 0) {
//            if (cur_block > 0) {
//              assert(*prev_block_deg == 0);
//              *prev_block_deg = PAR_DEGREE_TWO; // full block; write prev block's degree
//              block_offsets[cur_block-1] = current_offset; // write start of this block
//            }
//            prev_block_deg = (uintE*)(edgeArray + current_offset);
//            *prev_block_deg = 0;
//            current_offset += sizeof(uintE);
//            current_offset = compressFirstEdge(edgeArray, current_offset, u, v);
//            current_offset = compressWeight<W>(edgeArray, current_offset, w);
//            cur_block++;
//          } else {
//            current_offset = compressEdge(edgeArray, current_offset, v - last_ngh);
//            current_offset = compressWeight<W>(edgeArray, current_offset, w);
//          }
//          last_ngh = v;
//          proc++;
//        };
//
//        GA.V[i].mapOutNgh(i, map_f, false);
//
//        assert(prev_block_deg != nullptr);
//        assert(*prev_block_deg == 0);
//        *prev_block_deg = deg % PAR_DEGREE_TWO;
//
//        assert(current_offset == (byte_offsets[i+1] - byte_offsets[i]));

        if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
          cout << "nbytes = " << nbytes << ". Should be: " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << endl;
          exit(0);
        }
        assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
      }
    }, 1);
    cout << "Compressed" << endl;
//    exit(0);

//    parallel_for(size_t i=0; i<n; i++) {
//      size_t xr = 0;
//      auto map_f = [&] (uintE src, uintE ngh, const W& wgh, size_t off) {
//        xr ^= (src ^ ngh);
//        return true;
//      };
//      auto edge_start= edges.begin() + byte_offsets[i];
//      size_t deg = degrees[i];
//      if (deg > 0) {
//        bytepd_amortized::decode<W>(map_f, edge_start, i, deg, false);
//      }
//      xors[i] = xr;
//    }
//    cout << "output graph: output red = " << pbbs::reduce_xor(xors) << endl;

//    parallel_for(size_t i=0; i<n; i++) {
//      assert(degrees[i] == GA.V[i].getOutDegree());
//      uintE our_deg = pbbs::log2_up(degrees[i]);
//      bool selfl = false;
//      size_t pri = 0;
//      auto map_f = [&] (uintE src, uintE ngh, const W& wgh, size_t off) {
//        uintE ngh_deg = pbbs::log2_up(degrees[ngh]);
//        if (src == ngh) {
//          selfl = true;
//        }
//        if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) && hash_or_lt(src, ngh))) {
//          pri++;
//        }
//        return true;
//      };
//      auto edge_start= edges.begin() + byte_offsets[i];
//      size_t deg = degrees[i];
//      if (deg > 0) {
//        bytepd_amortized::decode<W>(map_f, edge_start, i, deg, false);
//      }
//      self_arr[i] = selfl;
//      xors[i] = pri;
//    }
//    cout << "output graph: priorities = " << pbbslib::reduce_add(xors) << endl;
//    cout << "output graph: self-loops = " << pbbslib::reduce_add(self_arr) << endl;


//    exit(0);

    long* sizes = pbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = total_space;
    out.write((char*)sizes,sizeof(long)*3); //write n, m and space used
    out.write((char*)byte_offsets.begin(),sizeof(uintT)*(n+1)); //write offsets
    out.write((char*)degrees.begin(),sizeof(uintE)*n);
    out.write((char*)edges.begin(),total_space); //write edges
    out.close();
  }
}; // namespace bytepd_amortized

namespace binary_format {

  template <template <class W> class vertex, class W>
  void write_graph_binary_format(symmetric_graph<vertex, W>& GA, ofstream& out, size_t n_batches=4) {
    size_t n = GA.n; size_t m = GA.m;
    using edge_type = std::tuple<uintE, W>;

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
        auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
          return wgh;
        };
        auto write_f = [&] (const uintE& ngh, const uintT& offset, const W& val) {
          edges[offset] = std::make_tuple(ngh, val);
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
auto converter(G& GA, commandLine P) {
  auto outfile = P.getOptionValue("-o", "");
  bool symmetric = P.getOptionValue("-s");
  cout << "outfile = " << outfile << endl;
  if (outfile == "") {
    cout << "specify a valid outfile" << endl;
    exit(0);
  }
  ofstream out(outfile.c_str(), ofstream::out | ios::binary);
  auto encoding = P.getOptionValue("-enc", "byte");

  if (encoding == "byte") {
    byte::write_graph_byte_format(GA, out, symmetric);
  } else if (encoding == "bytepd-amortized") {
    bytepd_amortized::write_graph_bytepd_amortized_format(GA, out, symmetric);
  } else if (encoding == "binary") {
    binary_format::write_graph_binary_format(GA, out);
  } else {
    cout << "Unknown encoding: " << encoding << endl;
    exit(0);
  }
  return static_cast<double>(0);
}


generate_main(converter, false);
