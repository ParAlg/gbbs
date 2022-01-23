#include "gbbs/encodings/byte_pd_amortized.h"
#include "gbbs/gbbs.h"
#include "gbbs/helpers/parse_command_line.h"
#include "gbbs/io.h"

#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>

// Provides utilities for converting between different compressed
// representations.

namespace gbbs {

template <class Graph>
void write_graph_bytepd_amortized_directed(Graph& GA, std::ofstream& out) {
  namespace encoding = bytepd_amortized;
  using W = typename Graph::weight_type;
  size_t n = GA.n;

  // out-edges
  // 1. Calculate total size
  {
    auto degrees = sequence<uintE>(n);
    auto byte_offsets = sequence<uintT>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      size_t total_bytes = 0;
      uintE last_ngh = 0;
      size_t deg = 0;
      uchar tmp[16];
      auto f = [&](uintE u, uintE v, W w) {
        long bytes = 0;
        if ((deg % PARALLEL_DEGREE) == 0) {
          bytes = encoding::compressFirstEdge(tmp, bytes, u, v);
          bytes = encoding::compressWeight<W>(tmp, bytes, w);
        } else {
          bytes = encoding::compressEdge(tmp, bytes, v - last_ngh);
          bytes = encoding::compressWeight<W>(tmp, bytes, w);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
        return false;
      };
      GA.get_vertex(i).out_neighbors().map(f, false);

      if (deg > 0) {
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
    });
    byte_offsets[n] = 0;
    size_t total_space = parlay::scan_inplace(make_slice(byte_offsets));
    std::cout << "total in-space = " << total_space << std::endl;

    // 2. Create compressed format in-memory
    auto edges = sequence<uchar>(total_space);
    parallel_for(0, n, [&](size_t i) {
      uintE deg = degrees[i];
      if (deg > 0) {
        auto it = GA.get_vertex(i).out_neighbors().get_iter();
        size_t nbytes = encoding::sequentialCompressEdgeSet<W>(
            edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it);
        if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
          std::cout << "nbytes = " << nbytes << ". Should be: "
                    << (byte_offsets[i + 1] - byte_offsets[i])
                    << " deg = " << deg << " i = " << i << std::endl;
          exit(0);
        }
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    });
    std::cout << "Compressed" << std::endl;

    long* sizes = gbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = total_space;
    out.write((char*)sizes, sizeof(long) * 3);  // write n, m and space used
    out.write((char*)byte_offsets.begin(),
              sizeof(uintT) * (n + 1));  // write offsets
    out.write((char*)degrees.begin(), sizeof(uintE) * n);
    out.write((char*)edges.begin(), total_space);  // write edges
  }

  {
    // in-edges
    // 1. Calculate total size
    auto degrees = sequence<uintE>(n);
    auto byte_offsets = sequence<uintT>(n + 1);
    parallel_for(0, n, [&](size_t i) {
      size_t total_bytes = 0;
      uintE last_ngh = 0;
      size_t deg = 0;
      uchar tmp[16];
      auto f = [&](uintE u, uintE v, W w) {
        long bytes = 0;
        if ((deg % PARALLEL_DEGREE) == 0) {
          bytes = encoding::compressFirstEdge(tmp, bytes, u, v);
          bytes = encoding::compressWeight<W>(tmp, bytes, w);
        } else {
          bytes = encoding::compressEdge(tmp, bytes, v - last_ngh);
          bytes = encoding::compressWeight<W>(tmp, bytes, w);
        }
        last_ngh = v;
        total_bytes += bytes;
        deg++;
        return false;
      };
      GA.get_vertex(i).in_neighbors().map(f, false);

      if (deg > 0) {
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
    });
    byte_offsets[n] = 0;
    size_t total_space = parlay::scan_inplace(make_slice(byte_offsets));
    std::cout << "total in-space = " << total_space << std::endl;

    // 2. Create compressed format in-memory
    auto edges = sequence<uchar>(total_space);
    parallel_for(0, n, [&](size_t i) {
      uintE deg = degrees[i];
      if (deg > 0) {
        auto it = GA.get_vertex(i).in_neighbors().get_iter();
        size_t nbytes = encoding::sequentialCompressEdgeSet<W>(
            edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it);
        if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
          std::cout << "nbytes = " << nbytes << ". Should be: "
                    << (byte_offsets[i + 1] - byte_offsets[i])
                    << " deg = " << deg << " i = " << i << std::endl;
          exit(0);
        }
        assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
      }
    });
    std::cout << "Compressed" << std::endl;
    long inTotalSpace[1];
    inTotalSpace[0] = total_space;
    out.write((char*)inTotalSpace, sizeof(long));  // in-edges total space
    out.write((char*)byte_offsets.begin(),
              sizeof(uintT) * (n + 1));  // write offsets
    out.write((char*)degrees.begin(), sizeof(uintE) * n);
    out.write((char*)edges.begin(), total_space);  // write edges
  }
}

template <class Graph>
void write_graph_bytepd_amortized_format(Graph& GA, std::ofstream& out,
                                         bool symmetric) {
  namespace encoding = bytepd_amortized;
  using W = typename Graph::weight_type;
  if (!symmetric) {
    write_graph_bytepd_amortized_directed(GA, out);
    return;
  }
  size_t n = GA.n;

  //    auto xors = sequence<size_t>(n);
  //    parallel_for(size_t i=0; i<n; i++) {
  //      size_t xr = 0;
  //      auto map_f = wrap_f<W>([&] (uintE src, uintE ngh) {
  //        xr ^= (src ^ ngh);
  //      });
  //      GA.V[i].mapOutNgh(i, map_f, false);
  //      xors[i] = xr;
  //    }
  //    std::cout << "input graph: output red = " << parlay::reduce_xor(xors) <<
  //    std::endl;
  //
  //  auto hash_or_lt = [&] (const uintE& src, const uintE& ngh) {
  //    uint32_t src_h = parlay::hash32(src);
  //    uint32_t ngh_h = parlay::hash32(ngh);
  //    return (src_h < ngh_h) || ((src_h == ngh_h) && src < ngh);
  //  };

  //    auto self_arr = sequence<size_t>(n);
  //    parallel_for(size_t i=0; i<n; i++) {
  //      uintE our_deg = parlay::log2_up(GA.V[i].getOutDegree());
  //      bool selfl = false;
  //      size_t pri = 0;
  //      auto map_f = wrap_f<W>([&] (uintE src, uintE ngh) {
  //        uintE ngh_deg = parlay::log2_up(GA.V[ngh].getOutDegree());
  //        if (src == ngh) {
  //          selfl = true;
  //        }
  //        if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) && hash_or_lt(src,
  //        ngh))) {
  //          pri++;
  //        }
  //      });
  //      GA.V[i].mapOutNgh(i, map_f, false);
  //      self_arr[i] = selfl;
  //      xors[i] = pri;
  //    }
  //    std::cout << "input graph: priorities = " << parlay::reduce(xors) <<
  //    std::endl;
  //    std::cout << "input graph: self-loops = " << parlay::reduce(self_arr) <<
  //    std::endl;

  //    parallel_for(size_t i=0; i<n; i++) {
  //      uintE our_deg = parlay::log2_up(GA.V[i].getOutDegree());
  //      bool selfl = false;
  //      size_t pri = 0;
  //      auto it = GA.V[i].getOutIter(i);
  //      size_t d = 0;
  //      size_t degree = GA.V[i].getOutDegree();
  //      uintE src = i;
  //      if (degree > 0) {
  //        uintE ngh = get<0>(it.cur());
  //        if (src == ngh) { selfl = true; }
  //        uintE ngh_deg = parlay::log2_up(GA.V[ngh].getOutDegree());
  //        if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) && hash_or_lt(src,
  //        ngh))) {
  //          pri++;
  //        }
  //        for (size_t i=1; i<degree; i++) {
  //          ngh = get<0>(it.next());
  //          ngh_deg = parlay::log2_up(GA.V[ngh].getOutDegree());
  //          if (src == ngh) { selfl = true; }
  //          if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) &&
  //          hash_or_lt(src, ngh))) {
  //            pri++;
  //          }
  //        }
  //      }
  //      self_arr[i] = selfl;
  //      xors[i] = pri;
  //    }
  //    std::cout << "input graph: priorities = " << parlay::reduce(xors) <<
  //    std::endl;
  //    std::cout << "input graph: self-loops = " << parlay::reduce(self_arr) <<
  //    std::endl;

  // 1. Calculate total size
  auto degrees = sequence<uintE>(n);
  auto byte_offsets = sequence<uintT>(n + 1);
  parallel_for(0, n, [&](size_t i) {
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    uchar tmp[16];
    auto f = [&](uintE u, uintE v, W w) {
      //        if (u == v) {
      //          return;
      //        }
      long bytes = 0;
      if ((deg % PARALLEL_DEGREE) == 0) {
        bytes = encoding::compressFirstEdge(tmp, bytes, u, v);
        bytes = encoding::compressWeight<W>(tmp, bytes, w);
      } else {
        bytes = encoding::compressEdge(tmp, bytes, v - last_ngh);
        bytes = encoding::compressWeight<W>(tmp, bytes, w);
      }
      last_ngh = v;
      total_bytes += bytes;
      deg++;
    };
    GA.get_vertex(i).out_neighbors().map(f, false);

    if (deg > 0) {
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
  });
  byte_offsets[n] = 0;
  size_t total_space = parlay::scan_inplace(make_slice(byte_offsets));
  std::cout << "total space = " << total_space << std::endl;
  auto deg_f = [&](size_t i) { return degrees[i]; };
  auto deg_im = parlay::delayed_seq<size_t>(n, deg_f);
  std::cout << "sum degs = " << parlay::reduce(deg_im) << std::endl;

  // 2. Create compressed format in-memory
  auto edges = sequence<uchar>(total_space);
  parallel_for(0, n, [&](size_t i) {
    uintE deg = degrees[i];
    if (deg > 0) {
      auto it = GA.get_vertex(i).out_neighbors().get_iter();
      size_t nbytes = encoding::sequentialCompressEdgeSet<W>(
          edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it);

      //        uchar* edgeArray = edges.begin() + byte_offsets[i];
      //        size_t degree = deg;
      //
      //        size_t current_offset = 0;
      //        size_t num_blocks = 1+(degree-1)/PARALLEL_DEGREE;
      //        uintE* vertex_ctr = (uintE*)edgeArray;
      //        *vertex_ctr = degree;
      //
      //        uintE* block_offsets = (uintE*)(edgeArray + sizeof(uintE));
      //        current_offset += sizeof(uintE) + (num_blocks-1)*sizeof(uintE);
      //        // virtual deg + block_offs
      //
      //        size_t proc = 0;
      //        size_t cur_block = 0;
      //        uintE* prev_block_deg = nullptr;
      //        uintE last_ngh;
      //        auto map_f = [&] (uintE u, uintE v, W w) {
      //          if ((proc % PARALLEL_DEGREE) == 0) {
      //            if (cur_block > 0) {
      //              assert(*prev_block_deg == 0);
      //              *prev_block_deg = PARALLEL_DEGREE; // full block; write
      //              prev block's degree
      //              block_offsets[cur_block-1] = current_offset; // write
      //              start of this block
      //            }
      //            prev_block_deg = (uintE*)(edgeArray + current_offset);
      //            *prev_block_deg = 0;
      //            current_offset += sizeof(uintE);
      //            current_offset = compressFirstEdge(edgeArray,
      //            current_offset, u, v);
      //            current_offset = compressWeight<W>(edgeArray,
      //            current_offset, w);
      //            cur_block++;
      //          } else {
      //            current_offset = compressEdge(edgeArray, current_offset, v -
      //            last_ngh);
      //            current_offset = compressWeight<W>(edgeArray,
      //            current_offset, w);
      //          }
      //          last_ngh = v;
      //          proc++;
      //        };
      //
      //        GA.V[i].mapOutNgh(i, map_f, false);
      //
      //        assert(prev_block_deg != nullptr);
      //        assert(*prev_block_deg == 0);
      //        *prev_block_deg = deg % PARALLEL_DEGREE;
      //
      //        assert(current_offset == (byte_offsets[i+1] - byte_offsets[i]));

      if (nbytes != (byte_offsets[i + 1] - byte_offsets[i])) {
        std::cout << "nbytes = " << nbytes
                  << ". Should be: " << (byte_offsets[i + 1] - byte_offsets[i])
                  << " deg = " << deg << " i = " << i << std::endl;
        exit(0);
      }
      assert(nbytes == (byte_offsets[i + 1] - byte_offsets[i]));
    }
  });
  std::cout << "Compressed" << std::endl;
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
  //    std::cout << "output graph: output red = " << parlay::reduce_xor(xors)
  //    << std::endl;

  //    parallel_for(size_t i=0; i<n; i++) {
  //      assert(degrees[i] == GA.V[i].getOutDegree());
  //      uintE our_deg = parlay::log2_up(degrees[i]);
  //      bool selfl = false;
  //      size_t pri = 0;
  //      auto map_f = [&] (uintE src, uintE ngh, const W& wgh, size_t off) {
  //        uintE ngh_deg = parlay::log2_up(degrees[ngh]);
  //        if (src == ngh) {
  //          selfl = true;
  //        }
  //        if ((ngh_deg > our_deg) || ((ngh_deg == our_deg) && hash_or_lt(src,
  //        ngh))) {
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
  //    std::cout << "output graph: priorities = " << parlay::reduce(xors) <<
  //    std::endl;
  //    std::cout << "output graph: self-loops = " << parlay::reduce(self_arr)
  //    << std::endl;

  //    exit(0);

  long* sizes = gbbs::new_array_no_init<long>(3);
  sizes[0] = GA.n;
  sizes[1] = GA.m;
  sizes[2] = total_space;
  out.write((char*)sizes, sizeof(long) * 3);  // write n, m and space used
  out.write((char*)byte_offsets.begin(),
            sizeof(uintT) * (n + 1));  // write offsets
  out.write((char*)degrees.begin(), sizeof(uintE) * n);
  out.write((char*)edges.begin(), total_space);  // write edges
  out.close();
}

template <class Graph>
double converter(Graph& GA, commandLine P) {
  auto outfile = P.getOptionValue("-o", "");
  bool symmetric = P.getOptionValue("-s");
  std::cout << "Outfile: " << outfile << std::endl;
  if (outfile == "") {
    std::cout << "Please specify an output file" << std::endl;
    exit(0);
  }
  std::ofstream out(outfile.c_str(), std::ofstream::out | std::ios::binary);
  auto encoding = P.getOptionValue("-enc", "bytepd-amortized");

  if (encoding == "bytepd-amortized") {
    write_graph_bytepd_amortized_format(GA, out, symmetric);
  } else {
    std::cout << "Unknown encoding: " << encoding << std::endl;
    exit(0);
  }
  std::cout << "Finished converting." << std::endl;
  exit(0);
  return 0;
}
}  // namespace gbbs

generate_main(gbbs::converter, false);
