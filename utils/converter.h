#pragma once

#include "ligra/ligra.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

#include "to_char_arr.h"

// Provides utilities for converting between different compressed
// representations.
constexpr int max_weight = 32;
size_t PAR_DEGREE_TWO = 1000;

namespace byte {
  template <class Graph>
  void write_graph_byte_format(Graph& GA, std::ofstream& out, bool symmetric) {
    using W = typename Graph::weight_type;
    size_t n = GA.n;

    // 1. Calculate total size
    auto degrees = pbbs::sequence<uintE>(n);
    auto byte_offsets = pbbs::sequence<uintT>(n+1);
    std::cout << "# calculating size" << std::endl;
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
    std::cout << "# total_space = " << total_space << std::endl;

    // 2. Create compressed format in-memory
    auto edges = pbbs::sequence<uchar>(total_space);
    parallel_for(0, n, [&] (size_t i) {
      uintE deg = degrees[i];
      assert(deg == GA.get_vertex(i).getOutDegree());
      if (deg > 0) {
        auto it = GA.get_vertex(i).getOutIter(i);
        size_t nbytes = byte::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it);
        if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
        std::cout << "# nbytes = " << nbytes << " but offs = " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << std::endl;
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

  template <class Graph>
  void write_graph_bytepd_amortized_directed(Graph& GA, std::ofstream& out) {
    size_t n = GA.n;
    using W = typename Graph::weight_type;

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
      std::cout << "# total in-space = " << total_space << std::endl;

      // 2. Create compressed format in-memory
      auto edges = pbbs::sequence<uchar>(total_space);
      parallel_for(0, n, [&] (size_t i) {
        uintE deg = degrees[i];
        if (deg > 0) {
          auto it = GA.get_vertex(i).getOutIter(i);
          size_t nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
          if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
          std::cout << "# nbytes = " << nbytes << ". Should be: " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << std::endl;
            exit(0);
          }
          assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
        }
      }, 1);
      std::cout << "# Compressed" << std::endl;

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
      std::cout << "# total in-space = " << total_space << std::endl;

      // 2. Create compressed format in-memory
      auto edges = pbbs::sequence<uchar>(total_space);
      parallel_for (0, n, [&] (size_t i) {
        uintE deg = degrees[i];
        if (deg > 0) {
          auto it = GA.get_vertex(i).getInIter(i);
          size_t nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
          if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
          std::cout << "# nbytes = " << nbytes << ". Should be: " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << std::endl;
            exit(0);
          }
          assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
        }
      }, 1);
      std::cout << "# Compressed" << std::endl;
      long inTotalSpace[1];
      inTotalSpace[0] = total_space;
      out.write((char*)inTotalSpace, sizeof(long)); // in-edges total space
      out.write((char*)byte_offsets.begin(),sizeof(uintT)*(n+1)); //write offsets
      out.write((char*)degrees.begin(),sizeof(uintE)*n);
      out.write((char*)edges.begin(),total_space); //write edges
    }
  }

  template <class Graph>
  void write_graph_bytepd_amortized_format(Graph& GA, std::ofstream& out, bool symmetric, size_t n_batches = 6) {
    if (!symmetric) {
      write_graph_bytepd_amortized_directed(GA, out);
      return;
    }
    size_t n = GA.n;
    using W = typename Graph::weight_type;

//    auto xors = pbbs::sequence<size_t>(n);
//    parallel_for(size_t i=0; i<n; i++) {
//      size_t xr = 0;
//      auto map_f = wrap_f<W>([&] (uintE src, uintE ngh) {
//        xr ^= (src ^ ngh);
//      });
//      GA.V[i].mapOutNgh(i, map_f, false);
//      xors[i] = xr;
//    }
//    cout << "input graph: output red = " << pbbs::reduce_xor(xors) << std::endl;
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
//    cout << "input graph: priorities = " << pbbslib::reduce_add(xors) << std::endl;
//    cout << "input graph: self-loops = " << pbbslib::reduce_add(self_arr) << std::endl;

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
//    cout << "input graph: priorities = " << pbbslib::reduce_add(xors) << std::endl;
//    cout << "input graph: self-loops = " << pbbslib::reduce_add(self_arr) << std::endl;


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
    std::cout << "# total space = " << total_space << std::endl;
    auto deg_im = pbbs::delayed_seq<size_t>(n, [&] (size_t i) { return degrees[i]; });
    std::cout << "# sum degs = " << pbbslib::reduce_add(deg_im) << std::endl;

    long* sizes = pbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = total_space;
    out.write((char*)sizes,sizeof(long)*3); //write n, m and space used
    out.write((char*)byte_offsets.begin(),sizeof(uintT)*(n+1)); //write offsets
    out.write((char*)degrees.begin(),sizeof(uintE)*n);



    // 2. Create compressed format in-memory

    size_t bs = pbbs::num_blocks(n, n_batches);

    for (size_t i=0; i<bs; i++) {
      size_t start = i*bs;
      size_t end = std::min(start+bs, n);
      if (start >= end) break;
      std::cout << "# writing vertices " << start << " to " << end << std::endl;

      // create slab of graph and write out
      size_t start_offset = byte_offsets[start];
      size_t end_offset = byte_offsets[end];
      size_t n_bytes = end_offset - start_offset;
      uchar* edges = pbbs::new_array_no_init<uchar>(n_bytes);
      parallel_for(start, end, [&] (size_t j) {
        size_t our_offset = byte_offsets[j] - start_offset;
        uintE deg = degrees[j];
        if (deg > 0) {
          auto it = GA.get_vertex(j).getOutIter(j);
          size_t nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges + our_offset, 0, deg, (uintE)j, it, PAR_DEGREE_TWO);

          if (nbytes != (byte_offsets[j+1] - byte_offsets[j])) {
          std::cout << "# nbytes = " << nbytes << ". Should be: " << (byte_offsets[j+1] - byte_offsets[j]) << " deg = " << deg << " j = " << j << std::endl;
            exit(0);
          }
          assert(nbytes == (byte_offsets[j+1] - byte_offsets[j]));
        }
      });
      out.write((char*)edges,n_bytes); //write edges
      std::cout << "# finished writing vertices " << start << " to " << end << std::endl;
      pbbs::free_array(edges);
    }
////    parallel_for (0, n, [&] (size_t i) {
//      uintE deg = degrees[i];
//      if (deg > 0) {
//        auto it = GA.get_vertex(i).getOutIter(i);
//        size_t nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
//
//        if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
//          cout << "nbytes = " << nbytes << ". Should be: " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << std::endl;
//          exit(0);
//        }
//        assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
//      }
//    }, 1);
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
//    cout << "output graph: output red = " << pbbs::reduce_xor(xors) << std::endl;

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
//    cout << "output graph: priorities = " << pbbslib::reduce_add(xors) << std::endl;
//    cout << "output graph: self-loops = " << pbbslib::reduce_add(self_arr) << std::endl;


//    exit(0);

    out.close();
  }

  template <class Graph>
  inline uintE* rankNodes(Graph& GA, size_t n) {
    uintE* r = pbbslib::new_array_no_init<uintE>(n);
  //  uintE* o = pbbslib::new_array_no_init<uintE>(n);
    sequence<uintE> o(n);

    timer t;
    t.start();
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) { o[i] = i; });
    pbbslib::sample_sort_inplace(o.slice(), [&](const uintE u, const uintE v) {
      return GA.get_vertex(u).getOutDegree() < GA.get_vertex(v).getOutDegree();
    });
    par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                    { r[o[i]] = i; });
    t.stop();
    debug(t.reportTotal("Rank time"););
    return r;
  }

  template <class Graph>
  void degree_reorder(Graph& GA, std::ofstream& out, bool symmetric, size_t n_batches = 6) {
    if (!symmetric) {
      assert(false);
      exit(0);
    }
    size_t n = GA.n;
    using W = typename Graph::weight_type;

    // mapping from v -> new id
    uintE* rank = rankNodes(GA, GA.n);

    auto inverse_rank = pbbs::sequence<uintE>(n);
    parallel_for(0, n, [&] (size_t i) {
      uintE rank_i = rank[i];
      inverse_rank[rank_i] = i;
    });

    // 1. Calculate total size
    auto degrees = pbbs::sequence<uintE>(n);
    auto byte_offsets = pbbs::sequence<uintT>(n+1);
    parallel_for(0, n, [&] (size_t i) {
      size_t total_bytes = 0;
      uintE last_ngh = 0;
      size_t deg = 0;
      uchar tmp[16];

      uintE stk[8192];
      uintE* nghs = (uintE*)stk;
      auto vtx = GA.get_vertex(i);
      if (vtx.getOutDegree() > 0) {
        if (vtx.getOutDegree() > 8192) {
          nghs = pbbs::new_array_no_init<uintE>(deg);
        }

        size_t k = 0;
        auto map_ngh_f = [&] (const uintE& u, const uintE& w, const W& wgh) {
          nghs[k++] = rank[w];
        };
        vtx.mapOutNgh(i, map_ngh_f, false);

        auto new_ngh_seq = pbbslib::make_sequence(nghs, deg);
        pbbs::sample_sort_inplace(new_ngh_seq, std::less<uintE>());

        uintE our_new_id = rank[i];

        for (size_t j=0; j<new_ngh_seq.size(); j++) {
          long bytes = 0;
          uintE ngh_id = new_ngh_seq[j];
          if ((deg % PAR_DEGREE_TWO) == 0) {
            bytes = compressFirstEdge(tmp, bytes, our_new_id, ngh_id);
          } else {
            bytes = compressEdge(tmp, bytes, ngh_id - last_ngh);
          }
          last_ngh = ngh_id;
          total_bytes += bytes;
          deg++;
        }

        if (deg > 0) {
          size_t n_chunks = 1+(deg-1)/PAR_DEGREE_TWO;
          // To account for the byte offsets
          total_bytes += (n_chunks-1)*sizeof(uintE);
          // To account for the per-block counters
          total_bytes += (n_chunks)*sizeof(uintE);
          // To account for the virtual degree
          total_bytes += sizeof(uintE);
        }
        if (vtx.getOutDegree() > 8192) {
          pbbs::free_array(nghs);
        }
      }
      degrees[i] = deg;
      byte_offsets[i] = total_bytes;
    }, 1);
    byte_offsets[n] = 0;
    size_t total_space = pbbslib::scan_add_inplace(byte_offsets.slice());
    std::cout << "# total space = " << total_space << std::endl;
    auto deg_im = pbbs::delayed_seq<size_t>(n, [&] (size_t i) { return degrees[i]; });
    std::cout << "# sum degs = " << pbbslib::reduce_add(deg_im) << std::endl;

    long* sizes = pbbs::new_array_no_init<long>(3);
    sizes[0] = GA.n;
    sizes[1] = GA.m;
    sizes[2] = total_space;
    out.write((char*)sizes,sizeof(long)*3); //write n, m and space used
    out.write((char*)byte_offsets.begin(),sizeof(uintT)*(n+1)); //write offsets
    out.write((char*)degrees.begin(),sizeof(uintE)*n);

    // 2. Create compressed format in-memory

    size_t bs = pbbs::num_blocks(n, n_batches);

    for (size_t i=0; i<bs; i++) {
      size_t start = i*bs;
      size_t end = std::min(start+bs, n);
      if (start >= end) break;
      std::cout << "# writing vertices " << start << " to " << end << std::endl;

      // create slab of graph and write out
      size_t start_offset = byte_offsets[start];
      size_t end_offset = byte_offsets[end];
      size_t n_bytes = end_offset - start_offset;
      uchar* edges = pbbs::new_array_no_init<uchar>(n_bytes);
      parallel_for(start, end, [&] (size_t j) {
        size_t our_offset = byte_offsets[j] - start_offset;
        uintE deg = degrees[j];
        if (deg > 0) {
          auto it = GA.get_vertex(j).getOutIter(j);
          size_t nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges + our_offset, 0, deg, (uintE)j, it, PAR_DEGREE_TWO);

          if (nbytes != (byte_offsets[j+1] - byte_offsets[j])) {
          std::cout << "# nbytes = " << nbytes << ". Should be: " << (byte_offsets[j+1] - byte_offsets[j]) << " deg = " << deg << " j = " << j << std::endl;
            exit(0);
          }
          assert(nbytes == (byte_offsets[j+1] - byte_offsets[j]));
        }
      });
      out.write((char*)edges,n_bytes); //write edges
      std::cout << "# finished writing vertices " << start << " to " << end << std::endl;
      pbbs::free_array(edges);
    }
////    parallel_for (0, n, [&] (size_t i) {
//      uintE deg = degrees[i];
//      if (deg > 0) {
//        auto it = GA.get_vertex(i).getOutIter(i);
//        long nbytes = bytepd_amortized::sequentialCompressEdgeSet<W>(edges.begin() + byte_offsets[i], 0, deg, (uintE)i, it, PAR_DEGREE_TWO);
//
//        if (nbytes != (byte_offsets[i+1] - byte_offsets[i])) {
//          cout << "nbytes = " << nbytes << ". Should be: " << (byte_offsets[i+1] - byte_offsets[i]) << " deg = " << deg << " i = " << i << std::endl;
//          exit(0);
//        }
//        assert(nbytes == (byte_offsets[i+1] - byte_offsets[i]));
//      }
//    }, 1);
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
//    cout << "output graph: output red = " << pbbs::reduce_xor(xors) << std::endl;

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
//    cout << "output graph: priorities = " << pbbslib::reduce_add(xors) << std::endl;
//    cout << "output graph: self-loops = " << pbbslib::reduce_add(self_arr) << std::endl;


//    exit(0);

    out.close();
  }


}; // namespace bytepd_amortized

namespace binary_format {

  template <class Graph>
  void write_graph_binary_format(Graph& GA, std::ofstream& out, size_t n_batches=4) {
    size_t n = GA.n; size_t m = GA.m;
    using W = typename Graph::weight_type;
    using edge_type = std::tuple<uintE, W>;

    // 1. Calculate total size
    auto offsets = pbbs::sequence<uintT>(n+1);
    std::cout << "# calculating size" << std::endl;
    parallel_for(0, n, [&] (size_t i) {
      offsets[i] = GA.get_vertex(i).getOutDegree();
    });
    offsets[n] = 0;
    size_t offset_scan = pbbslib::scan_add_inplace(offsets.slice());
    std::cout << "# offset_scan = " << offset_scan << " m = " << m << std::endl;
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
      std::cout << "# writing vertices " << start << " to " << end << std::endl;

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

      std::cout << "# finished writing vertices " << start << " to " << end << std::endl;
      edges_written += n_edges;

      pbbs::free_array(edges);
    }
    std::cout << "# Wrote " << edges_written << " edges in total, m = " << m << std::endl;
    assert(edges_written == m);
    out.close();
  }


}; // namespace binary_format

template <class Graph>
void edgearray(Graph& GA, std::ofstream& out) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto degs = pbbs::sequence<uintT>(n);
  par_for(0, n, [&] (size_t i) { degs[i] = GA.get_vertex(i).getOutDegree(); });
  pbbs::scan_inplace(degs.slice(), pbbs::addm<uintT>());

  auto edges = pbbs::sequence<std::tuple<uintE, uintE, W>>(m);

  parallel_for(0, n, [&] (size_t i) {
    size_t off = degs[i];
    size_t k = 0;
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      edges[off + k++] = std::make_tuple(u, v, wgh);
    };
    GA.get_vertex(i).mapOutNgh(i, map_f, false);
    assert(k == GA.get_vertex(i).getOutDegree());
  }, 1);

  writeArrayToStream(out, edges); // m edge triples
  out.close();
  std::cout << "# Wrote file." << std::endl;
}

template <class Graph>
auto converter(Graph& GA, commandLine P) {
  auto outfile = P.getOptionValue("-o", "");
  bool symmetric = P.getOptionValue("-s");
  std::cout << "# outfile = " << outfile << std::endl;
  if (outfile == "") {
    std::cout << "# specify a valid outfile" << std::endl;
    exit(0);
  }
  std::ofstream out(outfile.c_str(), std::ofstream::out | std::ios::binary);
  auto encoding = P.getOptionValue("-enc", "byte");

  PAR_DEGREE_TWO = P.getOptionLongValue("-bs", 256);

  if (encoding == "byte") {
    byte::write_graph_byte_format(GA, out, symmetric);
  } else if (encoding == "bytepd-amortized") {
    bytepd_amortized::write_graph_bytepd_amortized_format(GA, out, symmetric);
  } else if (encoding == "binary") {
    binary_format::write_graph_binary_format(GA, out);
  } else if (encoding == "degree") {
    bytepd_amortized::degree_reorder(GA, out, symmetric);
  } else if (encoding == "edgearray") {
    edgearray(GA, out);
  } else {
    std::cout << "# Unknown encoding: " << encoding << std::endl;
    exit(0);
  }
  return static_cast<double>(0);
}
