// Based on AsyncST/CC.C
#pragma once

#include <algorithm>
#include "ligra.h"
#include "utils/stats.h"
#include "../benchmark/MIS.h"
#include "pbbslib/random_shuffle.h"
#include <unordered_map>

template <class G, class S>
inline symmetric_graph<symmetricVertex, pbbs::empty> renumber(G& GA, S& perm) {
  using W = typename G::weight_type;
  uintE n = GA.n;

  // build the new graph
  auto offsets = sequence<size_t>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE our_id = perm[i];
    auto count_f = [&](uintE src, uintE ngh, const W& wgh) {
      uintE ngh_id = perm[ngh];
      return ngh_id < our_id;
    };
    offsets[our_id] = GA.get_vertex(i).countOutNgh(i, count_f);
  }, 1);

  size_t new_m = pbbslib::scan_add_inplace(offsets.slice());
  cout << "new_m = " << new_m << endl;

  auto new_edges = pbbs::sequence<uintE>(new_m);
  parallel_for(0, n, [&] (size_t i) {
    uintE new_id = perm[i];
    size_t offset = offsets[new_id];
    size_t directed_degree = ((new_id == n-1) ? new_m : offsets[new_id+1]) - offset;

    size_t k = 0;
    auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
      uintE id_v = perm[v];
      if (id_v < new_id) {
        assert(offset + k <= new_m);
        new_edges[offset + k] = id_v;
        k++;
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f, false); // run seq
    assert(directed_degree == k);
    auto seq = pbbslib::make_sequence(new_edges.begin() + offset, k);
    // sort by id
    pbbs::sample_sort_inplace(seq, std::less<uintE>());
  });

  auto v = pbbs::sequence<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE new_id = perm[i];
    size_t o = offsets[new_id];
    size_t degree = ((new_id == n-1) ? new_m : offsets[new_id+1]) - o;
    v[new_id].degree = degree;
    v[new_id].offset = o;
  });

  auto vv = v.to_array();
  auto edge_arr = (std::tuple<uintE, pbbs::empty>*)new_edges.to_array();
  return symmetric_graph<symmetricVertex, pbbs::empty>(vv, n, new_m, get_deletion_fn(vv, edge_arr), edge_arr);
}


template <class G, class S>
inline symmetric_graph<csv_bytepd_amortized, pbbs::empty> renumber_compressed(G& GA, S& perm) {
  using W = typename G::weight_type;
  uintE n = GA.n;

  // build the new graph
  auto byte_offsets = sequence<size_t>(n+1);
  auto degrees = sequence<size_t>(n+1);
  parallel_for(0, n, [&] (size_t i) {
    uintE our_id = perm[i];
    uintE our_stk[1000];
    uintE* tmp = (uintE*)our_stk;
    if (GA.get_vertex(i).getOutDegree() > 1000) {
      tmp = pbbs::new_array_no_init<uintE>(GA.get_vertex(i).getOutDegree());
    }
    size_t ct = 0;
    auto map_f = [&](uintE src, uintE ngh, const W& wgh) {
      uintE ngh_id = perm[ngh];
      if (ngh_id < our_id) {
        tmp[ct] = ngh_id;
        ct++;
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f, false);
    degrees[our_id] = ct;
    auto seq = pbbslib::make_sequence(tmp, ct);
    pbbs::sample_sort_inplace(seq, std::less<uintE>());

    // map over the new neighbors, compute their size
    size_t total_bytes = 0;
    uintE last_ngh = 0;
    size_t deg = 0;
    uchar bytetmp[16];
    uintE u = our_id;
    for (size_t j=0; j<ct; j++) {
      long bytes = 0;
      uintE v = tmp[j];
      if ((j % PARALLEL_DEGREE) == 0) {
        bytes = bytepd_amortized::compressFirstEdge(bytetmp, bytes, u, v);
      } else {
        bytes = bytepd_amortized::compressEdge(bytetmp, bytes, v - last_ngh);
        if (last_ngh == v) {
          cout << "last_ngh = " << last_ngh << " v = " << v << " j = " << j << " ct = " << ct << " source id = " << i << endl;
          cout << "duplicate!!!" << endl;
        }
      }
      last_ngh = v;
      total_bytes += bytes;
    }
    if (ct > 0) {
      size_t n_chunks = 1+(ct-1)/PARALLEL_DEGREE;
      // To account for the byte offsets
      total_bytes += (n_chunks-1)*sizeof(uintE);
      // To account for the per-block counters
      total_bytes += (n_chunks)*sizeof(uintE);
      // To account for the virtual degree
      total_bytes += sizeof(uintE);
    }

    byte_offsets[our_id] = total_bytes;

    if (GA.get_vertex(i).getOutDegree() > 1000) {
      pbbs::free_array(tmp);
    }
  }, 1);

  byte_offsets[n] = 0;
  degrees[n] = 0;
  size_t total_bytes = pbbslib::scan_add_inplace(byte_offsets.slice());
  size_t new_m = pbbslib::scan_add_inplace(degrees.slice());
  cout << "new total bytes = " << total_bytes << endl;
  cout << "new m = " << new_m << endl;

  uchar* edges = pbbs::new_array_no_init<uchar>(total_bytes);

  parallel_for(0, n, [&] (size_t i) {
    uintE our_id = perm[i];
    size_t our_offset = byte_offsets[our_id];
    size_t deg = degrees[our_id+1]-degrees[our_id];

    if (deg > 0) {
      uintE our_stk[1000];
      uintE* tmp = (uintE*)our_stk;
      if (deg > 1000) {
        tmp = pbbs::new_array_no_init<uintE>(deg);
      }
      size_t ct = 0;
      auto map_f = [&](uintE src, uintE ngh, const W& wgh) {
        uintE ngh_id = perm[ngh];
        if (ngh_id < our_id) {
          tmp[ct] = ngh_id;
          ct++;
        }
      };
      GA.get_vertex(i).mapOutNgh(i, map_f, false); // original id
      assert(ct == deg);

      auto seq = pbbslib::make_sequence(tmp, ct);
      pbbs::sample_sort_inplace(seq, std::less<uintE>());

      auto it = vertex_ops::get_iter<pbbs::empty>((std::tuple<uintE, pbbs::empty>*)tmp, ct);

      long nbytes = bytepd_amortized::sequentialCompressEdgeSet<pbbs::empty>(edges + our_offset, 0, ct, our_id, it, PARALLEL_DEGREE);

      long expected_bytes = (byte_offsets[our_id+1] - byte_offsets[our_id]);
      if (nbytes != expected_bytes) {
        cout << "nbytes = " << nbytes << ". Should be: " << expected_bytes << " deg = " << deg << " vtxid = " << our_id << endl;
        assert(nbytes == expected_bytes);
        // exit(0);
      }

      if (deg > 1000) {
        pbbs::free_array(tmp);
      }
    }
  }, 1);

  auto v = pbbs::sequence<vertex_data>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE new_id = perm[i];
    size_t o = byte_offsets[new_id];
    size_t degree = degrees[new_id+1]-degrees[new_id];
    v[new_id].degree = degree;
    v[new_id].offset = o;
  });
  auto vv = v.to_array();
  return symmetric_graph<csv_bytepd_amortized, pbbs::empty>(vv, n, new_m, get_deletion_fn(vv, edges), edges);
}

template <class G>
std::tuple<bool, size_t> YoshidaCache(G& GA, uintE u, std::unordered_map<uintE, bool>& cache) {
  using W = typename G::weight_type;
  bool in_mis = true;
  size_t total_work = 1;
  bool quit = false;
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) -> bool {
    bool ngh_in = false;
    size_t ngh_total_work = 0;
    auto lookup = cache.find(v);
    if (lookup == cache.end()) {
      std::tie(ngh_in, ngh_total_work) = YoshidaCache(GA, v, cache);
      total_work += ngh_total_work + 1; // recursive work done below, plus one for this neighbor
      cache.insert(std::make_pair(v, ngh_in));
    } else {
      ngh_in = lookup->second;
      total_work += 1; // 1 for this neighbor (pruned)
    }
    if (ngh_in) {
      in_mis = false;
      quit = true;
      return false; // break
    } else {
      return true; // keep going, ngh was not in MIS
    }
  };
  if (GA.get_vertex(u).getOutDegree() > 0) {
    GA.get_vertex(u).decodeOutNghSparseSeqCond(u, map_f);
  } // otherwise base case; return that we are in mis, and total_work is 1 (to check this vertex)
  return std::make_tuple(in_mis, total_work);
}

template <class G>
std::tuple<bool, size_t> Yoshida(G& GA, uintE u) {
  using W = typename G::weight_type;
  bool in_mis = true;
  size_t total_work = 1;
  auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) -> bool {
    assert(v < u);
    bool ngh_in = false;
    size_t ngh_total_work = 0;
    std::tie(ngh_in, ngh_total_work) = Yoshida(GA, v);
    total_work += ngh_total_work + 1; // recursive work done below, plus one for this neighbor
    if (ngh_in) {
      in_mis = false;
      return false; // break
    } else {
      return true; // keep going, ngh was not in MIS
    }
  };
  if (GA.get_vertex(u).getOutDegree() > 0) {
    GA.get_vertex(u).decodeOutNghSparseSeqCond(u, map_f);
  } // otherwise base case; return that we are in mis, and total_work is 1 (to check this vertex)
  return std::make_tuple(in_mis, total_work);
}

template <class G>
auto YoshidaMIS(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto in_mis = pbbs::sequence<bool>(n, false);

  auto perm = pbbslib::random_permutation<uintE>(n);

  auto GR = renumber_compressed(GA, perm);

  timer search_t; search_t.start();
  cout << "GR.n = " << GR.n << " GR.m = " << GR.m << endl;
  auto search_lengths = pbbs::sequence<size_t>(n, static_cast<size_t>(0));

  auto inv_map = pbbs::sequence<uintE>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE id = i;
    uintE new_id = perm[i];
    inv_map[new_id] = i;
  });

  auto xors = pbbs::sequence<uintE>(n);
  parallel_for(0, n, [&] (size_t i) {
      uintE xorr = 0;
      auto f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        xorr ^= v;
        return true;
      };
    GR.get_vertex(i).decodeOutNghSparseSeqCond(i, f);
    xors[i] = xorr;
  }, 1);
  cout << "xor of graph = " << pbbslib::reduce_xor(xors) << endl;

  size_t n_batches = 10;
  size_t bs = pbbs::num_blocks(n, 5);
  n_batches = pbbs::num_blocks(n, bs);
  for (size_t i=0; i<n_batches; i++) {
    size_t start = bs*i;
    size_t end = std::min(bs*(i+1), (size_t)n);
    cout << "running start = " << start << " to " << end << endl;
    parallel_for(start, end, [&] (size_t i) {
      bool is_in_mis = false;
      size_t search_length = 0;

      std::tie(is_in_mis, search_length) = Yoshida(GR, i);
      uintE orig_id = inv_map[i];
      in_mis[orig_id] = is_in_mis;
      search_lengths[orig_id] = search_length;
    }, 512);
  }

  search_t.stop(); search_t.reportTotal("yoshida search time");

  size_t total_work = pbbslib::reduce_add(search_lengths);
  cout << "total work = " << total_work << endl;

  size_t max_work = pbbslib::reduce_max(search_lengths);
  cout << "max search length = " << max_work << endl;

  pbbs::sample_sort_inplace(search_lengths.slice(), std::less<size_t>());
  cout << "median search length = " << search_lengths[n/2] << endl;

  cout << "avg search length = " << (static_cast<double>(total_work)/n) << endl;

//  MIS_rootset::verify_mis(GA, in_mis);
//  cout << "mis size = " << mis_size << endl;
}

template <class G>
auto YoshidaMIS_2(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto in_mis = pbbs::sequence<bool>(n, false);

  auto perm = pbbslib::random_permutation<uintE>(n);

  auto GR = renumber(GA, perm);

  timer search_t; search_t.start();
  cout << "GR.n = " << GR.n << " GR.m = " << GR.m << endl;
  auto search_lengths = pbbs::sequence<size_t>(n, static_cast<size_t>(0));

  auto inv_map = pbbs::sequence<uintE>(n);
  parallel_for(0, n, [&] (size_t i) {
    uintE id = i;
    uintE new_id = perm[i];
    inv_map[new_id] = i;
  });

  auto xors = pbbs::sequence<uintE>(n);
  parallel_for(0, n, [&] (size_t i) {
      uintE xorr = 0;
      auto f = [&] (const uintE& u, const uintE& v, const W& wgh) {
        xorr ^= v;
        return true;
      };
    GR.get_vertex(i).decodeOutNghSparseSeqCond(i, f);
    xors[i] = xorr;
  }, 1);
  cout << "xor of graph = " << pbbslib::reduce_xor(xors) << endl;

  size_t n_batches = 10;
  size_t bs = pbbs::num_blocks(n, 5);
  n_batches = pbbs::num_blocks(n, bs);
  for (size_t i=0; i<n_batches; i++) {
    size_t start = bs*i;
    size_t end = std::min(bs*(i+1), (size_t)n);
    cout << "running start = " << start << " to " << end << endl;
    parallel_for(start, end, [&] (size_t i) {
      bool is_in_mis = false;
      size_t search_length = 0;

      std::tie(is_in_mis, search_length) = Yoshida(GR, i);
      uintE orig_id = inv_map[i];
      in_mis[orig_id] = is_in_mis;
      search_lengths[orig_id] = search_length;
    }, 512);
  }

  search_t.stop(); search_t.reportTotal("yoshida search time");

  size_t total_work = pbbslib::reduce_add(search_lengths);
  cout << "total work = " << total_work << endl;

  size_t max_work = pbbslib::reduce_max(search_lengths);
  cout << "max search length = " << max_work << endl;

  pbbs::sample_sort_inplace(search_lengths.slice(), std::less<size_t>());
  cout << "median search length = " << search_lengths[n/2] << endl;

  cout << "avg search length = " << (static_cast<double>(total_work)/n) << endl;

//  MIS_rootset::verify_mis(GA, in_mis);
//  cout << "mis size = " << mis_size << endl;
}
