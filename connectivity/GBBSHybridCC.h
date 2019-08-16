// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "../benchmark/LDD.h"
#include "ligra.h"
#include "utils/contract.h"
#include "utils/stats.h"
#include "Async.h"
#include "NdOpt.h"
#include "Rem.h"

namespace gbbs_hybridcc {

//template <template <class W> class vertex, class W>
//auto UFContractOriginal(graph<vertex<W>>& GA, pbbs::sequence<uintE>& clusters) {
//  size_t n = GA.n;
//  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; }); // UINT_E_MAX; });
////  auto hooks = pbbs::sequence<uintE>(num_clusters*nd::line_width, [&] (size_t i) { return UINT_E_MAX; });
//  auto locks = pbbs::sequence<std::mutex>(n);
//
//  timer tt; tt.start();
//  // Most of the time is spent just fetching clusters[v]. Can verify with perf.
//  parallel_for(0, n, [&] (size_t i) {
//    uintE ic = clusters[i];
//    auto map_f = [&] (uintE u_orig, uintE v_orig, const W& wgh) {
//      if (u_orig < v_orig) {
//        uintE vc = clusters[v_orig];
//        if (ic != vc) {
//          remcc::unite(ic, vc, parents, locks);
////          async_cc::unite(ic, vc, parents);
////          nd::unite(ic, vc, parents, hooks);
////          while (1) {
////            uintE u = nd::find_wide(ic,parents);
////            uintE v = nd::find_wide(vc,parents);
////            if(u == v) break;
////            if(u > v) std::swap(u,v);
////            //if successful, store the ID of the edge used in hooks[u]
////            if(hooks[u*nd::line_width] == UINT_E_MAX && __sync_bool_compare_and_swap(&hooks[u*nd::line_width],UINT_E_MAX,i)) {
////              parents[u*nd::line_width]=v;
////              break;
////            }
////          }
//        }
//      }
//    };
//    GA.V[i].mapOutNgh(i, map_f); // in parallel
//  }, 128);
//  tt.stop(); tt.reportTotal("traverse edges time");
//
//  parallel_for(0, n, [&] (size_t i) {
////    hooks[i] = nd::find(i, parents);
//    parents[i] = async::find_compress(i,parents);
//  });
//
//  return parents;
////  return hooks;
//}


template <class G>
auto UFContracted(G& GA, pbbs::sequence<uintE>& clusters, size_t num_clusters) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  auto parents = pbbs::sequence<uintE>(num_clusters, [&] (size_t i) { return i; }); // UINT_E_MAX; });
//  auto hooks = pbbs::sequence<uintE>(num_clusters*nd::line_width, [&] (size_t i) { return UINT_E_MAX; });
  auto locks = pbbs::sequence<std::mutex>(num_clusters);
  size_t threshold = GA.n / 45000000;
  cout << "threshold = " << threshold << " num_clusters = " << num_clusters << endl;

  timer tt; tt.start();
  // Most of the time is spent just fetching clusters[v]. Can verify with perf.
  parallel_for(0, n, [&] (size_t i) {
    uintE ic = clusters[i];
    auto map_f = [&] (uintE u_orig, uintE v_orig, const W& wgh) {
      if (u_orig < v_orig) {
        uintE vc = clusters[v_orig];
        if (ic != vc) {
          remcc::unite(ic, vc, parents, locks);
//          async_cc::unite(ic, vc, parents);
//          nd::unite(ic, vc, parents, hooks);
//          while (1) {
//            uintE u = nd::find_wide(ic,parents);
//            uintE v = nd::find_wide(vc,parents);
//            if(u == v) break;
//            if(u > v) std::swap(u,v);
//            //if successful, store the ID of the edge used in hooks[u]
//            if(hooks[u*nd::line_width] == UINT_E_MAX && __sync_bool_compare_and_swap(&hooks[u*nd::line_width],UINT_E_MAX,i)) {
//              parents[u*nd::line_width]=v;
//              break;
//            }
//          }
        }
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  }, 128);
  tt.stop(); tt.reportTotal("traverse edges time");

  parallel_for(0, num_clusters, [&] (size_t i) {
//    hooks[i] = nd::find(i, parents);
    parents[i] = async::find_compress(i,parents);
  });

  return parents;
//  return hooks;
}

template <class G>
inline sequence<uintE> CC_impl(G& GA, double beta,
                                 size_t level, bool permute = false) {
  size_t n = GA.n;
  permute |= (level > 0);
  timer ldd_t;
  ldd_t.start();
  auto clusters = LDD(GA, beta, permute);
  ldd_t.stop();
  debug(ldd_t.reportTotal("ldd time"));

  timer relabel_t;
  relabel_t.start();
  size_t num_clusters = contract::RelabelIds(clusters);
  relabel_t.stop();
  debug(relabel_t.reportTotal("relabel time"));

  // contract based on UF. One pass through all edges. Apply UF to each IC edge.
  timer uf; uf.start();
  auto cluster_comps = UFContracted(GA, clusters, num_clusters);
//  auto cluster_comps = UFContractOriginal(GA, clusters);
  uf.stop(); debug(uf.reportTotal("union find time"));

  parallel_for(0, n, [&] (size_t i) {
    uintE c_i = clusters[i];
    uintE cc_i = cluster_comps[c_i];
    clusters[i] = cc_i;
  });
  cout << "returned clusters" << endl;

  return clusters;
}

template <class G>
inline sequence<uintE> CC(G& GA, double beta = 0.2, bool permute = false) {
  return CC_impl(GA, beta, 0, permute);
}

}  // namespace gbbs_hybrid
