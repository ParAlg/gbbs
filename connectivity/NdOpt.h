// Based on ndSTOpt.C
#pragma once

#include <algorithm>
#include "ligra.h"
#include "utils/stats.h"

namespace nd {
  constexpr const size_t line_width = 16;
  inline uintT find_wide(uintE i, pbbs::sequence<uintE>& parent) {
    //if (parent[i] < 0) return i;
    uintT j = i; //parent[i];
    if (parent[j*line_width] == UINT_E_MAX) return j;
    do j = parent[j*line_width];
    while (parent[j*line_width] < UINT_E_MAX);
    //note: path compression can happen in parallel in the same tree, so
    //only link from smaller to larger to avoid cycles
    uintE tmp;
    while((tmp=parent[i*line_width])<j){ parent[i*line_width]=j; i=tmp;}
    return j;
  }


  // Assumes root is negative
  // Not making parent array volatile improves
  // performance and doesn't affect correctness
  inline uintT find(uintE i, pbbs::sequence<uintE>& parent) {
    //if (parent[i] < 0) return i;
    uintT j = i; //parent[i];
    if (parent[j] == UINT_E_MAX) return j;
    do j = parent[j];
    while (parent[j] < UINT_E_MAX);
    //note: path compression can happen in parallel in the same tree, so
    //only link from smaller to larger to avoid cycles
    uintE tmp;
    while((tmp=parent[i])<j){ parent[i]=j; i=tmp;}
    return j;
  }



//  inline void unite(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents, pbbs::sequence<uintE>& hooks) {
//   while (1) {
//     uintE u = nd::find_wide(ic,parents);
//     uintE v = nd::find_wide(vc,parents);
//     if(u == v) break;
//     if(u > v) std::swap(u,v);
//     //if successful, store the ID of the edge used in hooks[u]
//     if(hooks[u*nd::line_width] == UINT_E_MAX && __sync_bool_compare_and_swap(&hooks[u*nd::line_width],UINT_E_MAX,1)) {
//       parents[u*nd::line_width]=v;
//       break;
//     }
//   }
//  }
} // namespace nd

template <class G>
pbbs::sequence<uintE> NdOptCC(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });
  auto hooks = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u_orig, uintE v_orig, const W& wgh) {
      if (u_orig < v_orig) {
//        nd::unite(u_orig, v_orig, parents, hooks);
        while (1) {
          uintE u = nd::find(u_orig,parents);
          uintE v = nd::find(v_orig,parents);
          if(u == v) break;
          if(u > v) std::swap(u,v);
          //if successful, store the ID of the edge used in hooks[u]
          if(hooks[u] == UINT_E_MAX && __sync_bool_compare_and_swap(&hooks[u],UINT_E_MAX,i)) {
            parents[u]=v;
            break;
          }
        }
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel over nghs
  }, 128);

  //uncomment this code for connected component labels instead of
  //spanning forest edges
  parallel_for(0, n, [&] (size_t i) {
    hooks[i] = nd::find(i,parents);
  });
  return hooks;
}


template <class G>
pbbs::sequence<std::pair<uintE, uintE>> NdOptSF(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });
  auto hooks = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u_orig, uintE v_orig, const W& wgh) {
      if (u_orig < v_orig) {
        while (1) {
          uintE u = nd::find(u_orig,parents);
          uintE v = nd::find(v_orig,parents);
          if(u == v) break;
          if(u > v) std::swap(u,v);
          //if successful, store the ID of the edge used in hooks[u]
          if(hooks[u] == UINT_E_MAX && __sync_bool_compare_and_swap(&hooks[u],UINT_E_MAX,i)) {
    	  parents[u]=v;
    	  break;
          }
        }
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel over nghs
  });

  using edge = std::pair<uintE, uintE>;
  auto edge_seq = pbbs::delayed_seq<edge>(n, [&] (size_t i) {
      return std::make_pair(i, hooks[i]);
  });
  auto sf_edges = pbbs::filter(edge_seq, [&] (const edge& e) { return e.second != UINT_E_MAX; });
  return sf_edges;

//  //get the IDs of the edges in the spanning forest
//  _seq<uintT> stIdx = sequence::filter((uintT*) hooks, n, notMax());
//
//  free(parents); free(hooks);
//  cout<<"nInSt = "<<stIdx.n<<endl;
//  return pair<uintT*,uintT>(stIdx.A, stIdx.n);
}


