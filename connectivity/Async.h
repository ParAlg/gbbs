// Based on AsyncST/CC.C
#pragma once

#include <algorithm>
#include "ligra.h"
#include "utils/stats.h"

namespace async {
inline uintE find(uintE i, pbbs::sequence<uintE>& parent) {
  while(i != parent[i])
    i = parent[i];
  return i;
}

inline uintE find_compress(uintE i, pbbs::sequence<uintE>& parent) {
  uintE j = i; //parent[i];
  if (parent[j] == j) return j;
  do j = parent[j];
  while (parent[j] != j);
  //note: path compression can happen in parallel in the same tree, so
  //only link from smaller to larger to avoid cycles
  uintE tmp;
  while((tmp=parent[i])<j){ parent[i]=j; i=tmp;}
  return j;
}

inline void find_check(uintE i, pbbs::sequence<uintE>& parent) {
  while(i < parent[i]) {
    i = parent[i];
  }
  if(i > parent[i]) { cout << "error\n"; exit(0); }
}

inline uintE find_split(uintE i, pbbs::sequence<uintE>& parent) {
  while(1) {
    uintE v = parent[i];
    uintE w = parent[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parent[i],v,w);
      i = v;
    }
  }
}

inline uintE find_halve(uintE i, pbbs::sequence<uintE>& parent) {
  while(1) {
    uintE v = parent[i];
    uintE w = parent[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parent[i],v,w);
      //i = w;
      i = parent[i];
    }
  }
}

inline void uniteEarly(uintE u, uintE v, pbbs::sequence<uintE>& parents) {
  while(1) {
    if(u == v) return;
    if(v < u) std::swap(u,v);
    if(pbbs::atomic_compare_and_swap(&parents[u],u,v)) { return; }
    uintE z = parents[u];
    uintE w = parents[z];
    pbbs::atomic_compare_and_swap(&parents[u],z,w);
    u = z;
  }
}
} // namespace async

namespace async_cc {
inline void unite(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents) {
  uintE u = u_orig;
  uintE v = v_orig;
  while(1) {
    u = async::find_compress(u,parents);
    v = async::find_compress(v,parents);
    if(u == v) break;
    else if (u < v && parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) {
      return;
    }
    else if (v < u && parents[v] == v && pbbs::atomic_compare_and_swap(&parents[v],v,u)) {
      return;
    }
  }
}
} // namespace async_cc

// returns a component labeling
template <class G>
pbbs::sequence<uintE> AsyncCC(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        async_cc::unite(u, v, parents);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  });

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = async::find_compress(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  // filter UINT_E_MAX from hooks if spanning forest edges are desired.
  return parents;
}

// does a compare-and-swap on the parents array.
// ndSTOpt does CAS on the hooks array.
namespace async_sf {
inline void unite(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents, pbbs::sequence<uintE>& hooks) {
  uintE u = u_orig;
  uintE v = v_orig;
  while(1) {
    u = async::find_compress(u,parents);
    v = async::find_compress(v,parents);
    if(u == v) break;
    else if (u < v && parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) {
      hooks[u] = v;
      return;
    }
    else if (v < u && parents[v] == v && pbbs::atomic_compare_and_swap(&parents[v],v,u)) {
      hooks[v] = u;
      return;
    }
  }
}

inline void uniteEarly(uintE u, uintE v, pbbs::sequence<uintE>& parents, pbbs::sequence<uintE>& hooks) {
  while(1) {
    if(u == v) return;
    if(v < u) std::swap(u,v);
    if(pbbs::atomic_compare_and_swap(&parents[u],u,v)) { hooks[u] = v; return;}
    uintE z = parents[u];
    uintE w = parents[z];
    pbbs::atomic_compare_and_swap(&parents[u],z,w);
    u = z;
  }
}
} // namespace async_sf

// returns edge pairs
template <class G>
pbbs::sequence<std::pair<uintE, uintE>> AsyncSF(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  auto hooks = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        async_sf::unite(u, v, parents, hooks);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  });

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = async::find_compress(i,parents);
  });
  ft.stop(); ft.reportTotal("find time");
  using edge = std::pair<uintE, uintE>;
  auto edge_seq = pbbs::delayed_seq<edge>(n, [&] (size_t i) {
      return std::make_pair(i, hooks[i]);
  });
  auto sf_edges = pbbs::filter(edge_seq, [&] (const edge& e) { return e.second != UINT_E_MAX; });
  return sf_edges;
}



