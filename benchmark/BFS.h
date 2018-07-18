#pragma once

#include "ligra.h"

template <class W>
struct BFS_F {
  uintE* Parents;
  BFS_F(uintE* _Parents) : Parents(_Parents) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (Parents[d] == UINT_E_MAX) {
      Parents[d] = s;
      return 1;
    } else
      return 0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    return (CAS(&Parents[d], UINT_E_MAX, s));
  }
  inline bool cond(const uintE& d) { return (Parents[d] == UINT_E_MAX); }
};

template <template <class W> class vertex, class W>
auto BFS(graph<vertex<W> >& GA, uintE src) {
  using w_vertex = vertex<W>;

  // Creates Parents array, initialized to all -1, except for src.
  auto Parents = array_imap<uintE>(GA.n, [&](size_t i) { return UINT_E_MAX; });
  Parents[src] = src;

  vertexSubset Frontier(GA.n, src);
  size_t reachable = 0;
  while (!Frontier.isEmpty()) {
    cout << Frontier.size() << endl;
    reachable += Frontier.size();
    vertexSubset output =
        edgeMap(GA, Frontier, BFS_F<W>(Parents.start()), -1, sparse_blocked);
    Frontier.del();
    Frontier = output;
  }
  Frontier.del();
  cout << "Reachable: " << reachable << endl;
  return Parents;
}
