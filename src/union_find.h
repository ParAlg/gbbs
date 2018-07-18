#pragma once

#include "lib/index_map.h"
#include "lib/speculative_for.h"

struct UnionFind {
  size_t n;
  intT* parents;
  UnionFind(size_t _n) :
    n(_n) {
    parents = newA(intT, n);
    parallel_for(size_t i=0; i<n; i++) {
      parents[i] = -1;
    }
  }

  intT find(int32_t i) {
    if (parents[i] < 0) return i;
    intT j = parents[i];
    if (parents[j] < 0) return j;
    do j = parents[j];
    while (parents[j] >= 0);
    intT tmp;
    while ((tmp = parents[i]) != j) {
      parents[i] = j;
      i = tmp;
    }
    return j;
  }

  void link(intT u, intT v) {
    parents[u] = v;
  }

  void del() {
    free(parents);
  }

};

// edges: <uintE, uintE, W>
template <class intT, class Edges, class ST, class UF>
struct UnionFindStep {
  using res = reservation<intT>;
  using storage = tuple<intT, intT>;
  Edges& E;
  res* R;
  ST& inST;
  UF& uf;
  storage* indices;

  size_t n;

  void del() {
    free(indices);
  }

  UnionFindStep(Edges& _E, res* _R, ST& ist, UF& _uf) :
    E(_E), R(_R), inST(ist), uf(_uf) {
    n = uf.n;
    indices = newA(storage, E.non_zeros);
  }

  bool reserve(intT i) {
    assert(i < E.non_zeros);
    auto e = E.E[i];
    intT u = uf.find(get<0>(e));
    intT v = uf.find(get<1>(e));
    if (u != v) {
      indices[i] = make_tuple(u, v);
      R[v].reserve(i);
      R[u].reserve(i);
      assert(u < n); assert(v < n);
      return 1; // active
    } else return 0; // done
  }

  bool commit(intT i) {
    assert(i < E.non_zeros);
    // read back u and v from 'reserve'
    auto st = indices[i];
    intT u = get<0>(st), v = get<1>(st);
    if (u >= n || v >= n) {
      cout << "u = " << u << " v = " << v << " i = " << i << endl;
      exit(0);
    }
//    assert(u < n); assert(v < n);
    if (R[v].checkReset(i)) {
      R[u].checkReset(i);
      uf.link(v, u);
      inST[i] = 1;
      return 1;}
    else if (R[u].checkReset(i)) {
      uf.link(u, v);
      inST[i] = 1;
      return 1; }
    else return 0;
  }
};

template <class intT, class Edges, class R, class ST, class UF>
auto make_uf_step(Edges& e, R r, ST& ist, UF& uf) {
  return UnionFindStep<intT, Edges, ST, UF>(e, r, ist, uf);
}

