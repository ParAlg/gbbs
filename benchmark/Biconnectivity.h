#pragma once

#include "CC.h"
#include "lib/dyn_arr.h"
#include "lib/random.h"
#include "lib/sample_sort.h"
#include "lib/sparse_table.h"
#include "lib/speculative_for.h"
#include "oldlib/benchIO.h"
#include "ligra.h"

namespace bc {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;
};  // namespace bc

using labels = tuple<uintE, uintE>;

// Used to compute the size of each subtree (leaffix)
struct AugF {
  uintE* sizes;
  intE* cts;
  AugF(uintE* _sizes, intE* _cts) : sizes(_sizes), cts(_cts) {}
  bool update(const uintE& s, const uintE& d) {
    sizes[d] += sizes[s];
    cts[d]--;
    if (!cts[d]) {
      return true;
    }
    return false;
  }
  bool updateAtomic(const uintE& s, const uintE& d) {
    writeAdd(&sizes[d], sizes[s]);
    intE res = writeAdd(&cts[d], -1);
    return (res == 0);
  }
  bool cond(const uintE& d) { return true; }
};

// Used to compute the (low, hi) values of each subtree (leaffix)
struct MinMaxF {
  labels* MM;  // min and max
  intE* cts;
  MinMaxF(labels* _MM, intE* _cts) : MM(_MM), cts(_cts) {}
  bool update(const uintE& s, const uintE& d) {
    labels lab_s = MM[s];
    uintE low = get<0>(lab_s), high = get<1>(lab_s);
    if (low < get<0>(MM[d])) {
      get<0>(MM[d]) = low;
    }
    if (high > get<1>(MM[d])) {
      get<1>(MM[d]) = high;
    }
    uintE ct = cts[d] - 1;
    cts[d] = ct;
    return ct == 0;
  }
  bool updateAtomic(const uintE& s, const uintE& d) {
    labels lab_s = MM[s];
    uintE low = get<0>(lab_s), high = get<1>(lab_s);
    if (low < get<0>(MM[d])) {
      writeMin(&get<0>(MM[d]), low);
    }
    if (high > get<1>(MM[d])) {
      writeMax(&get<1>(MM[d]), high);
    }
    intE res = writeAdd(&cts[d], -1);
    return (res == 0);
  }
  bool cond(const uintE& d) { return true; }
};

template <template <typename W> class vertex, class W, class Seq>
auto preorder_number(graph<vertex<W>>& GA, uintE* Parents, Seq& Sources) {
  using w_vertex = vertex<W>;
  size_t n = GA.n;
  using edge = tuple<uintE, uintE>;
  auto out_edges = array_imap<edge>(
      n, [](size_t i) { return make_tuple(UINT_E_MAX, UINT_E_MAX); });
  parallel_for(size_t i = 0; i < n; i++) {
    uintE p_i = Parents[i];
    if (p_i != i) {
      out_edges[i] = make_tuple(p_i, i);
    }
  }

  auto edges = pbbs::filter(
      out_edges, [](const edge& e) { return get<0>(e) != UINT_E_MAX; });
  out_edges.del();
  auto sort_tup = [](const edge& l, const edge& r) { return l < r; };
  pbbs::sample_sort(edges.start(), edges.size(), sort_tup);

  auto starts = array_imap<uintE>(n + 1, [](size_t i) { return UINT_E_MAX; });
  parallel_for(size_t i = 0; i < edges.size(); i++) {
    if (i == 0 || get<0>(edges[i]) != get<0>(edges[i - 1])) {
      starts[get<0>(edges[i])] = i;
    }
  }
  starts[n] = edges.size();

  timer seq;
  seq.start();
  for (long i = starts.size() - 1; i >= 0; i--) {
    if (starts[i] == UINT_E_MAX) {
      starts[i] = starts[i + 1];
    }
  }
  seq.stop();
  seq.reportTotal("seq time");

  auto nghs = array_imap<uintE>(edges.size(),
                                [&](size_t i) { return get<1>(edges[i]); });
  edges.del();

  timer augs;
  augs.start();
  // Create directed BFS tree
  using vtx = asymmetricVertex<pbbs::empty>;
  auto v = newA(vtx, n);
  parallel_for(size_t i = 0; i < n; i++) {
    uintE out_off = starts[i];
    uintE out_deg = starts[i + 1] - out_off;
    v[i].setOutDegree(out_deg);
    v[i].setOutNeighbors((tuple<uintE, pbbs::empty>*)(nghs.start() + out_off));
    uintE parent = Parents[i];
    if (parent != i) {
      v[i].setInDegree(1);
      v[i].setInNeighbors(((tuple<uintE, pbbs::empty>*)(Parents + i)));
    } else {
      v[i].setInDegree(0);
      v[i].setInNeighbors(nullptr);
    }
  }
  auto Tree = graph<vtx>(v, n, nghs.size(), []() {});

  // 1. Leaffix for Augmented Sizes
  auto aug_sizes = array_imap<uintE>(n, [](size_t i) { return 1; });
  auto cts =
      array_imap<intE>(n, [&](size_t i) { return Tree.V[i].getOutDegree(); });
  auto leaf_im = make_in_imap<bool>(n, [&](size_t i) {
    auto s_i = starts[i];
    auto s_n = starts[i + 1];
    size_t deg = s_n - s_i;
    return (Parents[i] != i) && deg == 0;
  });
  auto leafs = pbbs::pack_index<uintE>(leaf_im);

  auto vs = vertexSubset(n, leafs.size(), leafs.start());
  size_t rds = 0, tv = 0;
  while (!vs.isEmpty()) {
    rds++;
    size_t et;
    tv += vs.size();
    // histogram or write-add parents, produce next em.
    auto output = edgeMap(
        Tree, vs, wrap_em_f<pbbs::empty>(AugF(aug_sizes.start(), cts.start())),
        -1, in_edges | sparse_blocked);
    if (rds > 1) {
      vs.del();  // don't delete leafs.
    }
    vs = output;
    output.toSparse();
  }
  augs.stop();
  augs.reportTotal("aug size time");

  // Optional: Prefix sum over sources with aug-size (if we want distinct
  // preorder #s)
  auto s_copy = Sources.copy();

  timer pren;
  pren.start();
  auto PN = array_imap<uintE>(n);
  vs = vertexSubset(n, s_copy.size(), s_copy.get_array());
  parallel_for(size_t i = 0; i < Sources.size(); i++) {
    uintE v = s_copy[i];
    PN[v] = 0;
  }
  rds = 0;
  tv = 0;
  while (!vs.isEmpty()) {
    rds++;
    tv += vs.size();
    size_t et;
    auto offsets = array_imap<uintE>(vs.size(), [&](size_t i) {
      uintE v = vs.s[i];
      return Tree.V[v].getOutDegree();
    });
    auto tot = pbbs::scan_add(offsets, offsets);
    auto next_vs = array_imap<uintE>(tot);
    parallel_for(size_t i = 0; i < vs.size(); i++) {
      uintE v = vs.s[i];
      uintE off = offsets[i];
      uintE deg_v = Tree.V[v].getOutDegree();
      uintE preorder_number = PN[v] + 1;

      if (deg_v < 4000) {
        // Min and max in any vertex are [PN[v], PN[v] + aug_sizes[v])
        for (size_t j = 0; j < deg_v; j++) {
          uintE ngh = Tree.V[v].getOutNeighbor(j);
          PN[ngh] = preorder_number;
          preorder_number += aug_sizes[ngh];
          next_vs[off + j] = ngh;
        }
      } else {
        auto A = array_imap<uintE>(deg_v);
        parallel_for(size_t j = 0; j < deg_v; j++) {
          uintE ngh = Tree.V[v].getOutNeighbor(j);
          A[j] = aug_sizes[ngh];
        }
        pbbs::scan_add(A, A);
        parallel_for(size_t j = 0; j < deg_v; j++) {
          uintE ngh = Tree.V[v].getOutNeighbor(j);
          uintE pn = preorder_number + A[j];
          PN[ngh] = pn;
          next_vs[off + j] = ngh;
        }
      }
    }
    vs.del();
    vs = vertexSubset(n, next_vs.size(), next_vs.get_array());
  }
  pren.stop();
  pren.reportTotal("preorder number from sizes time");

  timer map_e;
  map_e.start();
  // Map all edges and update labels.
  auto MM = array_imap<labels>(n, [&](size_t i) {
    uintE pn_i = PN[i];
    return make_tuple(pn_i, pn_i + aug_sizes[i] - 1);
  });
  // Note that we have to exclude tree edges: it's all edges spliced in with
  // the exception of tree edges.
  auto map_f = wrap_f<W>([&](const uintE& u, const uintE& v) {
    if (u < v) {  // see if assigning roles helps temp loc later.
      if (u == Parents[v] || v == Parents[u]) {
        return;
      }
      uintE p_u = PN[u];
      uintE p_v = PN[v];
      if (p_u < get<0>(MM[v])) {
        writeMin(&get<0>(MM[v]), p_u);
      } else if (p_u > get<1>(MM[v])) {
        writeMax(&get<1>(MM[v]), p_u);
      }
      if (p_v < get<0>(MM[u])) {
        writeMin(&get<0>(MM[u]), p_v);
      } else if (p_v > get<1>(MM[u])) {
        writeMax(&get<1>(MM[u]), p_v);
      }
    }
  });
  parallel_for(size_t i = 0; i < n; i++) { GA.V[i].mapOutNgh(i, map_f); }
  map_e.stop();
  map_e.reportTotal("map edges time");

  timer leaff;
  leaff.start();
  // 1. Leaffix to update min/max
  parallel_for(size_t i = 0; i < n; i++) { cts[i] = Tree.V[i].getOutDegree(); }

  vs = vertexSubset(n, leafs.size(), leafs.get_array());
  rds = 0, tv = 0;
  while (!vs.isEmpty()) {
    rds++;
    tv += vs.size();
    // histogram or write-add parents, produce next em.
    auto output = edgeMap(
        Tree, vs, wrap_em_f<pbbs::empty>(MinMaxF(MM.start(), cts.start())), -1,
        in_edges | sparse_blocked);
    vs.del();
    vs = output;
  }
  // Delete tree
  free(v);
  nghs.del();
  leaff.stop();
  leaff.reportTotal("leaffix to update min max time");

  // Return the preorder numbers, the (min, max) for each subtree and the
  // augmented size for each subtree.
  return make_tuple(MM.get_array(), PN.get_array(), aug_sizes.get_array());
}

// Deterministic version
struct BC_BFS_F {
  uintE* Parents;
  BC_BFS_F(uintE* _Parents) : Parents(_Parents) {}
  inline bool update(uintE s, uintE d) {  // Update
    Parents[d] = s;
    return true;
  }
  inline bool updateAtomic(uintE s, uintE d) {  // Atomic version of Update
    return CAS(&Parents[d], UINT_E_MAX, s);
  }
  // Cond function checks if vertex has been visited yet
  inline bool cond(uintE d) { return (Parents[d] == UINT_E_MAX); }
};

template <template <typename W> class vertex, class W, class VS>
uintE* multi_bfs(graph<vertex<W>>& GA, VS& frontier) {
  size_t n = GA.n;
  auto Parents = array_imap<uintE>(n, [](size_t i) { return UINT_E_MAX; });
  frontier.toSparse();
  granular_for(i, 0, frontier.size(), (frontier.size() > 2000), {
    uintE v = frontier.s[i];
    Parents[v] = v;
  });
  while (!frontier.isEmpty()) {
    vertexSubset output =
        edgeMap(GA, frontier, wrap_em_f<W>(BC_BFS_F(Parents.start())), -1,
                sparse_blocked);
    frontier.del();
    frontier = output;
  }
  return Parents.get_array();
}

template <class Seq>
auto cc_sources(Seq& labels) {
  size_t n = labels.size();
  auto flags = array_imap<uintE>(n + 1, [&](size_t i) { return UINT_E_MAX; });
  parallel_for(size_t i = 0; i < n; i++) {
    uintE label = labels[i];
    writeMin(&flags[label], (uintE)i);
  }
  // Get min from each component
  return pbbs::filter(flags, [](uintE v) { return v != UINT_E_MAX; });
}

template <template <class W> class vertex, class W>
auto critical_connectivity(graph<vertex<W>>& GA, uintE* Parents, labels* MM_A,
                           uintE* PN_A, uintE* aug_sizes_A, char* out_f) {
  timer ccc;
  ccc.start();
  size_t n = GA.n;
  auto MM = make_array_imap(MM_A, n);
  auto PN = make_array_imap(PN_A, n);
  auto aug_sizes = make_array_imap(aug_sizes_A, n);

  parallel_for(size_t i = 0; i < n; i++) {
    uintE pi = Parents[i];
    if (pi != i) {
      labels clab = MM[i];
      uintE first_p = PN[pi];
      uintE last_p = first_p + aug_sizes[pi];  // not inclusive
      if ((first_p <= get<0>(clab)) && (get<1>(clab) < last_p)) {  // critical
        Parents[i] |= bc::TOP_BIT;
      }
    }
  }

  auto not_critical_edge = [&](const uintE& u, const uintE& v) {
    uintE e = Parents[u];
    uintE p_u = (e & bc::VAL_MASK);
    if (v == p_u) {
      return !(e & bc::TOP_BIT);
    }
    e = Parents[v];
    uintE p_v = (e & bc::VAL_MASK);
    if (p_v == u) {
      return !(e & bc::TOP_BIT);
    }
    return true;
  };

  timer ccpred;
  ccpred.start();
  // 1. Pack out all critical edges
  auto active = newA(bool, n);
  parallel_for(size_t i = 0; i < n; i++) { active[i] = true; }
  auto vs_active = vertexSubset(n, n, active);
  auto pack_predicate = wrap_f<W>([&](const uintE& src, const uintE& ngh) {
    return not_critical_edge(src, ngh);
  });
  edgeMapFilter(GA, vs_active, pack_predicate, pack_edges | no_output);
  vs_active.del();

  // 2. Run CC on the graph with the critical edges removed to compute
  // a unique label for each biconnected component
  auto cc = cc::CC(GA, 0.2, true);
  ccpred.stop();
  ccpred.reportTotal("cc pred time");

  // Note that counting components here will count initially isolated vertices
  // as distinct components.
  //  auto flags = array_imap<uintE>(n+1, [&] (size_t i) { return 0; });
  //  parallel_for(size_t i=0; i<n; i++) {
  //    if (!flags[cc[i]]) {
  //      flags[cc[i]] = 1;
  //    }
  //  }
  //  pbbs::scan_add(flags, flags);
  //  size_t n_cc = flags[n];
  //  cout << "num biconnected components, including isolated vertices = " <<
  //  flags[n] << endl;

  if (out_f) {
    cout << "Writing labels to file: " << out_f << endl;
    ofstream out(out_f, ofstream::out);
    if (!out.is_open()) {
      std::cout << "Unable to open file " << out_f << std::endl;
      exit(0);
    }

    auto tups = array_imap<pair<uintE, uintE>>(n);
    parallel_for(size_t i=0; i<n; i++) {
      tups[i] = make_pair(Parents[i] & bc::VAL_MASK, cc[i]);
    }

    benchIO::writeArrayToStream(out, tups.start(), n);
//    for (size_t i = 0; i < n; i++) {
//      out << (Parents[i] & bc::VAL_MASK) << " " << cc[i] << endl;
//    }
    out.close();
  }
  cout << "BC done" << endl;
  free(MM_A);
  free(PN_A);
  free(aug_sizes_A);
  ccc.stop();
  ccc.reportTotal("critical conn time");
  return make_tuple(Parents, cc);
}

// CC -> BFS from one source from each component = set of BFS trees in a single
// array
template <class vertex>
auto Biconnectivity(graph<vertex>& GA, char* out_f = 0) {
  size_t n = GA.n;

  timer fcc;
  fcc.start();
  uintE* Components = cc::CC(GA, 0.2, false);
  fcc.stop();
  fcc.reportTotal("first cc");
  auto c_im = make_in_imap<uintE>(n, [&](size_t i) { return Components[i]; });

  timer sc;
  sc.start();
  auto Sources = cc_sources(c_im);
  free(Components);

  auto Sources_copy = Sources.copy();
  auto Centers = vertexSubset(n, Sources.size(), Sources.get_array());
  auto Parents = multi_bfs(GA, Centers);
  sc.stop();
  sc.reportTotal("sc, multibfs time");

  // Returns ((min, max), preorder#, and augmented sizes) of each subtree.
  timer pn;
  pn.start();

  labels* min_max;
  uintE* preorder_num;
  uintE* aug_sizes;
  std::tie(min_max, preorder_num, aug_sizes) =
      preorder_number(GA, Parents, Sources_copy);
  pn.stop();
  pn.reportTotal("preorder time");

  critical_connectivity(GA, Parents, min_max, preorder_num, aug_sizes, out_f);
}
