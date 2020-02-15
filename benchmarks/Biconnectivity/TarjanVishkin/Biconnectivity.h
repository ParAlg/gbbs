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

#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"
#include "ligra/ligra.h"
#include "ligra/pbbslib/dyn_arr.h"
#include "ligra/pbbslib/sparse_table.h"
#include "pbbslib/random.h"
#include "pbbslib/sample_sort.h"
#include "ligra/speculative_for.h"

namespace bc {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;
};  // namespace bc

using labels = std::tuple<uintE, uintE>;

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
    pbbslib::write_add(&sizes[d], sizes[s]);
    intE new_value = pbbslib::fetch_and_add(&cts[d], -1) - 1;
    return (new_value == 0);
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
    uintE low = std::get<0>(lab_s), high = std::get<1>(lab_s);
    if (low < std::get<0>(MM[d])) {
      std::get<0>(MM[d]) = low;
    }
    if (high > std::get<1>(MM[d])) {
      std::get<1>(MM[d]) = high;
    }
    uintE ct = cts[d] - 1;
    cts[d] = ct;
    return ct == 0;
  }
  bool updateAtomic(const uintE& s, const uintE& d) {
    labels lab_s = MM[s];
    uintE low = std::get<0>(lab_s), high = std::get<1>(lab_s);
    if (low < std::get<0>(MM[d])) {
      pbbslib::write_min(&std::get<0>(MM[d]), low);
    }
    if (high > std::get<1>(MM[d])) {
      pbbslib::write_max(&std::get<1>(MM[d]), high);
    }
    intE new_value = pbbslib::fetch_and_add(&cts[d], -1) - 1;
    return (new_value == 0);
  }
  bool cond(const uintE& d) { return true; }
};

template <template <typename W> class vertex, class W, class Seq>
inline std::tuple<labels*, uintE*, uintE*> preorder_number(symmetric_graph<vertex, W>& GA,
                                                           uintE* Parents,
                                                           Seq& Sources) {
  size_t n = GA.n;
  using edge = std::tuple<uintE, uintE>;
  auto out_edges = sequence<edge>(
      n, [](size_t i) { return std::make_tuple(UINT_E_MAX, UINT_E_MAX); });
  par_for(0, n, [&] (size_t i) {
    uintE p_i = Parents[i];
    if (p_i != i) {
      out_edges[i] = std::make_tuple(p_i, i);
    }
  });

  auto edges = pbbslib::filter(
      out_edges, [](const edge& e) { return std::get<0>(e) != UINT_E_MAX; });
  out_edges.clear();
  auto sort_tup = [](const edge& l, const edge& r) { return l < r; };
  pbbslib::sample_sort_inplace(edges.slice(), sort_tup, true);

  auto starts = sequence<uintE>(n + 1, [](size_t i) { return UINT_E_MAX; });
  par_for(0, edges.size(), [&] (size_t i) {
    if (i == 0 || std::get<0>(edges[i]) != std::get<0>(edges[i - 1])) {
      starts[std::get<0>(edges[i])] = i;
    }
  });
  starts[n] = edges.size();

  timer seq;
  seq.start();
  for (long i = starts.size() - 1; i >= 0; i--) {
    if (starts[i] == UINT_E_MAX) {
      starts[i] = starts[i + 1];
    }
  }
  seq.stop();
  debug(seq.reportTotal("seq time"););

  auto nghs = sequence<uintE>(
      edges.size(), [&](size_t i) { return std::get<1>(edges[i]); });
  edges.clear();

  timer augs;
  augs.start();

  // Create directed BFS tree

  auto v_out = pbbslib::new_array_no_init<vertex_data>(n);
  auto v_in = pbbslib::new_array_no_init<vertex_data>(n);
  par_for(0, n, [&] (size_t i) {
    uintE out_off = starts[i];
    uintE out_deg = starts[i + 1] - out_off;
    v_out[i].offset = out_off;
    v_out[i].degree = out_deg;

    uintE parent = Parents[i];
    if (parent != i) {
      v_in[i].degree = 1;
      v_in[i].offset = i;
    } else {
      v_in[i].degree = 0;
      v_in[i].offset = i;
    }
  });
  auto Tree = asymmetric_graph<asymmetric_vertex, pbbslib::empty>(v_out, v_in, n, nghs.size(), [](){}, (std::tuple<uintE, pbbs::empty>*)nghs.begin(), (std::tuple<uintE, pbbs::empty>*)Parents);

  // 1. Leaffix for Augmented Sizes
  auto aug_sizes = sequence<uintE>(n, [](size_t i) { return 1; });
  auto cts =
      sequence<intE>(n, [&](size_t i) { return Tree.get_vertex(i).getOutDegree(); });
  auto leaf_im = pbbslib::make_sequence<bool>(n, [&](size_t i) {
    auto s_i = starts[i];
    auto s_n = starts[i + 1];
    size_t deg = s_n - s_i;
    return (Parents[i] != i) && deg == 0;
  });
  auto leafs = pbbslib::pack_index<uintE>(leaf_im);

  auto vs = vertexSubset(n, leafs.size(), leafs.begin());
  size_t rds = 0, tv = 0;
  while (!vs.isEmpty()) {
    rds++;
    tv += vs.size();
    // histogram or write-add parents, produce next em.
    auto output = edgeMap(
        Tree, vs, wrap_em_f<pbbslib::empty>(AugF(aug_sizes.begin(), cts.begin())),
        -1, in_edges | sparse_blocked | fine_parallel);
    if (rds > 1) {
      vs.del();  // don't delete leafs.
    }
    vs = output;
    output.toSparse();
  }
  augs.stop();
  debug(augs.reportTotal("aug size time"););

  // Optional: Prefix sum over sources with aug-size (if we want distinct
  // preorder #s)

  // Use copy constructor
  sequence<uintE> s_copy(Sources);
  size_t sources_size = Sources.size();

  timer pren;
  pren.start();
  auto PN = sequence<uintE>(n);
  vs = vertexSubset(n, sources_size, s_copy.to_array());
  par_for(0, Sources.size(), [&] (size_t i) {
    uintE v = vs.vtx(i);
    PN[v] = 0;
  });
  rds = 0;
  tv = 0;
  while (!vs.isEmpty()) {
    rds++;
    tv += vs.size();
    auto offsets = sequence<uintE>(vs.size(), [&](size_t i) {
      uintE v = vs.s[i];
      return Tree.get_vertex(v).getOutDegree();
    });
    auto tot = pbbslib::scan_add_inplace(offsets);
    auto next_vs = sequence<uintE>(tot);
    par_for(0, vs.size(), 1, [&] (size_t i) {
      uintE v = vs.s[i];
      uintE off = offsets[i];
      uintE deg_v = Tree.get_vertex(v).getOutDegree();
      uintE preorder_number = PN[v] + 1;

      // should be tuned
      if (deg_v < 4000) {
        // Min and max in any vertex are [PN[v], PN[v] + aug_sizes[v])
        for (size_t j = 0; j < deg_v; j++) {
          uintE ngh = Tree.get_vertex(v).getOutNeighbor(j);
          PN[ngh] = preorder_number;
          preorder_number += aug_sizes[ngh];
          next_vs[off + j] = ngh;
        }
      } else {
        auto A = sequence<uintE>(deg_v);
        par_for(0, deg_v, [&] (size_t j) {
          uintE ngh = Tree.get_vertex(v).getOutNeighbor(j);
          A[j] = aug_sizes[ngh];
        });
        pbbslib::scan_add_inplace(A);
        par_for(0, deg_v, [&] (size_t j) {
          uintE ngh = Tree.get_vertex(v).getOutNeighbor(j);
          uintE pn = preorder_number + A[j];
          PN[ngh] = pn;
          next_vs[off + j] = ngh;
        });
      }
    });
    vs.del();
    size_t next_vs_size = next_vs.size();
    vs = vertexSubset(n, next_vs_size, next_vs.to_array());
  }
  pren.stop();
  debug(pren.reportTotal("preorder number from sizes time"););

  timer map_e;
  map_e.start();
  // Map all edges and update labels.
  auto MM = sequence<labels>(n, [&](size_t i) {
    uintE pn_i = PN[i];
    return std::make_tuple(pn_i, pn_i + aug_sizes[i] - 1);
  });
  // Note that we have to exclude tree edges: it's all edges spliced in with
  // the exception of tree edges.
  auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
    if (u < v) {  // see if assigning roles helps temp loc later.
      if (u == Parents[v] || v == Parents[u]) {
        return;
      }
      uintE p_u = PN[u];
      uintE p_v = PN[v];
      if (p_u < std::get<0>(MM[v])) {
        pbbslib::write_min(&std::get<0>(MM[v]), p_u);
      } else if (p_u > std::get<1>(MM[v])) {
        pbbslib::write_max(&std::get<1>(MM[v]), p_u);
      }
      if (p_v < std::get<0>(MM[u])) {
        pbbslib::write_min(&std::get<0>(MM[u]), p_v);
      } else if (p_v > std::get<1>(MM[u])) {
        pbbslib::write_max(&std::get<1>(MM[u]), p_v);
      }
    }
  };
  par_for(0, n, 1, [&] (size_t i) { GA.get_vertex(i).mapOutNgh(i, map_f); });
  map_e.stop();
  debug(map_e.reportTotal("map edges time"););

  timer leaff;
  leaff.start();
  // 1. Leaffix to update min/max
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { cts[i] = Tree.get_vertex(i).getOutDegree(); });

  size_t leafs_size = leafs.size();
  vs = vertexSubset(n, leafs_size, leafs.to_array());
  rds = 0, tv = 0;
  while (!vs.isEmpty()) {
    rds++;
    tv += vs.size();
    // histogram or write-add parents, produce next em.
    auto output = edgeMap(
        Tree, vs, wrap_em_f<pbbslib::empty>(MinMaxF(MM.begin(), cts.begin())), -1,
        in_edges | sparse_blocked | fine_parallel);
    vs.del();
    vs = output;
  }
  // Delete tree
  pbbslib::free_array(v_out);
  pbbslib::free_array(v_in);
  nghs.clear();
  leaff.stop();
  debug(leaff.reportTotal("leaffix to update min max time"););

  // Return the preorder numbers, the (min, max) for each subtree and the
  // augmented size for each subtree.
  return std::make_tuple(MM.to_array(), PN.to_array(), aug_sizes.to_array());
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
    return pbbslib::CAS(&Parents[d], UINT_E_MAX, s);
  }
  // Cond function checks if vertex has been visited yet
  inline bool cond(uintE d) { return (Parents[d] == UINT_E_MAX); }
};

template <template <typename W> class vertex, class W, class VS>
inline uintE* multi_bfs(symmetric_graph<vertex, W>& GA, VS& frontier) {
  size_t n = GA.n;
  auto Parents = sequence<uintE>(n, [](size_t i) { return UINT_E_MAX; });
  frontier.toSparse();
  par_for(0, frontier.size(), 2000, [&] (size_t i) {
    uintE v = frontier.s[i];
    Parents[v] = v;
  });
  while (!frontier.isEmpty()) {
    vertexSubset output =
        edgeMap(GA, frontier, wrap_em_f<W>(BC_BFS_F(Parents.begin())), -1,
                sparse_blocked);
    frontier.del();
    frontier = output;
  }
  return Parents.to_array();
}

// Deterministic version
struct DET_BFS_F {
  bool* visited;
  uintE* Parents;
  DET_BFS_F(bool* visited, uintE* _Parents) : visited(visited), Parents(_Parents) {}
  inline bool update(uintE s, uintE d) {  // Update
    if (s < Parents[d]) {
      Parents[d] = s;
    }
    return false;
  }
  inline bool updateAtomic(uintE s, uintE d) {  // Atomic version of Update
    pbbslib::write_min(&Parents[d], s);
    return false;
  }
  // Cond function checks if vertex has been visited yet
  inline bool cond(uintE d) { return !visited[d]; }
};


// Deterministic version
struct DET_BFS_F_2 {
  bool* visited;
  uintE* Parents;
  DET_BFS_F_2(bool* visited, uintE* _Parents) : visited(visited), Parents(_Parents) {}
  inline bool update(uintE s, uintE d) {  // Update
    if (Parents[d] == s) {
      visited[d] = true;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uintE s, uintE d) {  // Atomic version of Update
    if (Parents[d] == s) {
      visited[d] = true;
      return true;
    }
    return false;
  }
  // Cond function checks if vertex has been visited yet
  inline bool cond(uintE d) { return !visited[d]; }
};

template <template <typename W> class vertex, class W, class VS>
uintE* deterministic_multi_bfs(symmetric_graph<vertex, W>& GA, VS& frontier) {
  size_t n = GA.n;
  auto visited = sequence<bool>(n, [] (size_t i) { return false; });
  auto Parents = sequence<uintE>(n, [](size_t i) { return UINT_E_MAX; });
  frontier.toSparse();
  par_for(0, frontier.size(), [&] (size_t i) {
    uintE v = frontier.s[i];
    Parents[v] = v;
    visited[v] = true;
  });
  while (!frontier.isEmpty()) {
    edgeMap(GA, frontier, wrap_em_f<W>(DET_BFS_F(visited.begin(), Parents.begin())), -1,
                sparse_blocked);
    vertexSubset output =
        edgeMap(GA, frontier, wrap_em_f<W>(DET_BFS_F_2(visited.begin(), Parents.begin())), -1,
                sparse_blocked);

    frontier.del();
    frontier = output;
  }
  return Parents.to_array();
}


template <class Seq>
inline sequence<uintE> cc_sources(Seq& labels) {
  size_t n = labels.size();
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return UINT_E_MAX; });
  par_for(0, n, [&] (size_t i) {
    uintE label = labels[i];
    pbbslib::write_min(&flags[label], (uintE)i);
  });
  // Get min from each component
  return pbbslib::filter(flags, [](uintE v) { return v != UINT_E_MAX; });
}

template <template <class W> class vertex, class W>
inline std::tuple<uintE*, uintE*> critical_connectivity(
    symmetric_graph<vertex, W>& GA, uintE* Parents, labels* MM_A, uintE* PN_A,
    uintE* aug_sizes_A, char* out_f) {
  timer ccc;
  ccc.start();
  size_t n = GA.n;
  auto MM = pbbslib::make_sequence<labels>(MM_A, n);
  auto PN = pbbslib::make_sequence<uintE>(PN_A, n);
  auto aug_sizes = pbbslib::make_sequence<uintE>(aug_sizes_A, n);

  par_for(0, n, [&] (size_t i) {
    uintE pi = Parents[i];
    if (pi != i) {
      labels clab = MM[i];
      uintE first_p = PN[pi];
      uintE last_p = first_p + aug_sizes[pi];  // not inclusive
      if ((first_p <= std::get<0>(clab)) &&
          (std::get<1>(clab) < last_p)) {  // critical
        Parents[i] |= bc::TOP_BIT;
      }
    }
  });

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
  auto active = pbbslib::new_array_no_init<bool>(n);
  par_for(0, n, [&] (size_t i) { active[i] = true; });
//  auto vs_active = vertexSubset(n, n, active);
  auto pack_predicate = [&](const uintE& src, const uintE& ngh, const W& wgh) -> int {
    return !not_critical_edge(src, ngh);
  };
  timer ft; ft.start();
  filter_edges(GA, pack_predicate);
//  edgeMapFilter(GA, vs_active, pack_predicate, pack_edges | no_output);
  ft.stop(); debug(ft.reportTotal("filter edges time"););
//  vs_active.del();

  // 2. Run CC on the graph with the critical edges removed to compute
  // a unique label for each biconnected component
  auto cc = workefficient_cc::CC(GA, 0.2, true);
  ccpred.stop();
  debug(ccpred.reportTotal("cc pred time"););

//  //Note that counting components here will count initially isolated vertices
//  //as distinct components.
//   auto flags = sequence<uintE>(n+1, [&] (size_t i) { return 0; });
//   par_for(0, n, [&] (size_t i) {
//     if (!flags[cc[i]]) {
//       flags[cc[i]] = 1;
//     }
//   });
//   flags[n] = 0;
//   pbbslib::scan_add_inplace(flags.slice());
//   size_t n_cc = flags[n];
//   std::cout << "num biconnected components, including isolated vertices = "
//   << flags[n] << "\n";

  if (out_f) {
    std::cout << "Writing labels to file: " << out_f << "\n";
//    std::ofstream out(out_f, std::ofstream::out);
//    if (!out.is_open()) {
//      std::cout << "Unable to open file " << out_f << "\n";
//      exit(0);
//    }

    auto tups = sequence<std::pair<uintE, uintE>>(n);
    par_for(0, n, [&] (size_t i) {
        tups[i] = std::make_pair(Parents[i] & bc::VAL_MASK, cc[i]); });

    auto C = pbbslib::sequence_to_string(tups);
    pbbs::char_seq_to_file(C, out_f);
    // benchIO::writeArrayToStream(out, tups.begin(), n);
    // for (size_t i = 0; i < n; i++) {
    //   out << (Parents[i] & bc::VAL_MASK) << " " << cc[i] << "\n";
    // }
    // out.close();
  }
  debug(std::cout << "Bicc done"
            << "\n";);
  pbbslib::free_array(MM_A);
  pbbslib::free_array(PN_A);
  pbbslib::free_array(aug_sizes_A);
  ccc.stop();
  debug(ccc.reportTotal("critical conn time"););
  return std::make_tuple(Parents, cc.to_array());
}

// CC -> BFS from one source from each component = set of BFS trees in a single
// array
template <template <class W> class vertex, class W>
inline std::tuple<uintE*, uintE*> Biconnectivity(symmetric_graph<vertex, W>& GA,
                                                 char* out_f = 0) {
  size_t n = GA.n;

  timer fcc;
  fcc.start();
  sequence<uintE> Components = workefficient_cc::CC(GA, 0.2, false);
  fcc.stop();
  debug(fcc.reportTotal("first cc"););

  timer sc;
  sc.start();
  auto Sources = cc_sources(Components);
  Components.clear();

//  auto Sources_copy = Sources.copy(Sources);
  auto Sources_copy = Sources; // Use copy constructor
  auto Centers = vertexSubset(n, Sources_copy.size(), Sources.to_array());
//  auto Parents = deterministic_multi_bfs(GA, Centers); // useful for debugging
  auto Parents = multi_bfs(GA, Centers);
  sc.stop();
  debug(sc.reportTotal("sc, multibfs time"););

  // Returns ((min, max), preorder#, and augmented sizes) of each subtree.
  timer pn;
  pn.start();

  labels* min_max;
  uintE* preorder_num;
  uintE* aug_sizes;
  std::tie(min_max, preorder_num, aug_sizes) =
      preorder_number(GA, Parents, Sources_copy);
  pn.stop();
  debug(pn.reportTotal("preorder time"););

  return critical_connectivity(GA, Parents, min_max, preorder_num, aug_sizes,
                               out_f);
}
