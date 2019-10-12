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

#include "ligra/bridge.h"
#include "ligra/ligra.h"
#include "pbbslib/random.h"
#include "union_find_rules.h"

#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"

#include <iostream>
#include <limits.h>
#include <vector>
#include <mutex>
#include <atomic>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <random>


namespace union_find {

/* ================================== CSR templates ================================== */

/* Used if the algorithm only requires a parents array */
template <class Find, class Unite, class G>
struct UnionFindTemplate {
  G& GA;
  Unite& unite;
  Find& find;
  bool use_hooks;
  UnionFindTemplate(G& GA, Unite& unite, Find& find, uint32_t neighbor_rounds = 2, bool use_hooks=false) : GA(GA), unite(unite), find(find), use_hooks(use_hooks) {}

  pbbs::sequence<uintE> components() {
    using W = typename G::weight_type;
    size_t n = GA.n;

    auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
    pbbs::sequence<uintE> hooks;
    if (use_hooks) {
      hooks = pbbs::sequence<uintE>(n, UINT_E_MAX);
    }

    timer ut; ut.start();
    parallel_for(0, n, [&] (size_t i) {
      auto map_f = [&] (uintE u, uintE v, const W& wgh) {
        if (u < v) {
          if (use_hooks) {
            unite(u, v, parents, hooks);
          } else {
            unite(u, v, parents);
          }
        }
      };
      GA.get_vertex(i).mapOutNgh(i, map_f, false); // in parallel
    }, 1);
    ut.stop(); debug(ut.reportTotal("union time"));

    timer ft; ft.start();
    parallel_for(0, n, [&] (size_t i) {
      parents[i] = find(i,parents);
    });
    ft.stop(); debug(ft.reportTotal("find time"););
    return parents;
  }
};

/* ****************************** Sampling ******************************/

// From gapbs/cc.c, Thanks to S. Beamer + M. Sutton for the well documented
// reference implementation of afforest.
template <class Seq>
std::pair<typename Seq::value_type, double>
sample_frequent_element(Seq& S, uint64_t num_samples=1024) {
  using T = typename Seq::value_type;
  std::unordered_map<T, int> sample_counts(32);
  using kvp_type = typename std::unordered_map<T, int>::value_type;
  // Sample elements from 'S'
  std::mt19937 gen;
  std::uniform_int_distribution<T> distribution(0, S.size() - 1);
  for (T i = 0; i < num_samples; i++) {
    T n = distribution(gen);
    sample_counts[S[n]]++;
  }
  // Find most frequent element in samples (estimate of most frequent overall)
  auto most_frequent = std::max_element(
    sample_counts.begin(), sample_counts.end(),
    [](const kvp_type& a, const kvp_type& b) { return a.second < b.second; });

  double frac_of_graph = static_cast<double>(most_frequent->second) / num_samples;
  std::cout
    << "Skipping largest intermediate component (ID: " << most_frequent->first
    << ", approx. " << (frac_of_graph * 100)
    << "% of the graph)" << std::endl;
  return std::make_pair(most_frequent->first, frac_of_graph);
}

// returns a component labeling
// Based on the implementation in gapbs/cc.c, Thanks to S. Beamer + M. Sutton
// for the well documented reference implementation of afforest.
template <class Find, class Unite, class G>
struct UnionFindSampleTemplate {
  G& GA;
  Find& find;
  Unite& unite;
  uint32_t neighbor_rounds;
  bool use_hooks;

  UnionFindSampleTemplate(G& GA, Unite& unite, Find& find,
      uint32_t neighbor_rounds = 2, bool use_hooks=false) :
   GA(GA), unite(unite), find(find), neighbor_rounds(neighbor_rounds), use_hooks(use_hooks) {}

  pbbs::sequence<uintE> components() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    cout << "neighbor_rounds = " << neighbor_rounds << endl;

    auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
    pbbs::sequence<uintE> hooks;
    if (use_hooks) {
      hooks = pbbs::sequence<uintE>(n, UINT_E_MAX);
    }

    pbbs::random rnd;
    timer st; st.start();

    // Using random neighbor---some overhead (and faster for some graphs), also
    // theoretically defensible
    for (uint32_t r=0; r<neighbor_rounds; r++) {
      parallel_for(0, n, [&] (size_t u) {
        auto u_rnd = rnd.fork(u);
        auto u_vtx = GA.get_vertex(u);
        if (u_vtx.getOutDegree() > 0) {
          size_t ngh_idx = u_rnd.rand() % u_vtx.getOutDegree();
          uintE ngh; W wgh;
          std::tie(ngh, wgh) = u_vtx.get_ith_out_neighbor(u, ngh_idx);
          if (use_hooks) {
            unite(u, ngh, parents, hooks);
          } else {
            unite(u, ngh, parents);
          }
        }
      }, 512);
      // compress nodes fully (turns out this is faster)
      parallel_for(0, n, [&] (size_t u) {
        parents[u] = find(u, parents);
      }, 512);
      rnd = rnd.next();
    }

//    // Using r'th neighbor---usually faster, but not so defensible theoretically.
//    for (uint32_t r=0; r<neighbor_rounds; r++) {
//      if (r == 0) {
//        /* First round: sample a directed forest and compress instead of using
//         * unite */
//        parallel_for(0, n, [&] (size_t u) {
//          auto u_vtx = GA.get_vertex(u);
//          if (u_vtx.getOutDegree() > r) {
//            uintE ngh; W wgh;
//            std::tie(ngh, wgh) = u_vtx.get_ith_out_neighbor(u, r);
//            if (u_vtx.getOutDegree() < GA.get_vertex(ngh).getOutDegree()) {
//              parents[u] = ngh;
//            }
//          }
//        }, 512);
//      } else {
//        /* Subsequent rounds: use unite */
//        parallel_for(0, n, [&] (size_t u) {
//          auto u_vtx = GA.get_vertex(u);
//          if (u_vtx.getOutDegree() > r) {
//            uintE ngh; W wgh;
//            std::tie(ngh, wgh) = u_vtx.get_ith_out_neighbor(u, r);
//            if (use_hooks) {
//              unite(u, ngh, parents, hooks);
//            } else {
//              unite(u, ngh, parents);
//            }
//          }
//        }, 512);
//      }
//      // compress nodes fully (turns out this is faster)
//      parallel_for(0, n, [&] (size_t u) {
//        parents[u] = find(u, parents);
//      }, 512);
//    }


    uintE frequent_comp; double pct;
    std::tie(frequent_comp, pct) = sample_frequent_element(parents);
    st.stop(); st.reportTotal("sample time");

    timer ut; ut.start();
    parallel_for(0, n, [&] (size_t u) {
      // Only process edges for vertices not linked to the main component
      // note that this is safe only for undirected graphs. For directed graphs,
      // the in-edges must also be explored for all vertices.
      if (parents[u] != frequent_comp) {
        auto map_f = [&] (uintE u, uintE v, const W& wgh) {
          if (u < v) {
            if (use_hooks) {
              unite(u, v, parents, hooks);
            } else {
              unite(u, v, parents);
            }
          }
        };
        GA.get_vertex(u).mapOutNgh(u, map_f); // in parallel
      }
    }, 1);
    ut.stop(); ut.reportTotal("union time");

    timer ft; ft.start();
    parallel_for(0, n, [&] (size_t i) {
      parents[i] = find(i,parents);
    });
    ft.stop(); ft.reportTotal("find time");
    return parents;
   }
};


template <class W>
struct BFS_ComponentLabel_F {
  uintE* Parents;
  uintE src;
  BFS_ComponentLabel_F(uintE* _Parents, uintE src) : Parents(_Parents), src(src) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (Parents[d] == UINT_E_MAX) {
      Parents[d] = src;
      return 1;
    } else {
      return 0;
    }
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    return (pbbslib::atomic_compare_and_swap(&Parents[d], UINT_E_MAX, src));
  }
  inline bool cond(const uintE& d) { return (Parents[d] == UINT_E_MAX); }
};

template <class Graph>
inline sequence<uintE> BFS_ComponentLabel(Graph& G, uintE src) {
  using W = typename Graph::weight_type;
  /* Creates Parents array, initialized to all -1, except for src. */
  auto Parents = sequence<uintE>(G.n, [&](size_t i) { return UINT_E_MAX; });
  Parents[src] = src;

  vertexSubset Frontier(G.n, src);
  size_t reachable = 0; size_t rounds = 0;
  while (!Frontier.isEmpty()) {
    reachable += Frontier.size();
    vertexSubset output =
        edgeMap(G, Frontier, BFS_ComponentLabel_F<W>(Parents.begin(), src), -1, sparse_blocked | dense_parallel);
    Frontier.del();
    Frontier = output;
    rounds++;
  }
  Frontier.del();
  std::cout << "Reachable: " << reachable << " #rounds = " << rounds << std::endl;
  return Parents;
}

template <class Find, class Unite, class G>
struct UnionFindSampledBFSTemplate {
  G& GA;
  Find& find;
  Unite& unite;
  uint32_t neighbor_rounds;
  bool use_hooks;

  UnionFindSampledBFSTemplate(G& GA, Unite& unite, Find& find,
      uint32_t neighbor_rounds = 2, bool use_hooks=false) :
   GA(GA), unite(unite), find(find), neighbor_rounds(neighbor_rounds), use_hooks(use_hooks) {}

  pbbs::sequence<uintE> components() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    cout << "neighbor_rounds = " << neighbor_rounds << endl;

    pbbs::sequence<uintE> parents;
    pbbs::sequence<uintE> hooks;
    if (use_hooks) {
      hooks = pbbs::sequence<uintE>(n, UINT_E_MAX);
    }

    pbbs::random rnd;
    timer st; st.start();

    /* Assume that G has a massive component with size at least 10% of the
     * graph. Run BFSs at most max_trials times until we find it. */

    uint32_t max_trials = 3;
    uintE skip_comp = UINT_E_MAX;
    bool found_massive_component = false;
    for (uint32_t r=0; r<max_trials; r++) {
      uintE src = rnd.rand() % n;
      auto bfs_parents = BFS_ComponentLabel(GA, src);

      uintE frequent_comp; double pct;
      parents = std::move(bfs_parents);
      std::tie(frequent_comp, pct) = sample_frequent_element(parents);
      if (pct > static_cast<double>(0.1)) {
        std::cout << "BFS covered: " << pct << " of graph" << std::endl;
        skip_comp = frequent_comp;
        found_massive_component = true;
        break;
      }
      std::cout << "BFS covered only: " << pct << " of graph." << std::endl;
      rnd = rnd.next();
    }

    parallel_for(0, n, [&] (size_t i) {
      if (!found_massive_component || parents[i] != skip_comp) {
        parents[i] = i;
      }
    });

    auto bfs_parents = parents; // copy

    st.stop(); st.reportTotal("sample time");

    timer ut; ut.start();
    parallel_for(0, n, [&] (size_t u) {
      // Only process edges for vertices not linked to the main component
      // note that this is safe only for undirected graphs. For directed graphs,
      // the in-edges must also be explored for all vertices.
      if (bfs_parents[u] != skip_comp) {
        auto map_f = [&] (uintE u, uintE v, const W& wgh) {
          if (u < v) {
            if (use_hooks) {
              unite(u, v, parents, hooks);
            } else {
              unite(u, v, parents);
            }
          }
        };
        GA.get_vertex(u).mapOutNgh(u, map_f); // in parallel
      }
    }, 1024);
    ut.stop(); ut.reportTotal("union time");

    timer ft; ft.start();
    parallel_for(0, n, [&] (size_t i) {
      parents[i] = find(i,parents);
    });
    ft.stop(); ft.reportTotal("find time");
    return parents;
   }
};

template <class Find, class Unite, class G>
struct UnionFindLDDTemplate {
  G& GA;
  Find& find;
  Unite& unite;
  uint32_t neighbor_rounds;
  bool use_hooks;

  UnionFindLDDTemplate(G& GA, Unite& unite, Find& find,
      uint32_t neighbor_rounds = 2, bool use_hooks=false) :
   GA(GA), unite(unite), find(find), neighbor_rounds(neighbor_rounds), use_hooks(use_hooks) { }

  pbbs::sequence<uintE> components() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    cout << "neighbor_rounds = " << neighbor_rounds << endl;

    pbbs::sequence<uintE> hooks;
    if (use_hooks) {
      hooks = pbbs::sequence<uintE>(n, UINT_E_MAX);
    }

    pbbs::random rnd;
    timer st; st.start();

    auto clusters = LDD(GA, 0.2, /* permute = */false);

   parallel_for(0, n, [&] (uintE u) {
     if (clusters[u] == u) { // root, ok
     } else {
       assert(clusters[clusters[u]] == clusters[u]);
     }
   });

    pbbs::sequence<uintE> parents(n);
    parallel_for(0, n, [&] (size_t i) {
      parents[i] = clusters[i];
    });

    uintE frequent_comp; double pct;
    std::tie(frequent_comp, pct) = sample_frequent_element(parents);
    st.stop(); st.reportTotal("sample time");

    timer ut; ut.start();
    parallel_for(0, n, [&] (size_t u) {
      // Only process edges for vertices not linked to the main component
      // note that this is safe only for undirected graphs. For directed graphs,
      // the in-edges must also be explored for all vertices.
      if (clusters[u] != frequent_comp) {
        auto map_f = [&] (uintE u, uintE v, const W& wgh) {
          if (u < v) {
            if (use_hooks) {
              unite(u, v, parents, hooks);
            } else {
              unite(u, v, parents);
            }
          }
        };
        GA.get_vertex(u).mapOutNgh(u, map_f); // in parallel
      }
    }, 512);
    ut.stop(); ut.reportTotal("union time");

    timer ft; ft.start();
    parallel_for(0, n, [&] (size_t i) {
      parents[i] = find(i,parents);
    });
    ft.stop(); ft.reportTotal("find time");
    return parents;
   }
};



/* ================================== COO templates ================================== */

template <class Find, class Unite, class W>
inline pbbs::sequence<uintE> UnionFindTemplate_coo(edge_array<W>& G, Unite& unite, Find& find) {
  size_t n = G.num_rows;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  pbbs::sequence<uintE> hooks;
  /* TODO */

  timer ut; ut.start();
  parallel_for(0, G.non_zeros, [&] (size_t i) {
    uintE u, v; W wgh;
    std::tie(u, v, wgh) = G.E[i];
    if (u < v) {
      unite(u, v, parents);
    }
  }, 512);
  ut.stop(); debug(ut.reportTotal("union time"));

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = find(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  return parents;
}


//template <template <class Find> class Unite,
//          template <class F, class U, class G> class UFTemplate,
//          class Graph>
//pbbs::sequence<uintE> select_framework_algorithm(Graph& G, std::string& find_arg, uint32_t sampling_rounds=2, bool use_hooks=false) {
//  if (find_arg == "find_compress") {
//    auto find = find_variants::find_compress;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_naive") {
//    auto find = find_variants::find_naive;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_split") {
//    auto find = find_variants::find_split;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_halve") {
//    auto find = find_variants::find_halve;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_atomic_split") {
//    auto find = find_variants::find_atomic_split;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_atomic_halve") {
//    auto find = find_variants::find_atomic_halve;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  }
//  return pbbs::sequence<uintE>();
//}


//template <template <class Find, class Compress> class Unite,
//          template <class F, class U, class G> class UFTemplate,
//          class Find,
//          class Graph>
//pbbs::sequence<uintE> select_framework_algorithm(Graph& G, Find& find, std::string& compress_arg, uint32_t sampling_rounds=2, bool use_hooks=false) {
//  if (compress_arg == "naive") {
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (compress_arg == "find_naive") {
//    auto find = find_variants::find_naive;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  }
//  return pbbs::sequence<uintE>();
//}
//
//
//
//template <template <class Find, class Compress> class Unite,
//          template <class F, class U, class G> class UFTemplate,
//          class Graph>
//pbbs::sequence<uintE> select_framework_algorithm(Graph& G, std::string& find_arg, uint32_t sampling_rounds=2, bool use_hooks=false) {
//  if (find_arg == "find_compress") {
//    auto find = find_variants::find_compress;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_naive") {
//    auto find = find_variants::find_naive;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_split") {
//    auto find = find_variants::find_split;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_halve") {
//    auto find = find_variants::find_halve;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_atomic_split") {
//    auto find = find_variants::find_atomic_split;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  } else if (find_arg == "find_atomic_halve") {
//    auto find = find_variants::find_atomic_halve;
//    auto unite = Unite<decltype(find)>(G.n, find);
//    auto q = UFTemplate<decltype(find), decltype(unite), Graph>(G, unite, find, sampling_rounds, use_hooks);
//    return q.components();
//  }
//  return pbbs::sequence<uintE>();
//}

}  // namespace union_find
