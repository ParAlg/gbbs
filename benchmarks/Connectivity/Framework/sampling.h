#pragma once

#include "utils.h"

/* ****************************** Sampling ******************************/

template <
    class Graph,
    class Sampler,
    class Algorithm,
    AlgorithmType algorithm_type,
    SamplingOption sampling_option>
  struct SamplingAlgorithmTemplate {
    Graph& G;
    Sampler& sampler;
    Algorithm& algorithm;
    SamplingAlgorithmTemplate(Graph& G, Sampler& sampler, Algorithm& algorithm) : G(G), sampler(sampler), algorithm(algorithm) {}

    pbbs::sequence<parent> components() {
      auto parents = sampler.initial_components();

      parent frequent_comp; double pct;
      std::tie(frequent_comp, pct) = sample_frequent_element(parents);

      algorithm.initialize(parents);

      if constexpr (algorithm_type == liu_tarjan_type) {
        parallel_for(0, G.n, [&] (size_t i) {
          if (parents[i] == frequent_comp) {
            parents[i] = largest_comp;
          }
        }, 2048);
        algorithm.template compute_components<sampling_option>(parents, largest_comp);
        parallel_for(0, G.n, [&] (size_t i) {
          if (parents[i] == largest_comp) {
            parents[i] = frequent_comp;
          }
        }, 2048);
      } else {
        algorithm.template compute_components<sampling_option>(parents, frequent_comp);
      }
      return parents;
    }
  };

template <
    class Graph,
    class Algorithm,
    AlgorithmType algorithm_type>
  struct NoSamplingAlgorithmTemplate {
    Graph& G;
    Algorithm& algorithm;
    NoSamplingAlgorithmTemplate(Graph& G, Algorithm& algorithm) : G(G), algorithm(algorithm) {}

    pbbs::sequence<parent> components() {
      size_t n = G.n;
      auto parents = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
      algorithm.initialize(parents);
      algorithm.template compute_components<no_sampling>(parents);
      return parents;
    }
  };

/* *********************** Specific sampling strategies ********************* */

// returns a component labeling
// Based on the implementation in gapbs/cc.c, Thanks to S. Beamer + M. Sutton
// for the well documented reference implementation of afforest.
template <class Find, class Unite, class G>
struct AfforestSamplingTemplate {
  G& GA;
  Find& find;
  Unite& unite;
  uint32_t neighbor_rounds;

  AfforestSamplingTemplate(
      G& GA,
      Find& find,
      Unite& unite,
      commandLine& P) :
   GA(GA), find(find), unite(unite) {
    neighbor_rounds = P.getOptionLongValue("-sample_rounds", 2L);
   }

  pbbs::sequence<parent> initial_components() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    cout << "neighbor_rounds = " << neighbor_rounds << endl;

    auto parents = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
    pbbs::sequence<uintE> hooks;

    pbbs::random rnd;
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
          unite(u, ngh, parents);
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
//            unite(u, ngh, parents);
//          }
//        }, 512);
//      }
//      // compress nodes fully (turns out this is faster)
//      parallel_for(0, n, [&] (size_t u) {
//        parents[u] = find(u, parents);
//      }, 512);
//    }

    return parents;
   }
};

template <class W>
struct BFS_ComponentLabel_F {
  parent* Parents;
  uintE src;
  BFS_ComponentLabel_F(parent* _Parents, uintE src) : Parents(_Parents), src(src) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (Parents[d] != src) {
      Parents[d] = src;
      return 1;
    } else {
      return 0;
    }
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    return (pbbs::atomic_compare_and_swap(&Parents[d], d, src));
  }
  inline bool cond(const uintE& d) { return (Parents[d] == d); }
};

/* Returns a mapping from either i --> i, if i is not reached by the BFS, or
 * i --> src, if i is reachable from src in the BFS */
template <class Graph>
inline sequence<parent> BFS_ComponentLabel(Graph& G, uintE src) {
  using W = typename Graph::weight_type;
  /* Creates Parents array, initialized to all -1, except for src. */
  auto Parents = sequence<parent>(G.n, [&](size_t i) { return i; });
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
  // std::cout << "Reachable: " << reachable << " #rounds = " << rounds << std::endl;
  return Parents;
}

template <class G>
struct BFSSamplingTemplate {
  G& GA;

  BFSSamplingTemplate(G& GA, commandLine& P) :
   GA(GA) {}

  pbbs::sequence<parent> initial_components() {
    using W = typename G::weight_type;
    size_t n = GA.n;

    pbbs::sequence<parent> parents;

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

      parent frequent_comp; double pct;
      parents = std::move(bfs_parents);
      std::tie(frequent_comp, pct) = sample_frequent_element(parents);
      if (pct > static_cast<double>(0.1)) {
        std::cout << "# BFS covered: " << pct << " of graph" << std::endl;
        skip_comp = frequent_comp;
        found_massive_component = true;
        break;
      }
      std::cout << "# BFS covered only: " << pct << " of graph." << std::endl;
      rnd = rnd.next();
    }

    return parents;
  }
};

//    parallel_for(0, n, [&] (size_t i) {
//      if (!found_massive_component || parents[i] != skip_comp) {
//        parents[i] = i;
//      }
//    });
//
//    auto bfs_parents = parents; // copy
//
//    st.stop(); st.reportTotal("sample time");
//
//    timer ut; ut.start();
//    parallel_for(0, n, [&] (size_t u) {
//      // Only process edges for vertices not linked to the main component
//      // note that this is safe only for undirected graphs. For directed graphs,
//      // the in-edges must also be explored for all vertices.
//      if (bfs_parents[u] != skip_comp) {
//        auto map_f = [&] (uintE u, uintE v, const W& wgh) {
//          if (u < v) {
//            unite(u, v, parents);
//          }
//        };
//        GA.get_vertex(u).mapOutNgh(u, map_f); // in parallel
//      }
//    }, 1024);
//    ut.stop(); ut.reportTotal("union time");
//
//    timer ft; ft.start();
//    parallel_for(0, n, [&] (size_t i) {
//      parents[i] = find(i,parents);
//    });
//    ft.stop(); ft.reportTotal("find time");
//    return parents;
//   }


template <class G>
struct LDDSamplingTemplate {
  G& GA;

  LDDSamplingTemplate(G& GA, commandLine& P) : GA(GA) { }

  pbbs::sequence<parent> initial_components() {
    using W = typename G::weight_type;
    size_t n = GA.n;

    timer lddt; lddt.start();
    auto clusters_in = LDD(GA, 0.2, /* permute = */false);
    lddt.stop(); lddt.reportTotal("## ldd time");
    auto s = clusters_in.to_array();
    auto clusters = pbbs::sequence((parent*)s, n);

    return clusters;
  }
};

//   parallel_for(0, n, [&] (uintE u) {
//     if (clusters[u] == u) { // root, ok
//     } else {
//       assert(clusters[clusters[u]] == clusters[u]);
//     }
//   });
//
//    pbbs::sequence<parent> parents(n);
//    parallel_for(0, n, [&] (size_t i) {
//      parents[i] = clusters[i];
//    });
//
//    uintE frequent_comp; double pct;
//    std::tie(frequent_comp, pct) = sample_frequent_element(parents);
//    st.stop(); st.reportTotal("sample time");
//
//    timer ut; ut.start();
//    parallel_for(0, n, [&] (size_t u) {
//      // Only process edges for vertices not linked to the main component
//      // note that this is safe only for undirected graphs. For directed graphs,
//      // the in-edges must also be explored for all vertices.
//      if (clusters[u] != frequent_comp) {
//        auto map_f = [&] (uintE u, uintE v, const W& wgh) {
//          if (u < v) {
//            unite(u, v, parents);
//          }
//        };
//        GA.get_vertex(u).mapOutNgh(u, map_f); // in parallel
//      }
//    }, 512);
//    ut.stop(); ut.reportTotal("union time");
//
//    timer ft; ft.start();
//    parallel_for(0, n, [&] (size_t i) {
//      parents[i] = find(i,parents);
//    });
//    ft.stop(); ft.reportTotal("find time");
//    return parents;
//   }
//};

