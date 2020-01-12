#pragma once

#include "benchmarks/Connectivity/connectit.h"
#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"

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
      timer sample_t; sample_t.start();
      auto parents = sampler.initial_components();
      sample_t.stop();
      sample_t.reportTotal("sample time");

      parent frequent_comp; double pct;
      std::tie(frequent_comp, pct) = connectit::sample_frequent_element(parents);

      algorithm.initialize(parents);

      /* relabel for liu_tarjan */
      if constexpr (algorithm_type == liu_tarjan_type || algorithm_type == label_prop_type) {
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
template <class G>
struct KOutSamplingTemplate {
  G& GA;
  uint32_t neighbor_rounds;

  KOutSamplingTemplate(
      G& GA,
      commandLine& P) :
   GA(GA) {
    neighbor_rounds = P.getOptionLongValue("-sample_rounds", 2L);
   }

  void link(uintE u, uintE v, pbbs::sequence<parent>& parents) {
    parent p1 = parents[u];
    parent p2 = parents[v];
    while (p1 != p2) {
      parent high = p1 > p2 ? p1 : p2;
      parent low = p1 + (p2 - high);
      parent p_high = parents[high];
      // Was already 'low' or succeeded in writing 'low'
      if ((p_high == low) ||
          (p_high == high && pbbs::atomic_compare_and_swap(&parents[high], high, low)))
        break;
      p1 = parents[parents[high]];
      p2 = parents[low];
    }
  }

  pbbs::sequence<parent> initial_components() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    cout << "# neighbor_rounds = " << neighbor_rounds << endl;

    auto parents = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
    pbbs::sequence<uintE> hooks;

    pbbs::random rnd;
    uintE granularity = 1024;
    for (uint32_t r=0; r<neighbor_rounds; r++) {
      if (r == 0) {
        /* First round: sample a directed forest and compress instead of using
         * unite */
        parallel_for(0, n, [&] (size_t u) {
          auto u_vtx = GA.get_vertex(u);
          if (u_vtx.getOutDegree() > r) {
            uintE ngh; W wgh;
            std::tie(ngh, wgh) = u_vtx.get_ith_out_neighbor(u, r);
            link(u, ngh, parents);
          }
        }, granularity);
      } else {
        parallel_for(0, n, [&] (size_t u) {
          auto u_rnd = rnd.fork(u);
          auto u_vtx = GA.get_vertex(u);
          if (u_vtx.getOutDegree() > 0) {
            uintE deglb = (1 << std::max((int)pbbs::log2_up(u_vtx.getOutDegree()), (int)1) - 1) - 1;
            uintE ngh_idx;
            if (deglb == 0) {
              ngh_idx = 0;
            } else {
              ngh_idx = u_rnd.rand() & deglb;
            }
            uintE ngh; W wgh;
            std::tie(ngh, wgh) = u_vtx.get_ith_out_neighbor(u, ngh_idx);
            link(u, ngh, parents);
          }
        }, granularity);
      }
      // compress nodes fully (turns out this is faster)
      parallel_for(0, n, [&] (size_t u) {
        while (parents[u] != parents[parents[u]]) {
          parents[u] = parents[parents[u]];
        }
      }, granularity);
      rnd = rnd.next();
    }
    return parents;
   }

  pbbs::sequence<parent> initial_components_pure() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    std::cout << "# neighbor_rounds = " << neighbor_rounds << std::endl;

    auto parents = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
    pbbs::sequence<uintE> hooks;

    pbbs::random rnd;
    uintE granularity = 1024;
    // Using random neighbor---some overhead (and faster for some graphs), also
    // theoretically defensible
    for (uint32_t r=0; r<neighbor_rounds; r++) {
      parallel_for(0, n, [&] (size_t u) {
        auto u_rnd = rnd.fork(u);
        auto u_vtx = GA.get_vertex(u);
        auto out_degree = u_vtx.getOutDegree();
        if (out_degree > 0) {
          uintE ngh_idx = u_rnd.rand() % out_degree;
          auto [ngh, wgh] = u_vtx.get_ith_out_neighbor(u, ngh_idx);
          link(u, ngh, parents);
        }
      }, granularity);
      // compress nodes fully
      parallel_for(0, n, [&] (size_t u) {
        while (parents[u] != parents[parents[u]]) {
          parents[u] = parents[parents[u]];
        }
      }, granularity);
      rnd = rnd.next();
    }
    return parents;
   }

  pbbs::sequence<parent> initial_components_afforest() {
    using W = typename G::weight_type;
    size_t n = GA.n;
    std::cout << "# neighbor_rounds = " << neighbor_rounds << std::endl;

    auto parents = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
    pbbs::sequence<uintE> hooks;

    uintE granularity = 1024;
    // Using random neighbor---some overhead (and faster for some graphs), also
    // theoretically defensible
    for (uint32_t r=0; r<neighbor_rounds; r++) {
      parallel_for(0, n, [&] (size_t u) {
        auto u_vtx = GA.get_vertex(u);
        auto out_degree = u_vtx.getOutDegree();
        if (out_degree > r) {
          auto [ngh, wgh] = u_vtx.get_ith_out_neighbor(u, r);
          link(u, ngh, parents);
        }
      }, granularity);
      // compress nodes fully
      parallel_for(0, n, [&] (size_t u) {
        while (parents[u] != parents[parents[u]]) {
          parents[u] = parents[parents[u]];
        }
      }, granularity);
    }
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
    return (pbbs::atomic_compare_and_swap(&Parents[d], static_cast<parent>(d), static_cast<parent>(src)));
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
    size_t n = GA.n;

    pbbs::sequence<parent> parents;

    pbbs::random rnd;
    timer st; st.start();

    /* Assume that G has a massive component with size at least 10% of the
     * graph. Run BFSs at most max_trials times until we find it. */

    uint32_t max_trials = 3;
    for (uint32_t r=0; r<max_trials; r++) {
      uintE src = rnd.rand() % n;
      auto bfs_parents = BFS_ComponentLabel(GA, src);

      parent frequent_comp; double pct;
      parents = std::move(bfs_parents);
      std::tie(frequent_comp, pct) = connectit::sample_frequent_element(parents);
      if (pct > static_cast<double>(0.1)) {
        std::cout << "# BFS covered: " << pct << " of graph" << std::endl;
        break;
      }
      std::cout << "# BFS covered only: " << pct << " of graph." << std::endl;
      rnd = rnd.next();
    }

    return parents;
  }
};


template <class G>
struct LDDSamplingTemplate {
  G& GA;

  LDDSamplingTemplate(G& GA, commandLine& P) : GA(GA) { }

  pbbs::sequence<parent> initial_components() {
    size_t n = GA.n;

    timer lddt; lddt.start();
    auto clusters_in = LDD(GA, 0.2, /* permute = */false);
    lddt.stop(); lddt.reportTotal("## ldd time");
    auto s = clusters_in.to_array();
    auto clusters = pbbs::sequence((parent*)s, n);

    return clusters;
  }
};

