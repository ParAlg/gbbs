#pragma once

namespace gbbs {
template <class W>
struct BFS_ComponentLabel_F {
  sequence<parent>& Parents;
  uintE src;
  BFS_ComponentLabel_F(sequence<parent>& _Parents, uintE src)
      : Parents(_Parents), src(src) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (Parents[d] != src) {
      Parents[d] = src;
      return 1;
    } else {
      return 0;
    }
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    if (Parents[d] != src &&
        gbbs::atomic_compare_and_swap(&Parents[d], static_cast<parent>(d),
                                      static_cast<parent>(src))) {
      return true;
    }
    return false;
  }
  inline bool cond(const uintE& d) { return (Parents[d] == d); }
};

/* Returns a mapping from either i --> i, if i is not reached by the BFS, or
 * i --> src, if i is reachable from src in the BFS */
template <class Graph>
inline auto BFS_ComponentLabel(Graph& G, uintE src) {
  using W = typename Graph::weight_type;

  auto Edges = sequence<edge>(G.n, empty_edge);
  auto Parents = sequence<parent>(G.n, [&](size_t i) { return i; });
  Parents[src] = src;

  vertexSubset Frontier(G.n, src);
  size_t reachable = 0;
  size_t rounds = 0;
  while (!Frontier.isEmpty()) {
    reachable += Frontier.size();
    vertexSubset output =
        edgeMap(G, Frontier, BFS_ComponentLabel_F<W>(Parents, src), -1,
                sparse_blocked | dense_parallel);
    Frontier = std::move(output);
    rounds++;
  }
  parallel_for(0, G.n, [&](size_t i) {
    uintE par = Parents[i];
    if (par != i) {
      Edges[i] = std::make_pair(i, Parents[i]);
    }
  });
  // std::cout << "Reachable: " << reachable << " #rounds = " << rounds <<
  // std::endl;
  return std::make_pair(std::move(Parents), std::move(Edges));
}

template <class G>
struct BFSSamplingTemplate {
  G& GA;

  BFSSamplingTemplate(G& GA, commandLine& P) : GA(GA) {}

  auto initial_spanning_forest() {
    size_t n = GA.n;

    sequence<parent> Parents;
    sequence<edge> Edges;

    parlay::random rnd;
    timer st;
    st.start();

    uint32_t max_trials = 3;
    for (uint32_t r = 0; r < max_trials; r++) {
      uintE src = rnd.rand() % n;
      auto[bfs_parents, bfs_edges] = BFS_ComponentLabel(GA, src);

      parent frequent_comp;
      double pct;
      Parents = std::move(bfs_parents);
      Edges = std::move(bfs_edges);

      std::tie(frequent_comp, pct) =
          connectit::sample_frequent_element(Parents);
      if (pct > static_cast<double>(0.1)) {
        std::cout << "# BFS from " << src << " covered: " << pct << " of graph"
                  << std::endl;
        break;
      }
      std::cout << "# BFS covered only: " << pct << " of graph." << std::endl;
      rnd = rnd.next();
    }

    return std::make_tuple(Parents, Edges);
  }
};
}  // namespace gbbs
