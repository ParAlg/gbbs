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

#include <vector>
#include "gbbs/bridge.h"
#include "gbbs/edge_map_reduce.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace bc {

using fType = double;

template <class W, class S, class V>
struct SSBetweennessCentrality_F {
  S& Scores;
  V& Visited;
  SSBetweennessCentrality_F(S& _Scores, V& _Visited)
      : Scores(_Scores), Visited(_Visited) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    fType oldV = Scores[d];
    Scores[d] += Scores[s];
    return oldV == 0.0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    fType to_add = Scores[s];
    fType n_val = gbbs::fetch_and_add(&Scores[d], to_add);
    return n_val == 0;
  }
  inline bool cond(uintE d) { return Visited[d] == 0; }
};

template <class W, class S, class V>
inline SSBetweennessCentrality_F<W, S, V> make_bc_f(S& scores, V& visited) {
  return SSBetweennessCentrality_F<W, S, V>(scores, visited);
}

// marks visited
template <class V>
struct SSBetweennessCentrality_Vertex_F {
  V& Visited;
  SSBetweennessCentrality_Vertex_F(V& _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 1;
    return 1;
  }
};

template <class V>
inline SSBetweennessCentrality_Vertex_F<V> make_bc_vertex_f(V& visited) {
  return SSBetweennessCentrality_Vertex_F<V>(visited);
}

// marks visited and adds to dependencies
template <class V, class D>
struct SSBetweennessCentrality_Back_Vertex_F {
  V& Visited;
  D &Dependencies, &NumPaths;
  SSBetweennessCentrality_Back_Vertex_F(V& _Visited, D& _Dependencies,
                                        D& _NumPaths)
      : Visited(_Visited), Dependencies(_Dependencies), NumPaths(_NumPaths) {}
  inline bool operator()(uintE i) {
    Visited[i] = 1;
    Dependencies[i] += NumPaths[i];
    return 1;
  }
};

template <class V, class D>
inline SSBetweennessCentrality_Back_Vertex_F<V, D> make_bc_back_vertex_f(
    V& visited, D& dependencies, D& num_paths) {
  return SSBetweennessCentrality_Back_Vertex_F<V, D>(visited, dependencies,
                                                     num_paths);
}

template <class Graph>
inline sequence<fType> SSBetweennessCentrality(Graph& G, const uintE& start) {
  using W = typename Graph::weight_type;
  size_t n = G.n;

  auto NumPaths =
      sequence<fType>::from_function(n, [](size_t i) { return 0.0; });
  NumPaths[start] = 1.0;

  auto Visited = sequence<bool>::from_function(n, [](size_t i) { return 0; });
  Visited[start] = 1;

  vertexSubset Frontier(n, start);

  std::vector<vertexSubset> Levels;

  long round = 0;
  while (!Frontier.isEmpty()) {
    debug(std::cout << "round = " << round << " fsize = " << Frontier.size()
                    << std::endl;);
    round++;
    //      vertexSubset output = edgeMap(G, Frontier,
    //      make_bc_f<W>(NumPaths,Visited), -1, sparse_blocked | dense_forward);
    vertexSubset output = edgeMap(G, Frontier, make_bc_f<W>(NumPaths, Visited),
                                  -1, sparse_blocked | fine_parallel);
    vertexMap(output, make_bc_vertex_f(Visited));  // mark visited
    Levels.push_back(std::move(Frontier));         // save frontier
    Frontier = std::move(output);
  }
  Levels.push_back(std::move(Frontier));

  auto Dependencies =
      sequence<fType>::from_function(n, [](size_t i) { return 0.0; });

  // Invert numpaths
  parallel_for(0, n, kDefaultGranularity,
               [&](size_t i) { NumPaths[i] = 1 / NumPaths[i]; });

  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { Visited[i] = 0; });
  Frontier = std::move(Levels[round - 1]);
  vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));

  timer bt;
  bt.start();
  for (long r = round - 2; r >= 0; r--) {
    //      edgeMap(G, Frontier, make_bc_f<W>(Dependencies,Visited), -1,
    //      no_output | in_edges | dense_forward);
    edgeMap(G, Frontier, make_bc_f<W>(Dependencies, Visited), -1,
            no_output | in_edges | fine_parallel);
    Frontier = std::move(Levels[r]);
    vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));
  }
  bt.stop();
  debug(bt.next("back total time"););

  // Update dependencies scores
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    Dependencies[i] = (Dependencies[i] - NumPaths[i]) / NumPaths[i];
  });
  return Dependencies;
}

template <class Graph, class E>
vertexSubset sparse_fa_dense_em(Graph& G, E& EM, vertexSubset& Frontier,
                                sequence<fType>& NumPaths,
                                sequence<fType>& Storage,
                                sequence<bool>& Visited, const flags fl) {
  using W = typename Graph::weight_type;
  size_t out_degrees = 0;
  if (Frontier.dense()) {
    auto degree_f = [&](size_t i) -> size_t {
      if (Frontier.d[i]) {
        return (fl & in_edges)
                   ? G.get_vertex(i).in_neighbors().get_virtual_degree()
                   : G.get_vertex(i).out_neighbors().get_virtual_degree();
      }
      return static_cast<size_t>(0);
    };
    auto degree_imap = parlay::delayed_seq<size_t>(Frontier.size(), degree_f);
    out_degrees = parlay::reduce(degree_imap);
  } else {
    auto degree_f = [&](size_t i) -> size_t {
      return (fl & in_edges)
                 ? G.get_vertex(i).in_neighbors().get_virtual_degree()
                 : G.get_vertex(i).out_neighbors().get_virtual_degree();
    };
    auto degree_imap = parlay::delayed_seq<size_t>(Frontier.size(), degree_f);
    out_degrees = parlay::reduce(degree_imap);
  }

  if (out_degrees > G.m / 20) {
    debug(std::cout << "dense, out_degrees = " << out_degrees << std::endl;);

    auto cond_f = [&](size_t i) { return (Visited[i] == false); };
    auto map_f = [&](const uintE& s, const uintE& d, const W& wgh) -> double {
      return NumPaths[d];
    };
    auto reduce_f = [&](double l, double r) { return l + r; };
    auto apply_f = [&](std::tuple<uintE, double> k)
        -> std::optional<std::tuple<uintE, gbbs::empty>> {
          const uintE& u = std::get<0>(k);
          const double& contribution = std::get<1>(k);
          if (contribution > 0) {
            Storage[u] = contribution;
            return std::optional<std::tuple<uintE, gbbs::empty>>(
                {u, gbbs::empty()});
          }
          return std::nullopt;
        };
    double id = 0.0;

    flags dense_fl = fl;
    dense_fl ^= in_edges;  // should be set if out_edges, unset if in_edges
    timer dt;
    dt.start();
    vertexSubset output = EM.template edgeMapReduce_dense<gbbs::empty, double>(
        Frontier, cond_f, map_f, reduce_f, apply_f, id, dense_fl);

    parallel_for(0, G.n, [&](size_t i) {
      if (Storage[i] != 0) {
        NumPaths[i] = Storage[i];
        Storage[i] = 0;
      }
    });

    dt.stop();
    dt.next("dense time");
    return output;
  } else {
    vertexSubset output = edgeMap(G, Frontier, make_bc_f<W>(NumPaths, Visited),
                                  -1, fl | sparse_blocked | fine_parallel);
    return output;
  }
}

template <class Graph>
inline sequence<fType> SSBetweennessCentrality_EM(Graph& G,
                                                  const uintE& start) {
  size_t n = G.n;
  auto EM = EdgeMap<fType, Graph>(G, std::make_tuple(UINT_E_MAX, (fType)0.0),
                                  (size_t)G.m / 1000);

  auto NumPaths = sequence<fType>(n, static_cast<fType>(0));
  auto Storage = sequence<fType>(n, static_cast<fType>(0));
  NumPaths[start] = 1.0;

  auto Visited = sequence<bool>::from_function(n, [](size_t i) { return 0; });
  Visited[start] = 1;

  vertexSubset Frontier(n, start);

  std::vector<vertexSubset> Levels;

  timer fwd;
  fwd.start();
  long round = 0;
  while (!Frontier.isEmpty()) {
    debug(std::cout << "round = " << round << " fsize = " << Frontier.size()
                    << std::endl;);
    round++;

    vertexSubset output =
        sparse_fa_dense_em(G, EM, Frontier, NumPaths, Storage, Visited, 0);

    vertexMap(output, make_bc_vertex_f(Visited));  // mark visited
    Levels.push_back(std::move(Frontier));         // save frontier
    Frontier = std::move(output);
  }
  Levels.push_back(std::move(Frontier));
  fwd.stop();
  debug(fwd.next("forward time"));

  for (size_t i = 0; i < 100; i++) {
    std::cout << NumPaths[i] << std::endl;
  }
  std::cout << "printed numpaths" << std::endl;

  auto Dependencies =
      sequence<fType>::from_function(n, [](size_t i) { return 0.0; });

  // Invert numpaths
  parallel_for(0, n, kDefaultGranularity,
               [&](size_t i) { NumPaths[i] = 1 / NumPaths[i]; });

  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { Visited[i] = 0; });
  Frontier = std::move(Levels[round - 1]);
  vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));

  timer bt;
  bt.start();
  for (long r = round - 2; r >= 0; r--) {
    //      edgeMap(G, Frontier, make_bc_f<W>(Dependencies,Visited), -1,
    //      no_output | in_edges | dense_forward);
    // edgeMap(G, Frontier, make_bc_f<W>(Dependencies, Visited), -1,
    //         no_output | in_edges | fine_parallel);

    sparse_fa_dense_em(G, EM, Frontier, Dependencies, Storage, Visited,
                       in_edges | no_output);

    Frontier = std::move(Levels[r]);
    vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));
  }
  bt.stop();
  debug(bt.next("back total time"););

  // Update dependencies scores
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    Dependencies[i] = (Dependencies[i] - NumPaths[i]) / NumPaths[i];
  });
  return Dependencies;
}

}  // namespace bc

namespace bc_bfs {

using fType = double;

template <class W>
struct BFS_F {
  sequence<uint8_t>& Visited;
  BFS_F(sequence<uint8_t>& Visited) : Visited(Visited) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    Visited[d] = 1; /* first visit */
    return 1;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
    return (gbbs::atomic_compare_and_swap(&Visited[d], (uint8_t)0,
                                          (uint8_t)1)); /* first visit */
  }
  inline bool cond(const uintE& d) { return (Visited[d] == 0); }
};

template <class Graph>
inline sequence<fType> SSBetweennessCentrality_BFS(Graph& G,
                                                   const uintE& start) {
  using W = typename Graph::weight_type;
  size_t n = G.n;

  auto NumPaths = sequence<fType>(n, static_cast<fType>(0));
  NumPaths[start] = 1.0;

  /* 0 = unvisited
   * 1 = first visit
   * 2 = finished */
  auto Visited =
      sequence<uint8_t>::from_function(n, [](size_t i) { return 0; });
  Visited[start] = 2;

  vertexSubset Frontier(n, start);
  std::vector<vertexSubset> Levels;

  /* Forward pass */
  timer fwd;
  fwd.start();
  long round = 0;
  {
    auto map_f = [&](const uintE& s, const uintE& d, const W& wgh) -> fType {
      if (Visited[d] == 2) {
        return NumPaths[d]; /* return values from finished vertices */
      }
      return (fType)0;
    };
    auto reduce_f = [&](const fType& l, const fType& r) -> fType {
      return l + r;
    };
    auto id = (fType)0;
    auto monoid_f = parlay::make_monoid(reduce_f, id);

    auto reduce_incident_edges = [&](vertexSubset& vs, flags fl) {
      vertexMap(vs, [&](const uintE& u) {
        NumPaths[u] =
            (fl & in_edges)
                ? G.get_vertex(u).in_neighbors().reduce(map_f, monoid_f)
                : G.get_vertex(u).out_neighbors().reduce(map_f, monoid_f);
      });
    };
    while (!Frontier.isEmpty()) {
      debug(std::cout << "round = " << round << " fsize = " << Frontier.size()
                      << std::endl;);
      round++;

      vertexSubset next_frontier = edgeMap(G, Frontier, BFS_F<W>(Visited), -1,
                                           sparse_blocked | dense_parallel);

      reduce_incident_edges(next_frontier, in_edges);

      vertexMap(next_frontier, [&](const uintE u) {
        Visited[u] = 2; /* finished */
      });

      Levels.push_back(std::move(Frontier));  // save frontier
      Frontier = std::move(next_frontier);
    }
  }
  Levels.push_back(std::move(Frontier));
  fwd.stop();
  fwd.next("forward time");

  /* Backwards pass */
  auto Dependencies =
      sequence<fType>::from_function(n, [](size_t i) { return 0.0; });
  // Invert numpaths
  parallel_for(0, n, kDefaultGranularity,
               [&](size_t i) { NumPaths[i] = 1 / NumPaths[i]; });

  parallel_for(0, n, kDefaultGranularity, [&](size_t i) { Visited[i] = 0; });
  Frontier = std::move(Levels[round - 1]);
  std::cout << "r-1 frontier, m = " << Frontier.m << std::endl;

  timer bt;
  bt.start();
  {
    auto finish_vertex_subset = [&](vertexSubset& vs) {
      vertexMap(Frontier, [&](const uintE& u) {
        Visited[u] = 2; /* finished */
        Dependencies[u] += NumPaths[u];
      });
    };
    finish_vertex_subset(Frontier);

    auto map_f = [&](const uintE& s, const uintE& d, const W& wgh) -> fType {
      if (Visited[d] == 2) {
        return Dependencies[d]; /* return values from finished vertices */
      }
      return (fType)0;
    };
    auto reduce_f = [&](const fType& l, const fType& r) -> fType {
      return l + r;
    };
    auto id = (fType)0;
    auto monoid_f = parlay::make_monoid(reduce_f, id);

    auto reduce_dependencies = [&](vertexSubset& vs) {
      vertexMap(vs, [&](const uintE& u) {
        Dependencies[u] =
            G.get_vertex(u).out_neighbors().reduce(map_f, monoid_f);
      });
    };
    for (long r = round - 2; r >= 0; r--) {
      Frontier = std::move(Levels[r]);

      reduce_dependencies(Frontier);

      finish_vertex_subset(Frontier);
    }
  }
  bt.stop();
  debug(bt.next("back total time"););

  // Update dependencies scores
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    Dependencies[i] = (Dependencies[i] - NumPaths[i]) / NumPaths[i];
  });
  return Dependencies;
}

}  // namespace bc_bfs
}  // namespace gbbs
