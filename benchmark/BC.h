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
#include "bridge.h"
#include "ligra.h"
#include "edge_map_reduce.h"

namespace bc {

using fType = double;

template <class W, class S, class V>
struct BC_F {
  S& Scores;
  V& Visited;
  BC_F(S& _Scores, V& _Visited) : Scores(_Scores), Visited(_Visited) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    fType oldV = Scores[d];
    Scores[d] += Scores[s];
    return oldV == 0.0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    fType to_add = Scores[s];
    fType n_val = pbbslib::fetch_and_add(&Scores[d], to_add);
    return n_val == 0;
  }
  inline bool cond(uintE d) { return Visited[d] == 0; }
};

template <class W, class S, class V>
inline BC_F<W, S, V> make_bc_f(S& scores, V& visited) {
  return BC_F<W, S, V>(scores, visited);
}

// marks visited
template <class V>
struct BC_Vertex_F {
  V& Visited;
  BC_Vertex_F(V& _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 1;
    return 1;
  }
};

template <class V>
inline BC_Vertex_F<V> make_bc_vertex_f(V& visited) {
  return BC_Vertex_F<V>(visited);
}

// marks visited and adds to dependencies
template <class V, class D>
struct BC_Back_Vertex_F {
  V& Visited;
  D& Dependencies, &NumPaths;
  BC_Back_Vertex_F(V& _Visited, D& _Dependencies, D& _NumPaths)
      : Visited(_Visited), Dependencies(_Dependencies), NumPaths(_NumPaths) {}
  inline bool operator()(uintE i) {
    Visited[i] = 1;
    Dependencies[i] += NumPaths[i];
    return 1;
  }
};

template <class V, class D>
inline BC_Back_Vertex_F<V, D> make_bc_back_vertex_f(V& visited, D& dependencies,
                                                    D& num_paths) {
  return BC_Back_Vertex_F<V, D>(visited, dependencies, num_paths);
}

template <template <class W> class vertex, class W>
inline sequence<fType> BC(graph<vertex<W>>& GA, const uintE& start) {
  size_t n = GA.n;

  auto NumPaths = sequence<fType>(n, [](size_t i) { return 0.0; });
  NumPaths[start] = 1.0;

  auto Visited = sequence<bool>(n, [](size_t i) { return 0; });
  Visited[start] = 1;

  vertexSubset Frontier(n, start);

  std::vector<vertexSubset> Levels;

  long round = 0;
  while (!Frontier.isEmpty()) {
    debug(cout << "round = " << round << " fsize = " << Frontier.size() << endl;);
    round++;
    //      vertexSubset output = edgeMap(GA, Frontier,
    //      make_bc_f<W>(NumPaths,Visited), -1, sparse_blocked | dense_forward);
    vertexSubset output = edgeMap(GA, Frontier, make_bc_f<W>(NumPaths, Visited),
                                  -1, sparse_blocked | fine_parallel);
    vertexMap(output, make_bc_vertex_f(Visited));  // mark visited
    Levels.push_back(Frontier);                    // save frontier
    Frontier = output;
  }
  Levels.push_back(Frontier);

  auto Dependencies = sequence<fType>(n, [](size_t i) { return 0.0; });

  // Invert numpaths
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { NumPaths[i] = 1 / NumPaths[i]; });

  Levels[round].del();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { Visited[i] = 0; });
  Frontier = Levels[round - 1];
  vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));

  timer bt;
  bt.start();
  for (long r = round - 2; r >= 0; r--) {
    //      edgeMap(GA, Frontier, make_bc_f<W>(Dependencies,Visited), -1,
    //      no_output | in_edges | dense_forward);
    edgeMap(GA, Frontier, make_bc_f<W>(Dependencies, Visited), -1,
            no_output | in_edges | fine_parallel);
    Frontier.del();
    Frontier = Levels[r];
     vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));
  }
  bt.stop();
  debug(bt.reportTotal("back total time"););

  Frontier.del();

  // Update dependencies scores
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    Dependencies[i] = (Dependencies[i] - NumPaths[i]) / NumPaths[i];
  });
  return Dependencies;
}

template <template <class W> class vertex, class W, class E>
vertexSubset sparse_fa_dense_em(graph<vertex<W>>& GA, E& EM, vertexSubset& Frontier, pbbs::sequence<fType>& NumPaths, pbbs::sequence<fType>& Storage, pbbs::sequence<bool>& Visited,  const flags fl) {
  size_t out_degrees = 0;
  if (Frontier.dense()) {
    auto degree_f = [&](size_t i) -> size_t {
      if (Frontier.d[i]) {
        return (fl & in_edges) ? GA.get_vertex(i).getInDegree() : GA.get_vertex(i).getOutDegree();
      }
      return static_cast<size_t>(0);
    };
    auto degree_imap = pbbslib::make_sequence<size_t>(Frontier.size(), degree_f);
    out_degrees = pbbslib::reduce_add(degree_imap);
  } else {
    auto degree_f = [&](size_t i) -> size_t {
      return (fl & in_edges) ? GA.get_vertex(Frontier.vtx(i)).getInDegree() : GA.get_vertex(Frontier.vtx(i)).getOutDegree();
    };
    auto degree_imap = pbbslib::make_sequence<size_t>(Frontier.size(), degree_f);
    out_degrees = pbbslib::reduce_add(degree_imap);
  }

  if (out_degrees > GA.m/20) {
    debug(cout << "dense, out_degrees = " << out_degrees << endl;);

    auto cond_f = [&] (size_t i) {
      return (Visited[i] == false);
    };
    auto map_f = [&] (const uintE& s, const uintE& d, const W& wgh) -> double {
      return NumPaths[d];
    };
    auto reduce_f = [&] (double l, double r) { return l + r; };
    auto apply_f = [&] (std::tuple<uintE, double> k) {
      const uintE& u = std::get<0>(k);
      const double& contribution = std::get<1>(k);
      if (contribution > 0) {
        Storage[u] = contribution;
        return Maybe<std::tuple<uintE, pbbs::empty>>(std::make_tuple(u, pbbs::empty()));
      }
      return Maybe<std::tuple<uintE, pbbs::empty>>();
    };
    double id = 0.0;

    flags dense_fl = fl;
    dense_fl ^= in_edges; // should be set if out_edges, unset if in_edges
    timer dt; dt.start();
    vertexSubset output = EM.template edgeMapReduce_dense<pbbs::empty, double>(Frontier, cond_f, map_f, reduce_f, apply_f, id, dense_fl);

    parallel_for(0, GA.n, [&] (size_t i) {
      if (Storage[i] != 0) {
        NumPaths[i] = Storage[i];
        Storage[i] = 0;
      }
    });

    dt.stop(); dt.reportTotal("dense time");
    return output;
  } else {
    vertexSubset output = edgeMap(GA, Frontier, make_bc_f<W>(NumPaths, Visited),
                                  -1, fl | sparse_blocked | fine_parallel);
    return output;
  }
}

template <template <class W> class vertex, class W>
inline sequence<fType> BC_EM(graph<vertex<W>>& GA, const uintE& start) {
  size_t n = GA.n;
  auto EM = EdgeMap<fType, vertex, W>(GA, std::make_tuple(UINT_E_MAX, (fType)0.0), (size_t)GA.m/1000);

  auto NumPaths = sequence<fType>(n, static_cast<fType>(0));
  auto Storage = sequence<fType>(n, static_cast<fType>(0));
  NumPaths[start] = 1.0;

  auto Visited = sequence<bool>(n, [](size_t i) { return 0; });
  Visited[start] = 1;

  vertexSubset Frontier(n, start);

  std::vector<vertexSubset> Levels;

  timer fwd; fwd.start();
  long round = 0;
  while (!Frontier.isEmpty()) {
    debug(cout << "round = " << round << " fsize = " << Frontier.size() << endl;);
    round++;

    vertexSubset output = sparse_fa_dense_em(GA, EM, Frontier, NumPaths, Storage, Visited, 0);

    vertexMap(output, make_bc_vertex_f(Visited));  // mark visited
    Levels.push_back(Frontier);                    // save frontier
    Frontier = output;
  }
  Levels.push_back(Frontier);
  fwd.stop(); debug(fwd.reportTotal("forward time"));

  for (size_t i=0; i<100; i++) {
    cout << NumPaths[i] << endl;
  }
  cout << "printed numpaths" << endl;

  auto Dependencies = sequence<fType>(n, [](size_t i) { return 0.0; });

  // Invert numpaths
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { NumPaths[i] = 1 / NumPaths[i]; });

  Levels[round].del();
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { Visited[i] = 0; });
  Frontier = Levels[round - 1];
  vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));

  timer bt;
  bt.start();
  for (long r = round - 2; r >= 0; r--) {
    //      edgeMap(GA, Frontier, make_bc_f<W>(Dependencies,Visited), -1,
    //      no_output | in_edges | dense_forward);
    // edgeMap(GA, Frontier, make_bc_f<W>(Dependencies, Visited), -1,
    //         no_output | in_edges | fine_parallel);

    sparse_fa_dense_em(GA, EM, Frontier, Dependencies, Storage, Visited, in_edges | no_output);

    Frontier.del();
    Frontier = Levels[r];
    vertexMap(Frontier, make_bc_back_vertex_f(Visited, Dependencies, NumPaths));
  }
  bt.stop();
  debug(bt.reportTotal("back total time"););

  Frontier.del();

  // Update dependencies scores
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    Dependencies[i] = (Dependencies[i] - NumPaths[i]) / NumPaths[i];
  });
  return Dependencies;



}

}  // namespace bc
