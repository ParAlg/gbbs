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

#include <cmath>
#include "gbbs/bucket.h"
#include "gbbs/gbbs.h"

namespace gbbs {
namespace widestpath {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

struct Visit_F {
  sequence<uintE>& width;
  Visit_F(sequence<uintE>& _width) : width(_width) {}

  inline std::optional<uintE> update(const uintE& s, const uintE& d,
                                     const intE& w) {
    uintE oval = width[d];
    uintE bottleneck = oval | TOP_BIT;
    uintE n_width = std::min((width[s] | TOP_BIT), (w | TOP_BIT));
    if (n_width > bottleneck) {
      if (!(oval & TOP_BIT)) {  // First visitor
        width[d] = n_width;
        return std::optional<uintE>(oval);
      }
      width[d] = n_width;
    }
    return std::nullopt;
  }

  inline std::optional<uintE> updateAtomic(const uintE& s, const uintE& d,
                                           const intE& w) {
    uintE oval = width[d];
    uintE bottleneck = oval | TOP_BIT;
    uintE n_width = std::min((width[s] | TOP_BIT), (w | TOP_BIT));
    if (n_width > bottleneck) {
      if (!(oval & TOP_BIT) &&
          gbbs::atomic_compare_and_swap(&(width[d]), oval,
                                        n_width)) {  // First visitor
        return std::optional<uintE>(oval);
      }
      gbbs::write_max(&(width[d]), n_width);
    }
    return std::nullopt;
  }

  inline bool cond(const uintE& d) const { return true; }
};

}  // namespace widestpath

template <class Graph>
inline sequence<uintE> SSWidestPath(Graph& G, uintE src,
                                    size_t num_buckets = 128,
                                    bool largemem = false,
                                    bool no_blocked = false) {
  using W = typename Graph::weight_type;
  timer t;
  t.start();

  timer mw;
  mw.start();
  W max_weight = (W)0;
  parallel_for(0, G.n, 1, [&](size_t i) {
    auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
      if (wgh > max_weight) {
        gbbs::write_max(&max_weight, wgh);
      }
    };
    G.get_vertex(i).out_neighbors().map(map_f);
  });
  mw.stop();
  mw.next("max weight time");
  std::cout << "max_weight = " << max_weight << std::endl;

  timer init;
  init.start();
  size_t n = G.n;

  auto width =
      sequence<uintE>::from_function(n, [&](size_t i) { return (uintE)0; });
  width[src] = INT_E_MAX;

  auto get_bkt = [&](const W& _width) -> uintE {
    return max_weight - _width + 1;
  };
  auto get_ring = parlay::delayed_seq<uintE>(n, [&](const size_t& v) -> uintE {
    auto d = width[v];
    if (d == 0) {
      return UINT_E_MAX;
    }
    if (d == INT_E_MAX) {
      return 0;
    }
    return get_bkt(d);
  });
  std::cout << "creating bucket" << std::endl;
  auto b = make_vertex_buckets(n, get_ring, increasing, num_buckets);
  std::cout << "created bucket" << std::endl;

  auto apply_f = [&](const uintE v, uintE& old_width) -> void {
    uintE new_width = width[v] & widestpath::VAL_MASK;
    width[v] = new_width;  // Remove the TOP_BIT in the distance.
    // Compute the previous bucket and new bucket for the vertex.
    uintE prev_bkt = get_bkt(old_width);
    uintE new_bkt = get_bkt(new_width);
    uintE dest_bkt = b.get_bucket(prev_bkt, new_bkt);
    old_width = dest_bkt;  // write back
  };

  init.stop();
  init.next("init time");
  timer bt, emt;
  auto bkt = b.next_bucket();
  size_t rd = 0;
  flags fl = dense_forward;
  if (!largemem) fl |= no_dense;
  if (!no_blocked) fl |= sparse_blocked;
  while (bkt.id != b.null_bkt) {
    auto active = vertexSubset(n, std::move(bkt.identifiers));
    emt.start();
    auto em_f = wrap_with_default<W, W>(widestpath::Visit_F(width), (W)1);
    auto res = edgeMapData<uintE>(G, active, em_f, G.m / 20, fl);
    vertexMap(res, apply_f);
    // update buckets with vertices that just moved
    emt.stop();
    bt.start();
    if (res.dense()) {
      b.update_buckets(res.get_fn_repr(), n);
    } else {
      b.update_buckets(res.get_fn_repr(), res.size());
    }
    bkt = b.next_bucket();
    bt.stop();
    rd++;
  }
  bt.next("bucket time");
  emt.next("edge map time");
  std::cout << "n rounds = " << rd << "\n";

  auto dist_im_f = [&](size_t i) {
    return ((width[i] == INT_E_MAX) || (width[i] == (uintE)(-1))) ? 0
                                                                  : width[i];
  };  // noop?
  auto dist_im = parlay::delayed_seq<size_t>(n, dist_im_f);
  std::cout << "max dist = " << parlay::reduce_max(dist_im)
            << " xor = " << parlay::reduce_xor(dist_im) << "\n";
  return width;
}

struct SSWidestPathBF_F {
  intE* width;
  intE* Visited;
  SSWidestPathBF_F(intE* _width, intE* _Visited)
      : width(_width), Visited(_Visited) {}
  inline bool update(const uintE& s, const uintE& d, const intE& edgeLen) {
    intE n_width = std::min(width[s], edgeLen);  // new width of this path.
    if (width[d] > n_width) {  // update width to d if it becomes higher
      width[d] = n_width;
      if (Visited[d] == 0) {
        Visited[d] = 1;
        return 1;
      }
    }
    return 0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d,
                           const intE& edgeLen) {
    intE n_width = std::min(width[s], edgeLen);
    return (gbbs::write_max(&width[d], n_width) &&
            gbbs::atomic_compare_and_swap(&Visited[d], 0, 1));
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

// reset visited vertices
struct SSWidestPath_BF_Vertex_F {
  intE* Visited;
  SSWidestPath_BF_Vertex_F(intE* _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 0;
    return 1;
  }
};

template <class Graph>
inline sequence<intE> SSWidestPathBF(Graph& G, const uintE& start) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto Visited = sequence<int>(n, 0);
  auto width = sequence<intE>(n, static_cast<intE>(-1));
  width[start] = INT_E_MAX;  // width(s) = \infty

  vertexSubset Frontier(n, start);
  size_t round = 0;
  while (!Frontier.isEmpty()) {
    // Check for a negative weight cycle
    if (round == n) {
      parallel_for(0, n, kDefaultGranularity,
                   [&](size_t i) { width[i] = -(INT_E_MAX / 2); });
      break;
    }
    auto em_f = wrap_with_default<W, intE>(
        SSWidestPathBF_F(width.begin(), Visited.begin()), (intE)1);
    auto output =
        edgeMap(G, Frontier, em_f, G.m / 10, sparse_blocked | dense_forward);
    vertexMap(output, SSWidestPath_BF_Vertex_F(Visited.begin()));
    std::cout << output.size() << "\n";
    Frontier = std::move(output);
    round++;
  }
  auto dist_im_f = [&](size_t i) {
    return ((width[i] == INT_E_MAX) || (width[i] == static_cast<intE>(-1)))
               ? 0
               : width[i];
  };  // noop?
  auto dist_im = parlay::delayed_seq<size_t>(n, dist_im_f);
  std::cout << "max dist = " << parlay::reduce_max(dist_im)
            << " xor = " << parlay::reduce_xor(dist_im) << "\n";
  std::cout << "n rounds = " << round << "\n";
  return width;
}
}  // namespace gbbs
