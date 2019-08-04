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
#include "bucket.h"
#include "ligra.h"

namespace widestpath {
constexpr uintE TOP_BIT = ((uintE)INT_E_MAX) + 1;
constexpr uintE VAL_MASK = INT_E_MAX;

template <class W, class GW>
struct Visit_F {
  sequence<uintE>& width;
  GW& get_weight;
  Visit_F(sequence<uintE>& _width, GW& get_weight) : width(_width), get_weight(get_weight) {}

  inline Maybe<uintE> update(const uintE& s, const uintE& d, const W& w) {
    uintE edgeLen = get_weight(s, d, w);
    uintE oval = width[d];
    uintE bottleneck = oval | TOP_BIT;
    uintE n_width = std::min((width[s] | TOP_BIT), (edgeLen | TOP_BIT));
    if (n_width > bottleneck) {
      if (!(oval & TOP_BIT)) {  // First visitor
        width[d] = n_width;
        return Maybe<uintE>(oval);
      }
      width[d] = n_width;
    }
    return Maybe<uintE>();
  }

  inline Maybe<uintE> updateAtomic(const uintE& s, const uintE& d,
                                   const W& w) {
    uintE edgeLen = get_weight(s, d, w);
    uintE oval = width[d];
    uintE bottleneck = oval | TOP_BIT;
    uintE n_width = std::min((width[s] | TOP_BIT), (edgeLen | TOP_BIT));
    if (n_width > bottleneck) {
      if (!(oval & TOP_BIT) &&
          pbbslib::atomic_compare_and_swap(&(width[d]), oval, n_width)) {  // First visitor
        return Maybe<uintE>(oval);
      }
      pbbslib::write_max(&(width[d]), n_width);
    }
    return Maybe<uintE>();
  }

  inline bool cond(const uintE& d) const { return true; }
};

template <class W, class GW>
Visit_F<W, GW> make_visit_f(sequence<uintE>& dists, GW& get_weight) {
  return Visit_F<W, GW>(dists, get_weight);
}

}  // namespace widestpath

template <class G, class GW>
inline sequence<uintE> WidestPath(G& GA, GW& get_weight, uintE src,
                              size_t num_buckets = 128, bool largemem = false,
                              bool no_blocked = false) {
  using W = typename G::weight_type;
  timer t;
  t.start();

  timer mw; mw.start();
  W max_weight = 0;
  parallel_for(0, GA.n, [&] (size_t i) {
    auto map_f = [&] (const uintE& u, const uintE& v, const W& w) {
      auto wgh = get_weight(u, v, w);
      if (wgh > max_weight) {
        pbbslib::write_max(&max_weight, wgh);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f);
  }, 1);
  mw.stop(); mw.reportTotal("max weight time");
  cout << "max_weight = " << max_weight << endl;

  timer init;
  init.start();
  size_t n = GA.n;

  auto width = sequence<uintE>(n, [&](size_t i) { return 0; });
  width[src] = INT_E_MAX;

  auto get_bkt = [&](const uintE& width) -> const uintE {
    return max_weight - width + 1;
  };
  auto get_ring = pbbslib::make_sequence<uintE>(n, [&](const size_t& v) -> const uintE {
    auto d = width[v];
    if (d == 0) { return UINT_E_MAX; }
    if (d == INT_E_MAX) { return 0; }
    return get_bkt(d);
  });
  cout << "creating bucket" << endl;
  auto b = make_vertex_buckets(n, get_ring, increasing, num_buckets);
  cout << "created bucket" << endl;

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
  init.reportTotal("init time");
  timer bt, emt;
  auto bkt = b.next_bucket();
  size_t rd = 0;
  flags fl = dense_forward;
  if (!largemem) fl |= no_dense;
  if (!no_blocked) fl |= sparse_blocked;
  while (bkt.id != b.null_bkt) {
    auto active = vertexSubset(n, bkt.identifiers);
    emt.start();
    auto res = edgeMapData<uintE>(GA, active, widestpath::make_visit_f<W>(width, get_weight), GA.m / 20, fl);
    vertexMap(res, apply_f);
    // update buckets with vertices that just moved
    emt.stop();
    bt.start();
    if (res.dense()) {
      b.update_buckets(res.get_fn_repr(), n);
    } else {
      b.update_buckets(res.get_fn_repr(), res.size());
    }
    res.del();
    active.del();
    bkt = b.next_bucket();
    bt.stop();
    rd++;
  }
  bt.reportTotal("bucket time");
  emt.reportTotal("edge map time");
  std::cout << "n rounds = " << rd << "\n";

  auto dist_im_f = [&](size_t i) { return ((width[i] == INT_E_MAX) || (static_cast<intE>(width[i]) == static_cast<intE>(-1))) ? 0 : width[i]; }; // noop?
  auto dist_im = pbbslib::make_sequence<size_t>(n, dist_im_f);
  std::cout << "max dist = " << pbbslib::reduce_max(dist_im) << " xor = " << pbbslib::reduce_xor(dist_im) << "\n";
  for (size_t i=0; i<100; i++) {
    cout << dist_im[i] << endl;
  }
  return width;
}

struct WidestPathBF_F {
  intE* width;
  intE* Visited;
  WidestPathBF_F(intE* _width, intE* _Visited) : width(_width), Visited(_Visited) {}
  inline bool update(const uintE& s, const uintE& d, const intE& edgeLen) {
    intE n_width = std::min(width[s], edgeLen); // new width of this path.
    if (width[d] > n_width) { // update width to d if it becomes higher
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
    return (pbbslib::write_max(&width[d], n_width) && pbbslib::atomic_compare_and_swap(&Visited[d], 0, 1));
  }
  inline bool cond(uintE d) { return cond_true(d); }
};

// reset visited vertices
struct WidestPath_BF_Vertex_F {
  intE* Visited;
  WidestPath_BF_Vertex_F(intE* _Visited) : Visited(_Visited) {}
  inline bool operator()(uintE i) {
    Visited[i] = 0;
    return 1;
  }
};

template <class G>
inline sequence<intE> WidestPathBF(G& GA, const uintE& start) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  auto Visited = sequence<int>(n, 0);
  auto width = sequence<intE>(n, static_cast<intE>(-1));
  width[start] = INT_E_MAX; // width(s) = \infty

  vertexSubset Frontier(n, start);
  size_t round = 0;
  while (!Frontier.isEmpty()) {
    // Check for a negative weight cycle
    if (round == n) {
      par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                      { width[i] = -(INT_E_MAX / 2); });
      break;
    }
    auto em_f =
        wrap_with_default<W, intE>(WidestPathBF_F(width.begin(), Visited.begin()), (intE)1);
    auto output =
        edgeMap(GA, Frontier, em_f, GA.m / 10, sparse_blocked | dense_forward);
    vertexMap(output, WidestPath_BF_Vertex_F(Visited.begin()));
    std::cout << output.size() << "\n";
    Frontier.del();
    Frontier = output;
    round++;
  }
  auto dist_im_f = [&](size_t i) { return ((width[i] == INT_E_MAX) || (width[i] == static_cast<intE>(-1))) ? 0 : width[i]; }; // noop?
  auto dist_im = pbbslib::make_sequence<size_t>(n, dist_im_f);
  std::cout << "max dist = " << pbbslib::reduce_max(dist_im) << " xor = " << pbbslib::reduce_xor(dist_im) << "\n";
  std::cout << "n rounds = " << round << "\n";
  return width;
}
