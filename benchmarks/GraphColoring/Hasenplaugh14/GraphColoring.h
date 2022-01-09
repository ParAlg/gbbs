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

#include "gbbs/gbbs.h"

namespace gbbs {
namespace coloring {
template <class Graph, class Seq>
inline uintE color(Graph& G, uintE v, Seq& colors) {
  using W = typename Graph::weight_type;
  uintE deg = G.get_vertex(v).out_degree();
  if (deg > 0) {
    bool* bits;
    bool s_bits[1000];
    if (deg > 1000)
      bits = gbbs::new_array_no_init<bool>(deg);
    else
      bits = (bool*)s_bits;

    parallel_for(0, deg, kDefaultGranularity, [&](size_t i) { bits[i] = 0; });
    auto map_f = [&](uintE src, uintE ngh, const W& wgh) {
      uintE color = colors[ngh];
      if (color < deg) {
        bits[color] = 1;
      }
    };
    G.get_vertex(v).out_neighbors().map(map_f);
    auto im_f = [&](size_t i) {
      return (bits[i] == 0) ? (uintE)i : UINT_E_MAX;
    };
    auto im = parlay::delayed_seq<uintE>(deg, im_f);
    uintE color = parlay::reduce(im, parlay::minm<uintE>());
    if (deg > 1000) {
      gbbs::free_array(bits, deg);
    }
    return (color == UINT_E_MAX) ? (deg + 1) : color;
  }
  return 0;
}
}  // namespace coloring

template <class W>
struct coloring_f {
  intE* p;
  coloring_f(intE* _p) : p(_p) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (s == d) {
      std::cout << "error"
                << "\n";
      exit(-1);
    }
    p[d]--;
    return p[d] == 0;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (s == d) {
      std::cout << "error"
                << "\n";
      exit(-1);
    }
    return (gbbs::fetch_and_add(&p[d], -1) == 1);
  }
  inline bool cond(uintE d) { return (p[d] > 0); }
};

template <class Graph>
inline sequence<uintE> Coloring(Graph& G, bool lf = false) {
  using W = typename Graph::weight_type;
  timer initt;
  initt.start();
  const size_t n = G.n;

  // For each vertex count the number of out-neighbors with log-degree >= us
  auto priorities = sequence<intE>(n);
  auto colors =
      sequence<uintE>::from_function(n, [](size_t i) { return UINT_E_MAX; });

  if (lf) {
    std::cout << "### Running LF"
              << "\n";
    // LF heuristic
    auto P = parlay::random_permutation<uintE>(n);
    parallel_for(0, n, 1, [&](size_t i) {
      uintE our_deg = G.get_vertex(i).out_degree();
      uintE i_p = P[i];
      auto count_f = [&](uintE src, uintE ngh, const W& wgh) {
        uintE ngh_deg = G.get_vertex(ngh).out_degree();
        return (ngh_deg > our_deg) || ((ngh_deg == our_deg) && P[ngh] < i_p);
      };
      priorities[i] = G.get_vertex(i).out_neighbors().count(count_f);
    });
  } else {
    std::cout << "### Running LLF"
              << "\n";
    // LLF heuristic
    auto P = parlay::random_permutation<uintE>(n);
    parallel_for(0, n, 1, [&](size_t i) {
      uintE our_deg = parlay::log2_up(G.get_vertex(i).out_degree());
      uintE i_p = P[i];
      // breaks ties using P
      auto count_f = [&](uintE src, uintE ngh, const W& wgh) {
        uintE ngh_deg = parlay::log2_up(G.get_vertex(ngh).out_degree());
        return (ngh_deg > our_deg) || ((ngh_deg == our_deg) && P[ngh] < i_p);
      };
      priorities[i] = G.get_vertex(i).out_neighbors().count(count_f);
    });
  }

  auto zero_map_f = [&](size_t i) { return priorities[i] == 0; };
  auto zero_map = parlay::delayed_seq<bool>(n, zero_map_f);
  auto roots = vertexSubset(n, parlay::pack_index<uintE>(zero_map));
  debug(initt.next("init time"););

  size_t finished = 0, rounds = 0;
  timer color_t;
  timer em_t;
  while (finished != n) {
    assert(roots.size() > 0);
    finished += roots.size();
    roots.toSparse();

    // color the rootset
    color_t.start();
    parallel_for(0, roots.size(), 1, [&](size_t i) {
      uintE v = roots.vtx(i);
      colors[v] = coloring::color(G, v, colors);
    });
    color_t.stop();

    // compute the new rootset
    em_t.start();
    auto new_roots = edgeMap(G, roots, coloring_f<W>(priorities.begin()), -1,
                             sparse_blocked);
    em_t.stop();
    roots = std::move(new_roots);
    rounds++;
  }
  std::cout << "### Total rounds = " << rounds << "\n";
  debug(color_t.next("coloring time"); em_t.next("edge map time"););
  return colors;
}

template <class Graph, class Seq>
inline void verify_coloring(Graph& G, Seq& colors) {
  using W = typename Graph::weight_type;
  size_t n = G.n;
  auto ok = sequence<bool>(n);
  parallel_for(0, n, [&](size_t i) {
    uintE src_color = colors[i];
    auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      uintE ngh_color = colors[ngh];
      return src_color == ngh_color;
    };
    size_t ct = G.get_vertex(i).out_neighbors().count(pred);
    ok[i] = (ct > 0);
  });
  auto im_f = [&](size_t i) { return (size_t)ok[i]; };
  auto im = parlay::delayed_seq<size_t>(n, im_f);
  size_t ct = parlay::reduce(im);
  std::cout << "ct = " << ct << "\n";
  if (ct > 0) {
    std::cout << "Invalid coloring"
              << "\n";
  } else {
    std::cout << "Valid coloring"
              << "\n";
  }
}
}  // namespace gbbs
