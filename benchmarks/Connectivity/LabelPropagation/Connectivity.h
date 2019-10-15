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

#include "ligra/ligra.h"
#include "pbbslib/random_shuffle.h"
#include "benchmarks/Connectivity/common.h"

namespace labelprop_cc {

  template <class W>
  struct LabelProp_F {
    pbbs::sequence<parent>& components;
    pbbs::sequence<bool>& changed;
    LabelProp_F(pbbs::sequence<parent>& components, pbbs::sequence<bool>& changed) : components(components), changed(changed) {}
    inline bool update(const uintE& s, const uintE& d, const W& w) {
      if (components[s] < components[d]) {
        components[d] = components[s];
        if (!changed[d]) {
          changed[d] = true;
          return 1;
        }
      }
      return 0;
    }
    inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
      if (components[s] < components[d]) {
        pbbs::write_min<uintE>(&components[d], components[s], std::less<uintE>());
        if (!changed[d]) {
          return pbbs::atomic_compare_and_swap(&changed[d], false, true);
        }
      }
      return 0;
    }
    inline bool cond(const uintE& d) { return true; }
  };

  template <class Graph>
  void CC_impl(Graph& G, pbbs::sequence<parent>& components) {
    using W = typename Graph::weight_type;
    size_t n = G.n;
    auto all = pbbs::sequence<bool>(n, true);
    auto vs = vertexSubset(n, n, all.to_array());
    size_t rounds = 0;
    auto changed = pbbs::sequence<bool>(n, false);
    size_t vertices_processed = 0;
    while (!vs.isEmpty()) {
      // std::cout << vs.size() << std::endl;
      vertices_processed += vs.size();
      timer tt; tt.start();
      auto next_vs = edgeMap(G, vs, LabelProp_F<W>(components, changed), -1, dense_forward);
      tt.stop(); tt.reportTotal("edge map time");
      vs.del();
      vs = next_vs;
      vertexMap(vs, [&] (const uintE u) { changed[u] = false; });
      rounds++;
    }
    std::cout << "# LabelProp: ran " << rounds << " many rounds." << std::endl;
    std::cout << "# processed " << vertices_processed << " many vertices" << std::endl;
  }

  template <bool use_permutation, class Graph>
  inline sequence<parent> CC(Graph& G) {
    size_t n = G.n;
    pbbs::sequence<parent> components;
    if constexpr (use_permutation) {
      components = pbbs::random_permutation<uintE>(n);
    } else {
      components = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
    }
    CC_impl(G, components);
    return components;
  }

  template <class Graph>
  struct LPAlgorithm {
    Graph& GA;
    LPAlgorithm(Graph& GA) : GA(GA) {}

    template <SamplingOption sampling_option>
    void compute_components(pbbs::sequence<parent>& parents, uintE frequent_comp = UINT_E_MAX) {
      using W = typename Graph::weight_type;
      size_t n = GA.n;


      auto vs = vertexSubset(n);
//      if constexpr (sampling_option == sample_bfs) { /* provides component */
//        auto in_imap = pbbslib::make_sequence<uintE>(n, [&] (size_t i) {
//          return i;
//        });
//        auto out = pbbs::filter(in_imap, [&] (const uintE& u) {
//          return parents[u] != frequent_comp;
//        });
//        size_t out_size = out.size();
//        vs = vertexSubset(n, out_size, ((std::tuple<uintE, pbbs::empty>*)out.to_array()));
//      } else {
        auto all = pbbs::sequence<bool>(n, true);
        vs = vertexSubset(n, n, all.to_array());
//      }

      size_t rounds = 0;
      auto changed = pbbs::sequence<bool>(n, false);
      size_t vertices_processed = 0;
      while (!vs.isEmpty()) {
        std::cout << "vs size = " << vs.size() << std::endl;
        vertices_processed += vs.size();
        auto next_vs = edgeMap(GA, vs, LabelProp_F<W>(parents, changed), -1, dense_forward);
        vs.del();
        vs = next_vs;
        vertexMap(vs, [&] (const uintE u) { changed[u] = false; });
        rounds++;
      }
      std::cout << "# LabelProp: ran " << rounds << " many rounds." << std::endl;
      std::cout << "# processed " << vertices_processed << " many vertices" << std::endl;
    }
  };



}  // namespace labelprop_cc
