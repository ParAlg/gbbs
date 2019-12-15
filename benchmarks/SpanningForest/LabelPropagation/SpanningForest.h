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
#include "benchmarks/SpanningForest/common.h"

namespace labelprop_sf {

  struct label_and_edge {
    uintE label;
    uintE edge_id;

    label_and_edge() {}
    label_and_edge(uintE label, uintE edge_id) : label(label), edge_id(edge_id) {}

    uintE get_label() const { return label; }
    uintE get_edge_id() const { return edge_id; }

    uintE set_label(uintE l) { label = l; }
    uintE set_edge_id(uintE id) { edge_id = id; }
  };

  bool less_label_and_edge(const label_and_edge& l, const label_and_edge& r) {
    return l.get_label() < r.get_label();
  }

  template <class W>
  struct LabelProp_F {
    pbbs::sequence<label_and_edge>& labels;
    pbbs::sequence<bool>& changed;
    LabelProp_F(pbbs::sequence<label_and_edge>& labels, pbbs::sequence<bool>& changed) : labels(labels), changed(changed) {}
    inline bool update(const uintE& s, const uintE& d, const W& w) {
      if (labels[s].get_label() < labels[d].get_label()) {
        labels[d].set_label(labels[s].get_label());
        labels[d].set_edge_id(s); /* use (d--s) edge in the forest */
        if (!changed[d]) {
          changed[d] = true;
          return 1;
        }
      }
      return 0;
    }
    inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) {
      if (labels[s].get_label() < labels[d].get_label()) {
        label_and_edge c(labels[s].get_label(), s);
        pbbs::write_min<label_and_edge>(&labels[d], c, less_label_and_edge);
        if (!changed[d]) {
          return pbbs::atomic_compare_and_swap(&changed[d], false, true);
        }
      }
      return 0;
    }
    inline bool cond(const uintE& d) { return true; }
  };

  template <class Graph>
  pbbs::sequence<edge> SpanningForest_impl(Graph& G, pbbs::sequence<label_and_edge>& labels) {
    using W = typename Graph::weight_type;
    size_t n = G.n;
    auto all = pbbs::sequence<bool>(n, true);
    auto vs = vertexSubset(n, n, all.to_array());
    size_t rounds = 0;
    auto changed = pbbs::sequence<bool>(n, false);
    size_t vertices_processed = 0;
    while (!vs.isEmpty()) {
      vertices_processed += vs.size();
      timer tt; tt.start();
      auto next_vs = edgeMap(G, vs, LabelProp_F<W>(labels, changed), -1, dense_forward);
      tt.stop(); tt.reportTotal("edge map time");
      vs.del();
      vs = next_vs;
      vertexMap(vs, [&] (const uintE u) { changed[u] = false; });
      rounds++;
    }
    std::cout << "# LabelProp: ran " << rounds << " many rounds." << std::endl;
    std::cout << "# processed " << vertices_processed << " many vertices" << std::endl;
    auto all_edges = pbbs::delayed_seq<edge>(n, [&] (uintE i) {
      return std::make_pair(i, labels[i].get_edge_id());
    });
    return pbbs::filter(all_edges, [&] (const edge& e) {
      return e.first != e.second && e.second != UINT_E_MAX;
    });
  }

  template <bool use_permutation, class Graph>
  inline pbbs::sequence<edge> SpanningForest(Graph& G) {
    size_t n = G.n;
    pbbs::sequence<label_and_edge> labels;
    if constexpr (use_permutation) {
      auto perm = pbbs::random_permutation<uintE>(n);
      parallel_for(0, n, [&] (size_t i) {
        labels[i].set_label(perm[i]);
      });
    } else {
      labels = pbbs::sequence<label_and_edge>(n, [&] (size_t i) { return label_and_edge(i, UINT_E_MAX); });
    }
    return SpanningForest_impl(G, labels);
  }

}  // namespace labelprop_sf
