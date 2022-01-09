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
#include "gbbs/helpers/resizable_table.h"

// The include below is currently not useful, as the majority of out/in-degree
// one vertices are removed in a single round of peeling (so multiple rounds are
// not necessary, at least on the graphs we tested on).
//#include "third_party/gbbs/src/chains.h"

namespace gbbs {
constexpr size_t TOP_BIT = ((size_t)LONG_MAX) + 1;
constexpr size_t VAL_MASK = LONG_MAX;
using label_type = size_t;

using K = uintE;
using V = uintE;
using T = std::tuple<K, V>;

// hash32 is sufficient
struct hash_kv {
  uint64_t operator()(const K& k) { return parlay::hash64(k); }
};

template <class W, class Seq, class Tab>
struct Search_F {
  Tab& tab;
  Seq& labels;
  bool* bits;
  Search_F(Tab& _tab, Seq& _labels, bool* _bits)
      : tab(_tab), labels(_labels), bits(_bits) {}
  inline bool update(const uintE& s, const uintE& d, const W& wgh) {
    return updateAtomic(s, d, wgh);
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    // check that d is in the same subproblem as s.
    if (labels[s] == labels[d]) {
      bool labels_changed = false;
      auto s_iter = tab.get_iter(s);
      // iterate through s's labels.
      if (s_iter.init()) {
        {
          auto table_entry = s_iter.next();
          uintE s_label = std::get<1>(table_entry);
          labels_changed |= tab.insert(std::make_tuple(d, s_label));
        }

        while (s_iter.has_next()) {
          auto table_entry = s_iter.next();
          uintE s_label = std::get<1>(table_entry);
          labels_changed |= tab.insert(std::make_tuple(d, s_label));
        }
      }
      if (labels_changed) {
        // d should be included in next frontier;
        // CAS to make sure only one ngh from this frontier adds it.
        return gbbs::atomic_compare_and_swap(&bits[d], false, true);
      }
    }
    return false;
  }

  inline bool cond(uintE d) {
    // only visit vertices that are not already in an
    // StronglyConnectedComponents.
    return !(labels[d] & TOP_BIT);
  }
};

template <class W, class Seq, class Tab>
inline Search_F<W, Seq, Tab> make_search_f(Tab& tab, Seq& labels, bool* bits) {
  return Search_F<W, Seq, Tab>(tab, labels, bits);
}

template <class Graph, class Seq, class VS>
inline gbbs::resizable_table<K, V, hash_kv> multi_search(Graph& GA, Seq& labels,
                                                         bool* bits,
                                                         VS& frontier,
                                                         size_t label_start,
                                                         const flags fl = 0) {
  using W = typename Graph::weight_type;

  // table stores (vertex, label) pairs
  T empty = std::make_tuple(UINT_E_MAX, UINT_E_MAX);
  size_t backing_size = 1 << parlay::log2_up(frontier.size() * 2);
  auto table =
      gbbs::resizable_table<K, V, hash_kv>(backing_size, empty, hash_kv());

  frontier.toSparse();
  parallel_for(0, frontier.size(), kDefaultGranularity, [&](size_t i) {
    uintE v = frontier.s[i];
    // each center initially just stores itself.
    table.insert(std::make_tuple(v, label_start + i));
  });
  table.update_nelms();

  size_t rd = 0;
  while (!frontier.isEmpty()) {
    frontier.toSparse();

    auto im_f = [&](size_t i) {
      uintE v = frontier.s[i];
      size_t n_labels = table.num_appearances(v);
      auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
        // can only add labels to vertices in our subproblem
        return labels[ngh] == labels[v];
      };
      size_t effective_degree =
          (fl & in_edges) ? GA.get_vertex(v).in_neighbors().count(pred)
                          : GA.get_vertex(v).out_neighbors().count(pred);
      return effective_degree * n_labels;
    };
    auto im = parlay::delayed_seq<size_t>(frontier.size(), im_f);

    size_t sum = parlay::reduce(im);
    table.maybe_resize(sum);

    parallel_for(0, frontier.size(), kDefaultGranularity, [&](size_t i) {
      uintE v = frontier.s[i];
      bits[v] = 0;  // reset flag
    });

    vertexSubset output = edgeMap(
        GA, frontier, make_search_f<W>(table, labels, bits), -1, fl | no_dense);
    table.update_nelms();
    frontier = std::move(output);
    rd++;
  }
  return table;
}

template <class V, class L>
struct First_Search {
  V& visited;
  L& labels;
  First_Search(V& _visited, L& _labels) : visited(_visited), labels(_labels) {}
  inline bool update(uintE s, uintE d) {
    visited[d] = true;
    return true;
  }
  inline bool updateAtomic(uintE s, uintE d) {
    return gbbs::atomic_compare_and_swap(&visited[d], false, true);
  }
  inline bool cond(uintE d) { return !(labels[d] & TOP_BIT) && !visited[d]; }
};

template <class V, class L>
inline First_Search<V, L> make_first_search(V& visited, L& labels) {
  return First_Search<V, L>(visited, labels);
}

template <class Graph, class L>
inline sequence<bool> first_search(Graph& GA, L& labels, uintE start,
                                   size_t label_start, const flags fl = 0) {
  using W = typename Graph::weight_type;
  size_t n = GA.n;

  auto Flags = sequence<bool>(n, false);
  Flags[start] = true;

  auto frontier = vertexSubset(n, start);
  size_t rd = 0;
  while (!frontier.isEmpty()) {
    vertexSubset output = edgeMap(
        GA, frontier, wrap_em_f<W>(make_first_search(Flags, labels)), -1, fl);
    frontier = std::move(output);
    rd++;
  }
  return Flags;
}

template <class Graph>
inline sequence<label_type> StronglyConnectedComponents(Graph& GA,
                                                        double beta = 1.5) {
  timer initt;
  initt.start();
  size_t n = GA.n;
  // Everyone's initial label is 0 (all in the same subproblem)
  auto labels =
      sequence<label_type>::from_function(n, [](size_t) { return 0; });
  auto ba = sequence<bool>(n, false);
  auto bits = ba.begin();

  auto v_im = parlay::delayed_seq<uintE>(n, [](size_t i) { return i; });
  auto zero = parlay::filter(v_im, [&](size_t i) {
    return (GA.get_vertex(i).out_degree() == 0) ||
           (GA.get_vertex(i).in_degree() == 0);
  });
  auto NZ = parlay::filter(v_im, [&](size_t i) {
    return (GA.get_vertex(i).out_degree() > 0) &&
           (GA.get_vertex(i).in_degree() > 0);
  });

  auto P = parlay::random_shuffle(NZ);
  std::cout << "Filtered: " << zero.size()
            << " vertices. Num remaining = " << P.size() << "\n";

  // Assign labels from [0...zero.size())
  parallel_for(0, zero.size(), kDefaultGranularity,
               [&](size_t i) { labels[zero[i]] = 1 + (i | TOP_BIT); });

  size_t step_size = 1, cur_offset = 0, finished = 0, cur_round = 0;
  double step_multiplier = beta;
  size_t label_offset = zero.size() + 1;

  initt.stop();
  initt.next("init");

  // Run the first search (BFS)
  {
    timer hd;
    hd.start();
    auto deg_im_f = [&](size_t i) {
      return std::make_tuple(i, GA.get_vertex(i).out_degree());
    };
    auto deg_im = parlay::delayed_seq<std::tuple<uintE, uintE>>(n, deg_im_f);
    auto red_f = [](const std::tuple<uintE, uintE>& l,
                    const std::tuple<uintE, uintE>& r) {
      return (std::get<1>(l) > std::get<1>(r)) ? l : r;
    };
    auto id = std::make_tuple<uintE, uintE>(0, 0);
    auto monoid = parlay::make_monoid(red_f, id);
    std::tuple<uintE, uintE> sAndD = parlay::reduce(deg_im, monoid);
    uintE start = std::get<0>(sAndD);

    if (!(labels[start] & TOP_BIT)) {
      auto in_visits = first_search(GA, labels, start, label_offset, in_edges);
      auto out_visits = first_search(GA, labels, start, label_offset);
      size_t label = label_offset;
      parallel_for(0, n, [&](size_t i) {
        bool inv = in_visits[i];
        bool outv = out_visits[i];
        if (inv && outv) {
          labels[i] = label | TOP_BIT;  // In the Big SCC
        } else if (inv || outv) {
          labels[i] = label;  // Reachabel from the Big SCC, but not in it.
        }
      });
      label_offset += 1;
      hd.stop();
      hd.next("big scc time");
    }
  }

  auto Q = parlay::filter(P, [&](uintE v) { return !(labels[v] & TOP_BIT); });
  std::cout << "After first round, Q = " << Q.size()
            << " vertices remain. Total done = " << (n - Q.size()) << "\n";

  while (finished < Q.size()) {
    timer rt;
    rt.start();
    // Running searches between P[cur_offset, end)
    size_t end = std::min(cur_offset + step_size, Q.size());
    size_t vs_size = end - cur_offset;
    finished += vs_size;
    step_size = ceil(step_size * step_multiplier);
    cur_round++;
    size_t round_offset = cur_offset;
    cur_offset += vs_size;

    auto centers_pre_filter = parlay::delayed_seq<uintE>(
        vs_size, [&](size_t i) { return Q[round_offset + i]; });
    auto centers = parlay::filter(
        centers_pre_filter, [&](uintE v) { return !(labels[v] & TOP_BIT); });

    std::cout << "round = " << cur_round << " n_centers = " << centers.size()
              << " originally was " << vs_size
              << " centers. Total vertices remaining = "
              << (Q.size() - finished) << "\n";

    if (centers.size() == 0) continue;

    size_t cur_label_offset = label_offset;
    label_offset += centers.size();

    if (cur_round == 1) {
      timer ft;
      ft.start();
      uintE start = centers[0];
      auto in_visits =
          first_search(GA, labels, start, cur_label_offset, in_edges);
      auto out_visits = first_search(GA, labels, start, cur_label_offset);

      size_t label = cur_label_offset;

      parallel_for(0, n, [&](size_t i) {
        bool inv = in_visits[i];
        bool outv = out_visits[i];
        if (inv && outv) {
          labels[i] = label | TOP_BIT;
        } else if (inv || outv) {
          labels[i] = label;
        }
      });
      ft.stop();
      ft.next("first round time");
      continue;
    }

    timer ins;
    ins.start();
    auto centers_2 = centers;
    auto in_f = vertexSubset(n, std::move(centers));
    auto in_table =
        multi_search(GA, labels, bits, in_f, cur_label_offset, in_edges);
    std::cout << "Finished in search"
              << "\n";
    ins.stop();
    ins.next("insearch time");

    timer outs;
    outs.start();
    auto out_f = vertexSubset(n, std::move(centers_2));
    auto out_table = multi_search(GA, labels, bits, out_f, cur_label_offset);
    std::cout << "in_table, m = " << in_table.m << " ne = " << in_table.ne
              << "\n";
    std::cout << "out_table, m = " << out_table.m << " ne = " << out_table.ne
              << "\n";
    outs.stop();
    outs.next("outsearch time");

    auto& smaller_t = (in_table.m <= out_table.m) ? in_table : out_table;
    auto& larger_t = (in_table.m > out_table.m) ? in_table : out_table;

    // intersect the tables
    auto map_f = [&](const std::tuple<K, V>& kev) {
      uintE v = std::get<0>(kev);
      size_t label = std::get<1>(kev);
      if (larger_t.contains(v, label)) {
        // in 'label' scc
        // Max visitor from this StronglyConnectedComponents acquires it.
        gbbs::write_max(&labels[v], label | TOP_BIT);
      } else {
        gbbs::write_max(&labels[v], label);
      }
    };
    smaller_t.map(map_f);

    // set the subproblems
    auto sp_map = [&](const std::tuple<K, V>& kev) {
      uintE v = std::get<0>(kev);
      size_t label = std::get<1>(kev);
      // note that if v is already in an StronglyConnectedComponents (from (1)),
      // the gbbs::write_max will
      // read, compare and fail, as the top bit is already set.
      gbbs::write_max(&labels[v], label);
    };
    larger_t.map(sp_map);

    rt.stop();
    rt.next("Round time");
  }

  parallel_for(0, labels.size(),
               [&](size_t i) { labels[i] = (labels[i] & VAL_MASK) - 1; });
  return labels;
}

template <class Seq>
inline size_t num_done(Seq& labels) {
  auto im_f = [&](size_t i) { return ((size_t)((labels[i] & TOP_BIT) > 0)); };
  auto im = parlay::delayed_seq<size_t>(labels.size(), im_f);

  return parlay::reduce(im);
}

template <class Seq>
inline size_t num_scc(Seq& labels) {
  size_t n = labels.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  parallel_for(0, n, kDefaultGranularity, [&](size_t i) {
    size_t label = labels[i] & VAL_MASK;
    if (!flags[label]) {
      flags[label] = 1;
    }
  });
  parlay::scan_inplace(flags);
  size_t n_scc = flags[n];
  std::cout << "n_scc = " << flags[n] << "\n";
  return n_scc;
}

template <class Seq>
inline void scc_stats(Seq& labels) {
  size_t n = labels.size();
  auto flags =
      sequence<uintE>::from_function(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    size_t label = labels[i] & VAL_MASK;
    flags[label]++;
  }
  size_t maxv = parlay::reduce_max(flags);
  std::cout << "Largest StronglyConnectedComponents has " << maxv << " vertices"
            << "\n";
  for (size_t i = 0; i < n; i++) {
    if (flags[i] == maxv) {
      std::cout << "max flag = " << i << "\n";
    }
  }
}
}  // namespace gbbs
