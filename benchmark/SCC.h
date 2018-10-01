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

#include "lib/dyn_arr.h"
#include "lib/random_shuffle.h"
#include "lib/resizable_table.h"
#include "lib/sparse_table.h"

#include "chains.h"
#include "ligra.h"

constexpr size_t TOP_BIT = ((size_t)LONG_MAX) + 1;
constexpr size_t VAL_MASK = LONG_MAX;

using K = uintE;
using V = uintE;
using T = tuple<K, V>;

using label_type = size_t;

struct hash_kv {
  uint64_t operator()(const K& k) { return pbbs::hash64(k); }
};

template <class W, class Seq, class Tab>
struct Search_F {
  Seq& labels;
  bool* bits;
  Tab& tab;
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
        auto table_entry = s_iter.next();
        uintE s_label = get<1>(table_entry);
        labels_changed |= tab.insert(make_tuple(d, s_label));

        while (s_iter.has_next()) {
          auto table_entry = s_iter.next();
          uintE s_label = get<1>(table_entry);
          labels_changed |= tab.insert(make_tuple(d, s_label));
        }
      }
      if (labels_changed) {
        // d should be included in next frontier;
        // CAS to make sure only one ngh from this frontier adds it.
        return CAS(&bits[d], false, true);
      }
    }
    return false;
  }

  inline bool cond(uintE d) {
    // only visit vertices that are not already in an SCC.
    return !(labels[d] & TOP_BIT);
  }
};

template <class W, class Seq, class Tab>
auto make_search_f(Tab& tab, Seq& labels, bool* bits) {
  return Search_F<W, Seq, Tab>(tab, labels, bits);
}

template <template <class W> class vertex, class W, class Seq, class VS>
auto multi_search(graph<vertex<W>>& GA, Seq& labels, bool* bits, VS& frontier,
                  size_t label_start, const flags fl = 0) {
  size_t n = GA.n;
  // table stores (vertex, label) pairs
  T empty = make_tuple(UINT_E_MAX, UINT_E_MAX);
  size_t backing_size = 1 << pbbs::log2_up(frontier.size() * 2);
  auto table_backing = array_imap<T>(backing_size);
  auto table = resizable_table<K, V, hash_kv>(backing_size, empty, hash_kv(),
                                              table_backing.get_array(), true);
  frontier.toSparse();
  parallel_for_bc(i, 0, frontier.size(),
                  (frontier.size() > pbbs::kSequentialForThreshold), {
                    uintE v = frontier.s[i];
                    table.insert(make_tuple(v, label_start + i));
                  });
  table.update_nelms();

  size_t rd = 0;
  while (!frontier.isEmpty()) {
    frontier.toSparse();

    auto im = make_in_imap<size_t>(frontier.size(), [&](size_t i) {
      uintE v = frontier.s[i];
      uintE v_lab = labels[v];
      size_t n_labels = table.num_appearances(v);
      auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
        // can only add labels to vertices in our subproblem
        return labels[ngh] == labels[v];
      };
      size_t effective_degree = (fl & in_edges) ? GA.V[v].countInNgh(v, pred)
                                                : GA.V[v].countOutNgh(v, pred);
      return effective_degree * n_labels;
    });

    size_t sum = pbbs::reduce_add(im);
    table.maybe_resize(sum);

    parallel_for_bc(i, 0, frontier.size(),
                    (frontier.size() > pbbs::kSequentialForThreshold), {
                      uintE v = frontier.s[i];
                      bits[v] = 0;  // reset flag
                    });

    vertexSubset output = edgeMap(
        GA, frontier, make_search_f<W>(table, labels, bits), -1, fl | no_dense);
    table.update_nelms();
    frontier.del();
    frontier = output;
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
    return CAS(&visited[d], false, true);
  }
  inline bool cond(uintE d) { return !(labels[d] & TOP_BIT) && !visited[d]; }
};

template <class V, class L>
auto make_first_search(V& visited, L& labels) {
  return First_Search<V, L>(visited, labels);
}

template <template <class W> class vertex, class W, class L>
auto first_search(graph<vertex<W>>& GA, L& labels, uintE start,
                  size_t label_start, const flags fl = 0) {
  size_t n = GA.n;

  auto Flags = array_imap<bool>(n, false);
  Flags[start] = true;

  auto frontier = vertexSubset(n, start);
  size_t rd = 0;
  while (!frontier.isEmpty()) {
    vertexSubset output = edgeMap(
        GA, frontier, wrap_em_f<W>(make_first_search(Flags, labels)), -1, fl);
    frontier.del();
    frontier = output;
    rd++;
  }
  return Flags.get_array();
}

template <class vertex>
auto SCC(graph<vertex>& GA, double beta = 1.1) {
  timer initt;
  initt.start();
  size_t n = GA.n;
  // Everyone's initial label is 0 (all in the same subproblem)
  auto labels = array_imap<label_type>(n, [](size_t) { return 0; });
  auto ba = array_imap<bool>(n, false);
  auto bits = ba.get_array();

  auto v_im = make_in_imap<uintE>(n, [](size_t i) { return i; });
  auto zero_pred = [&](size_t i) {
    return (GA.V[i].getOutDegree() == 0) || (GA.V[i].getInDegree() == 0);
  };
  auto not_zero_pred = [&](size_t i) {
    return (GA.V[i].getOutDegree() > 0) && (GA.V[i].getInDegree() > 0);
  };
  auto zero = pbbs::filter(v_im, zero_pred);
  auto P = pbbs::filter(v_im, not_zero_pred);
  pbbs::random_shuffle(P);
  cout << "Filtered: " << zero.size()
       << " vertices. Num remaining = " << P.size() << endl;

  // Assign labels from [0...zero.size())
  parallel_for_bc(i, 0, zero.size(),
                  (zero.size() > pbbs::kSequentialForThreshold),
                  { labels[zero[i]] = i | TOP_BIT; });

  size_t step_size = 1, cur_offset = 0, finished = 0, cur_round = 0;
  double step_multiplier = beta;
  size_t label_offset = zero.size() + 1;

  auto done_im = make_in_imap<size_t>(
      n, [&](size_t i) { return (bool)(labels[i] & TOP_BIT); });
  initt.stop();
  initt.reportTotal("init");

  {
    timer hd;
    hd.start();
    auto deg_im = make_in_imap<tuple<uintE, uintE>>(
        n, [&](size_t i) { return make_tuple(i, GA.V[i].getOutDegree()); });
    tuple<uintE, uintE> sAndD = pbbs::reduce(
        deg_im, [](const tuple<uintE, uintE>& l, const tuple<uintE, uintE>& r) {
          return (get<1>(l) > get<1>(r)) ? l : r;
        });
    uintE start = get<0>(sAndD);
    if (!(labels[start] & TOP_BIT)) {
      auto in_visits = first_search(GA, labels, start, label_offset, in_edges);
      auto out_visits = first_search(GA, labels, start, label_offset);
      size_t label = label_offset;
      parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
        bool inv = in_visits[i];
        bool outv = out_visits[i];
        if (inv && outv) {
          labels[i] = label | TOP_BIT;
        } else if (inv || outv) {
          labels[i] = label;
        }
      });
      free(in_visits);
      free(out_visits);
      label_offset += 1;
      hd.stop();
      hd.reportTotal("big scc time");
    }
  }

  auto Q = pbbs::filter(P, [&](uintE v) { return !(labels[v] & TOP_BIT); });
  cout << "After first round, Q = " << Q.size()
       << " vertices remain. Total done = " << (n - Q.size()) << endl;

  while (finished < Q.size()) {
    timer rt;
    rt.start();
    // Running searches between P[cur_offset, end)
    size_t end = min(cur_offset + step_size, Q.size());
    size_t vs_size = end - cur_offset;
    finished += vs_size;
    step_size = ceil(step_size * step_multiplier);
    cur_round++;
    size_t round_offset = cur_offset;
    cur_offset += vs_size;

    auto centers_pre_filter = make_in_imap<uintE>(
        vs_size, [&](size_t i) { return Q[round_offset + i]; });
    auto centers = pbbs::filter(
        centers_pre_filter, [&](uintE v) { return !(labels[v] & TOP_BIT); });

    cout << "round = " << cur_round << " n_centers = " << centers.size()
         << " originally was " << vs_size
         << " centers. Total vertices remaining = " << (Q.size() - finished)
         << endl;

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

      parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
        bool inv = in_visits[i];
        bool outv = out_visits[i];
        if (inv && outv) {
          labels[i] = label | TOP_BIT;
        } else if (inv || outv) {
          labels[i] = label;
        }
      });
      free(in_visits);
      free(out_visits);
      ft.stop();
      ft.reportTotal("first round time");
      continue;
    }

    auto centers_2 = centers.copy();
    auto in_f = vertexSubset(n, centers.size(), centers.get_array());
    auto in_table =
        multi_search(GA, labels, bits, in_f, cur_label_offset, in_edges);
    cout << "Finished in search" << endl;

    auto out_f = vertexSubset(n, centers_2.size(), centers_2.get_array());
    auto out_table = multi_search(GA, labels, bits, out_f, cur_label_offset);
    cout << "in_table, m = " << in_table.m << " ne = " << in_table.ne << endl;
    cout << "out_table, m = " << out_table.m << " ne = " << out_table.ne
         << endl;

    auto& smaller_t = (in_table.m <= out_table.m) ? in_table : out_table;
    auto& larger_t = (in_table.m > out_table.m) ? in_table : out_table;

    // intersect the tables
    auto map_f = [&](const tuple<K, V>& kev) {
      uintE v = get<0>(kev);
      size_t label = get<1>(kev);
      if (larger_t.contains(v, label)) {
        // in 'label' scc
        // Max visitor from this SCC acquires it.
        writeMax(&labels[v], label | TOP_BIT);
      } else {
        writeMax(&labels[v], label);
      }
    };
    smaller_t.map(map_f);

    // set the subproblems
    auto sp_map = [&](const tuple<K, V>& kev) {
      uintE v = get<0>(kev);
      size_t label = get<1>(kev);
      // note that if v is already in an SCC (from (1)), the writeMax will
      // read, compare and fail, as the top bit is already set.
      writeMax(&labels[v], label);
    };
    larger_t.map(sp_map);

    in_table.del();
    out_table.del();
    rt.stop();
    rt.reportTotal("Round time");
  }

  return labels;
}

template <class Seq>
size_t num_done(Seq& labels) {
  auto im = make_in_imap<size_t>(labels.size(), [&](size_t i) {
    return ((size_t)((labels[i] & TOP_BIT) > 0));
  });
  return pbbs::reduce_add(im);
}

template <class Seq>
size_t num_scc(Seq& labels) {
  size_t n = labels.size();
  auto flags = array_imap<uintE>(n + 1, [&](size_t i) { return 0; });
  parallel_for_bc(i, 0, n, (n > pbbs::kSequentialForThreshold), {
    if (labels[i] == 0) {
      cout << "unlabeled" << endl;
      exit(0);
    }
    size_t label = labels[i] & VAL_MASK;
    if (!flags[label]) {
      flags[label] = 1;
    }
  });
  pbbs::scan_add(flags, flags);
  size_t n_scc = flags[n];
  cout << "n_scc = " << flags[n] << endl;
  return n_scc;
}

template <class Seq>
void scc_stats(Seq& labels) {
  size_t n = labels.size();
  auto flags = array_imap<uintE>(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    size_t label = labels[i] & VAL_MASK;
    flags[label]++;
  }
  size_t maxv = pbbs::reduce_max(flags);
  cout << "Largest SCC has " << maxv << " vertices" << endl;
  for (size_t i = 0; i < n; i++) {
    if (flags[i] == maxv) {
      cout << "max flag = " << i << endl;
    }
  }
}
