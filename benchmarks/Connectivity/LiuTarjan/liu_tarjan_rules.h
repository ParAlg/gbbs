#pragma once

#include "ligra/maybe.h"

namespace liu_tarjan {

  namespace primitives {
    /* send minimum end to maximum end */
    bool connect(uintE u, uintE v, pbbs::sequence<uintE>& P) {
      uintE min_v = std::min(u, v);
      uintE max_v = std::max(u, v);
      return pbbs::write_min<uintE>(&P[max_v], min_v, std::less<uintE>());
    }

    bool parent_connect(uintE u, uintE v, pbbs::sequence<uintE>& P) {
      uintE p_u = P[u];
      uintE p_v = P[v];
      auto min_v = std::min(p_u, p_v);
      auto max_v = std::max(p_u, p_v);
      return pbbs::write_min<uintE>(&P[max_v], min_v, std::less<uintE>());
    }

    void extended_connect(uintE v, uintE w, pbbs::sequence<uintE>& P) {
      uintE x = P[v];
      uintE y = P[w];
      if (y < x) { /* send y to {v, x}*/
        pbbs::write_min<uintE>(&P[v], y, std::less<uintE>());
        pbbs::write_min<uintE>(&P[x], y, std::less<uintE>());
      } else { /* send x to {y, w} */
        pbbs::write_min<uintE>(&P[y], x, std::less<uintE>());
        pbbs::write_min<uintE>(&P[w], x, std::less<uintE>());
      }
    }

    void simple_update(uintE u, pbbs::sequence<uintE>& P) {
      abort(); // should not be called, since this fn is a noop
    }

    void root_update(uintE u, pbbs::sequence<uintE>& P) {
      /* find root */
      uintE r = u;
      while (P[r] != P[P[r]]) {
        r = P[r];
      }
      if (r != u) {
        P[u] = r; // shortcut
      }
    }

    void shortcut(uintE u, pbbs::sequence<uintE>& P) {
      if (P[u] != P[P[u]]) {
        P[u] = P[P[u]];
      }
    }

    void root_shortcut(uintE u, pbbs::sequence<uintE>& P) {
      while (P[u] != P[P[u]]) {
        P[u] = P[P[u]];
      }
    }

    Maybe<std::tuple<uintE, uintE>>
    alter(uintE u, uintE v, pbbs::sequence<uintE>& P) {
      uintE x = P[u];
      uintE y = P[v];
      if (x == y) {
        return Maybe<std::tuple<uintE, uintE>>();
      } else {
        return Maybe<std::tuple<uintE, uintE>>(std::make_tuple(x, y));
      }
    }
  } // namespace primitives


//  /* Only suitable for edge_array (COO), or for a graph_type permitting
//   * modifications of elements (alter) */
//  template <
//    class Graph,
//    class Connect,
//    class Update,
//    class Shortcut,
//    class Alter,
//    bool  should_edge_filter>
//  pbbs::sequence<uintE> liu_tarjan_algorithm(Graph& G, Connect& connect, Update& update, Shortcut& shortcut, Alter& alter) {
//    size_t n = G.n;
//    using W = typename G::weight_type;
//
//    auto P = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
//    auto next_P = P; // copy
//
//    auto parents_changed = true;
//    while (parents_changed) {
//      parents_changed = false;
//      // Parent-Connect
//      auto map_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//        connect(u, v, P, next_P);
//      };
//      G.map_edges(map_f);
//
//      // Update
//      parallel_for(0, n, [&] (size_t u) {
//        bool changed = update(u, P, next_P);
//        if (changed && !parents_changed) {
//          parents_changed = true;
//        }
//      });
//
//      // Shortcut
//      parallel_for(0, n, [&] (size_t u) {
//        shortcut(u, P);
//      });
//
//      auto alter_f = [&] (const uintE& u, const uintE& v, const W& wgh) {
//        uintE au, av;
//        std::tie(au, av) = alter(u, v, P);
//        return std::make_tuple(u, v, wgh);
//      };
//
//      G.alter_edges(alter_f);
//
//      /* Can optionally filter */
//      if constexpr (should_edge_filter) {
//        auto p = [&] (const uintE& u, const uintE& v, const W& wgh) {
//          return u != v;
//        };
//        G.filter_edges(p);
//      }
//
//      return P;
//    }
//  }

} // namespace liu_tarjan
