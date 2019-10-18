#pragma once

#include "ligra/maybe.h"
#include "benchmarks/Connectivity/common.h"

namespace lt {

  bool lt_less(uintE u, uintE v) {
    if (u == v) {
      return false;
    } else if (u == UINT_E_MAX) {
      return true;
    } else if (v == UINT_E_MAX) {
      return false;
    }
    return u < v;
  }

  bool lt_greater(uintE u, uintE v) {
    if (u == v) {
      return false;
    } else if (u == UINT_E_MAX) {
      return false;
    } else if (v == UINT_E_MAX) {
      return true;
    }
    return u > v;
  }

  uintE lt_min(uintE u, uintE v) {
    if (u == v) {
      return u;
    } else if (u == UINT_E_MAX) {
      return u;
    } else if (v == UINT_E_MAX) {
      return v;
    }
    return std::min(u, v);
  }

  uintE lt_max(uintE u, uintE v) {
    if (u == v) {
      return u;
    } else if (u == UINT_E_MAX) {
      return v;
    } else if (v == UINT_E_MAX) {
      return u;
    }
    return std::max(u, v);
  }

  namespace primitives {

    /* send minimum end to maximum end */
    bool connect(uintE u, uintE v, pbbs::sequence<parent>& P) {
      uintE min_v = std::min(u, v);
      uintE max_v = std::max(u, v);
      return pbbs::write_min<uintE>(&P[max_v], min_v, std::less<uintE>());
    }

    bool parent_connect(uintE u, uintE v, pbbs::sequence<parent>& P) {
      uintE p_u = P[u];
      uintE p_v = P[v];
      auto min_v = lt_min(p_u, p_v);
      auto max_v = lt_max(p_u, p_v);
      if (min_v != max_v) {
        return pbbs::write_min<uintE>(&P[max_v], min_v, lt_less);
      }
      return false;
    }

    bool extended_connect(uintE v, uintE w, pbbs::sequence<parent>& P) {
      uintE x = P[v];
      uintE y = P[w];
      bool updated = false;
      if (lt_less(y, x)) { /* send y to {v, x}*/
        updated |= pbbs::write_min(&P[v], y, lt_less);
        updated |= pbbs::write_min(&P[x], y, lt_less);
      } else if (lt_greater(y, x)) { /* send x to {y, w} */
        updated |= pbbs::write_min(&P[y], x, lt_less);
        updated |= pbbs::write_min(&P[w], x, lt_less);
      }
      return updated;
    }

    void simple_update(uintE u, pbbs::sequence<parent>& P) {
      abort(); // should not be called, since this fn is a noop
    }

    void root_update(uintE u, pbbs::sequence<parent>& P) {
      /* find root */
      uintE pu = P[u];
      while (pu != largest_comp && pu != P[pu]) {
        pu = P[pu];
      }
      uintE r = pu;

      if (r != u) {
        P[u] = r; // shortcut
      }
    }

    void shortcut(uintE u, pbbs::sequence<parent>& P) {
      uintE pu = P[u];
      if (pu != largest_comp && pu != P[pu]) {
        P[u] = P[pu];
      }
    }

    void root_shortcut(uintE u, pbbs::sequence<parent>& P) {
      uintE pu = P[u];
      while (pu != largest_comp && pu != P[pu]) {
        P[u] = P[pu];
        pu = P[u];
      }
    }

    std::tuple<uintE, uintE>
    alter(uintE u, uintE v, pbbs::sequence<parent>& P) {
      uintE x = P[u];
      uintE y = P[v];
      return std::make_tuple(x, y);
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
//  pbbs::sequence<parent> liu_tarjan_algorithm(Graph& G, Connect& connect, Update& update, Shortcut& shortcut, Alter& alter) {
//    size_t n = G.n;
//    using W = typename G::weight_type;
//
//    auto P = pbbs::sequence<parent>(n, [&] (size_t i) { return i; });
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

} // namespace lt
