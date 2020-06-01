#pragma once

#include "benchmarks/SpanningForest/common.h"

namespace lt {

  bool lt_less(uintE u, uintE v) {
    if (u == UINT_E_MAX) {
      return (v != UINT_E_MAX);
    } else if (v == UINT_E_MAX) {
      return false;
    }
    return u < v;
  }

  bool lt_greater(uintE u, uintE v) {
    if (u == UINT_E_MAX) {
      return (v == UINT_E_MAX);
    } else if (v == UINT_E_MAX) {
      return true;
    }
    return u > v;
  }

  uintE lt_min(uintE u, uintE v) {
    if (u == UINT_E_MAX) {
      return u;
    } else if (v == UINT_E_MAX) {
      return v;
    }
    return std::min(u, v);
  }

  uintE lt_max(uintE u, uintE v) {
    if (u == UINT_E_MAX) {
      return v;
    } else if (v == UINT_E_MAX) {
      return u;
    }
    return std::max(u, v);
  }

  struct message {
    uintE component;
    uintE padding;
    std::pair<uintE, uintE> original_edge;

    message(uintE component, std::pair<uintE, uintE> edge) : component(component), original_edge(edge) {}
    message() {}
  };

  bool message_less(message l, message r) {
    return lt_less(l.component, r.component);
  }

  namespace primitives {

    // For each edge e, request e.v.p from e.v and e.w.p from e.w; send the minimum of
    // the received vertices to the maximum of the received vertices.
    bool parent_connect(uintE u, uintE v, pbbs::sequence<parent>& P, pbbs::sequence<message>& messages) {
      uintE p_u = P[u];
      uintE p_v = P[v];
      auto min_v = lt_min(p_u, p_v);
      auto max_v = lt_max(p_u, p_v);
      if (min_v != max_v) {
        return pbbs::write_min<message>(&messages[max_v], message(min_v, std::make_pair(u,v)), message_less);
      }
      return false;
    }

    // For each edge e, request e.v.p from e.v and e.w.p from e.w; let the received
    // values be x and y, respectively; if y < x then send y to v and to x
    // else send x to w and to y.
    bool extended_connect(uintE v, uintE w, pbbs::sequence<parent>& P, pbbs::sequence<message>& messages) {
      uintE x = P[v];
      uintE y = P[w];
      bool updated = false;
      if (lt_less(y, x)) { /* send y to {v, x}*/
        updated |= pbbs::write_min(&messages[v], message(y, std::make_pair(v, w)), message_less);
        updated |= pbbs::write_min(&messages[x], message(y, std::make_pair(v, w)), message_less);
      } else if (lt_greater(y, x)) { /* send x to {y, w} */
        updated |= pbbs::write_min(&messages[y], message(x, std::make_pair(v, w)), message_less);
        updated |= pbbs::write_min(&messages[w], message(x, std::make_pair(v, w)), message_less);
      }
      return updated;
    }

    void root_update(uintE u, pbbs::sequence<parent>& P, pbbs::sequence<message>& messages, pbbs::sequence<edge>& Edges) {
      // find root
      uintE p_u = P[u];
      if (u != p_u) {
        return;
      }
      // otherwise u is a root. replace parent with min of self, and messages[u]
      P[u] = lt_min(p_u, messages[u].component);
      if (p_u != P[u]) {
        assert(Edges[p_u] == empty_edge);
        assert(messages[u].original_edge != empty_edge);
        Edges[p_u] = messages[u].original_edge;
      }
    }

    void shortcut(uintE u, pbbs::sequence<parent>& P) {
      uintE p_u = P[u];
      if (p_u != spanning_forest::largest_comp && p_u != P[p_u]) {
        P[u] = P[p_u];
      }
    }

    void root_shortcut(uintE u, pbbs::sequence<parent>& P) {
      uintE p_u = P[u];
      while (p_u != spanning_forest::largest_comp && p_u != P[p_u]) {
        P[u] = P[p_u];
        p_u = P[u];
      }
    }

    std::tuple<uintE, uintE>
    alter(uintE u, uintE v, pbbs::sequence<parent>& P) {
      uintE x = P[u];
      uintE y = P[v];
      return std::make_tuple(x, y);
    }

  } // namespace primitives

} // namespace lt
