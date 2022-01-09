#pragma once

#include "benchmarks/Connectivity/common.h"
#include "gbbs/bridge.h"
#include "jayanti.h"

#include <mutex>

namespace gbbs {
namespace find_variants {
inline uintE find_naive(uintE i, sequence<parent>& parents) {
  while (i != parents[i]) {
    i = parents[i];
  }
  return i;
}

inline uintE find_compress(uintE i, sequence<parent>& parents) {
  parent j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
  } while (parents[j] != j);
  parent tmp;
  while ((tmp = parents[i]) > j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

inline uintE find_atomic_split(uintE i, sequence<parent>& parents) {
  while (1) {
    parent v = parents[i];
    parent w = parents[v];
    if (v == w) {
      return v;
    } else {
      gbbs::atomic_compare_and_swap(&parents[i], v, w);
      // i = its parents
      i = v;
    }
  }
}

inline uintE find_atomic_halve(uintE i, sequence<parent>& parents) {
  while (1) {
    parent v = parents[i];
    parent w = parents[v];
    if (v == w) {
      return v;
    } else {
      gbbs::atomic_compare_and_swap(&parents[i], (parent)v, (parent)w);
      // i = its grandparent
      i = parents[i];
    }
  }
}
}  // namespace find_variants

namespace splice_variants {

/* Used in Rem-CAS variants for splice */
inline uintE split_atomic_one(uintE i, uintE x, sequence<parent>& parents) {
  parent v = parents[i];
  parent w = parents[v];
  if (v == w)
    return v;
  else {
    gbbs::atomic_compare_and_swap(&parents[i], v, w);
    i = v;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
inline uintE halve_atomic_one(uintE i, uintE x, sequence<parent>& parents) {
  parent v = parents[i];
  parent w = parents[v];
  if (v == w)
    return v;
  else {
    gbbs::atomic_compare_and_swap(&parents[i], v, w);
    i = w;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
inline uintE splice_atomic(uintE u, uintE v, sequence<parent>& parents) {
  parent z = parents[u];
  gbbs::atomic_compare_and_swap(&parents[u], z, parents[v]);
  return z;
}
}  // namespace splice_variants

namespace unite_variants {

template <class Find>
struct Unite {
  Find& find;
  Unite(Find& find) : find(find) {}

  inline uintE operator()(uintE u_orig, uintE v_orig,
                          sequence<parent>& parents) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, parents);
      v = find(v, parents);
      if (u == v)
        break;
      else if (u > v && parents[u] == u &&
               gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
        return u;
      } else if (v > u && parents[v] == v &&
                 gbbs::atomic_compare_and_swap(&parents[v], v, u)) {
        return v;
      }
    }
    return UINT_E_MAX;
  }
};

template <class Splice, class Compress, FindOption find_option>
struct UniteRemLock {
  uintE n;
  sequence<bool> locks;
  Compress& compress;
  Splice& splice;
  UniteRemLock(Compress& compress, Splice& splice, uintE n)
      : n(n), compress(compress), splice(splice) {
    locks = sequence<bool>(n, false);
  }

  inline void fence() { std::atomic_thread_fence(std::memory_order_seq_cst); }

  bool acquire_lock(uintE u) {
    if (!locks[u] && gbbs::atomic_compare_and_swap(&locks[u], false, true)) {
      return true;
    }
    return false;
  }

  void release_lock(uintE u) { locks[u] = false; }

  inline uintE operator()(uintE u_orig, uintE v_orig,
                          sequence<parent>& parents) {
    parent rx = u_orig;
    parent ry = v_orig;
    while (parents[rx] != parents[ry]) {
      /* link from high -> low */
      if (parents[rx] < parents[ry]) std::swap(rx, ry);
      if (rx == parents[rx]) {
        if (acquire_lock(rx)) {
          parent py = parents[ry];
          if (rx == parents[rx] && rx > py) {
            parents[rx] = py;
          }
          release_lock(rx);
        }
        /* link high -> low */
      } else {
        rx = splice(rx, ry, parents);
      }
    }
    if
      constexpr(find_option != find_naive) { /* aka find_none */
        compress(u_orig, parents);
        compress(v_orig, parents);
      }
    return UINT_E_MAX;  // TODO (ret-value)
  }
};

template <class Splice, class Compress, FindOption find_option>
struct UniteRemCAS {
  Compress& compress;
  Splice& splice;
  UniteRemCAS(Compress& compress, Splice& splice)
      : compress(compress), splice(splice) {}

  inline uintE operator()(uintE x, uintE y, sequence<parent>& parents) {
    uintE rx = x;
    uintE ry = y;
    while (parents[rx] != parents[ry]) {
      /* link high -> low */
      parent p_ry = parents[ry];
      parent p_rx = parents[rx];
      if (p_rx < p_ry) {
        std::swap(rx, ry);
        std::swap(p_rx, p_ry);
      }
      if (rx == parents[rx] &&
          gbbs::atomic_compare_and_swap(&parents[rx], rx, p_ry)) {
        if
          constexpr(find_option != find_naive) { /* aka find_none */
            compress(x, parents);
            compress(y, parents);
          }
        return rx;
      } else {
        // failure: locally compress by splicing and try again
        rx = splice(rx, ry, parents);
      }
    }
    return UINT_E_MAX;
  }
};

template <class Find, FindOption find_option>
struct UniteEarly {
  Find& find;
  UniteEarly(Find& find) : find(find) {}
  inline uintE operator()(uintE u, uintE v, sequence<parent>& parents) {
    [[maybe_unused]] uintE u_orig = u, v_orig = v;
    uintE ret = UINT_E_MAX;
    while (u != v) {
      /* link high -> low */
      if (v > u) std::swap(u, v);
      if (parents[u] == u && gbbs::atomic_compare_and_swap(&parents[u], u, v)) {
        ret = u;
        break;
      }
      parent z = parents[u];
      parent w = parents[z];
      gbbs::atomic_compare_and_swap(&parents[u], z, w);
      u = w;
    }
    if
      constexpr(find_option != find_naive) {
        u = find(u_orig, parents); /* force */
        v = find(v_orig, parents); /* force */
      }
    return ret;
  }
};

template <class Find>
struct UniteND {
  Find find;
  sequence<uintE> hooks;
  UniteND(size_t n, Find& find) : find(find) {
    hooks = sequence<uintE>(n, UINT_E_MAX);
  }

  inline uintE operator()(uintE u_orig, uintE v_orig,
                          sequence<parent>& parents) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, parents);
      v = find(v, parents);
      if (u == v) break;
      /* link high -> low */
      if (u < v) std::swap(u, v);
      if (hooks[u] == UINT_E_MAX &&
          gbbs::atomic_compare_and_swap(&hooks[u], UINT_E_MAX, v)) {
        parents[u] = v;
        return u;
      }
    }
    return UINT_E_MAX;
  }
};

}  // namespace unite_variants
}  // namespace gbbs
