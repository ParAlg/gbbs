#pragma once

#include "benchmarks/Connectivity/common.h"
#include "gbbs/bridge.h"
#include "jayanti.h"

#include <mutex>

namespace gbbs {
namespace find_variants {
inline uintE find_naive(uintE i, sequence<parent>& Parents) {
  uintE pathlen = 1;
  while (i != Parents[i]) {
    i = Parents[i];
    pathlen++;
  }
  report_pathlen(pathlen);
  return i;
}

inline uintE find_compress(uintE i, sequence<parent>& Parents) {
  uintE pathlen = 1;
  parent j = i;
  if (Parents[j] == j) return j;
  do {
    j = Parents[j];
    pathlen++;
  } while (Parents[j] != j);
  parent tmp;
  while ((tmp = Parents[i]) > j) {
    Parents[i] = j;
    i = tmp;
  }
  report_pathlen(pathlen);
  return j;
}

inline uintE find_atomic_split(uintE i, sequence<parent>& Parents) {
  uintE pathlen = 1;
  while (1) {
    parent v = Parents[i];
    parent w = Parents[v];
    if (v == w) {
      report_pathlen(pathlen);
      return v;
    } else {
      gbbs::atomic_compare_and_swap(&Parents[i], v, w);
      // i = its Parents
      i = v;
    }
    pathlen++;
  }
}

inline uintE find_atomic_halve(uintE i, sequence<parent>& Parents) {
  uintE pathlen = 1;
  while (1) {
    parent v = Parents[i];
    parent w = Parents[v];
    if (v == w) {
      report_pathlen(pathlen);
      return v;
    } else {
      gbbs::atomic_compare_and_swap(&Parents[i], (parent)v, (parent)w);
      // i = its grandparent
      i = Parents[i];
    }
    pathlen++;
  }
}
}  // namespace find_variants

namespace splice_variants {

/* Used in Rem-CAS variants for splice */
inline uintE split_atomic_one(uintE i, uintE x, sequence<parent>& Parents) {
  parent v = Parents[i];
  parent w = Parents[v];
  if (v == w)
    return v;
  else {
    gbbs::atomic_compare_and_swap(&Parents[i], v, w);
    i = v;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
inline uintE halve_atomic_one(uintE i, uintE x, sequence<parent>& Parents) {
  parent v = Parents[i];
  parent w = Parents[v];
  if (v == w)
    return v;
  else {
    gbbs::atomic_compare_and_swap(&Parents[i], v, w);
    i = w;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
inline uintE splice_atomic(uintE u, uintE v, sequence<parent>& Parents) {
  parent z = Parents[u];
  gbbs::atomic_compare_and_swap(&Parents[u], z, Parents[v]);
  return z;
}
}  // namespace splice_variants

namespace unite_variants {

template <class Find>
struct Unite {
  Find& find;
  Unite(Find& find) : find(find) {}

  inline void operator()(uintE u_orig, uintE v_orig, sequence<parent>& Parents,
                         sequence<edge>& Edges) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, Parents);
      v = find(v, Parents);
      if (u == v)
        break;
      else if (u > v && Parents[u] == u &&
               gbbs::atomic_compare_and_swap(&Parents[u], u, v)) {
        Edges[u] = std::make_pair(u_orig, v_orig);
        break;
      } else if (v > u && Parents[v] == v &&
                 gbbs::atomic_compare_and_swap(&Parents[v], v, u)) {
        Edges[v] = std::make_pair(u_orig, v_orig);
        break;
      }
    }
  }
};

template <class Compress, class Splice>
struct UniteRemLock {
  uintE n;
  std::mutex* locks;
  Compress& compress;
  Splice& splice;
  UniteRemLock(Compress& compress, Splice& splice, uintE n)
      : n(n), compress(compress), splice(splice) {
    locks = gbbs::new_array<std::mutex>(n);
  }

  ~UniteRemLock() { gbbs::free_array(locks, n); }

  inline void operator()(uintE u_orig, uintE v_orig, sequence<parent>& Parents,
                         sequence<edge>& Edges) {
    parent rx = u_orig;
    parent ry = v_orig;
    uintE pathlen = 1;
    while (Parents[rx] != Parents[ry]) {
      /* link from high -> low */
      if (Parents[rx] < Parents[ry]) std::swap(rx, ry);
      if (rx == Parents[rx]) {
        locks[rx].lock();
        parent py = Parents[ry];
        /* link high -> low */
        if (rx == Parents[rx] && rx > py) {
          Edges[rx] = std::make_pair(u_orig, v_orig);
          Parents[rx] = py;
        }
        locks[rx].unlock();
      } else {
        rx = splice(rx, ry, Parents);
      }
      pathlen++;
    }
    compress(u_orig, Parents);
    compress(v_orig, Parents);
    report_pathlen(pathlen);
  }
};

template <class Splice, class Compress, FindOption find_option>
struct UniteRemCAS {
  Compress& compress;
  Splice& splice;
  UniteRemCAS(Compress& compress, Splice& splice)
      : compress(compress), splice(splice) {}

  inline void operator()(uintE x, uintE y, sequence<parent>& Parents,
                         sequence<edge>& Edges) {
    uintE rx = x;
    uintE ry = y;
    uintE pathlen = 1;
    while (Parents[rx] != Parents[ry]) {
      /* link high -> low */
      parent p_ry = Parents[ry];
      if (Parents[rx] < p_ry) {
        std::swap(rx, ry);
      }
      if (rx == Parents[rx] &&
          gbbs::atomic_compare_and_swap(&Parents[rx], rx, p_ry)) {
        Edges[rx] = std::make_pair(x, y);
        // success
        if
          constexpr(find_option != find_naive) { /* aka find_none */
            compress(x, Parents);
            compress(y, Parents);
          }
        break;
      } else {
        // failure: locally compress by splicing and try again
        rx = splice(rx, ry, Parents);
      }
      pathlen++;
    }
    report_pathlen(pathlen);
  }
};

struct UniteEarly {
  UniteEarly() {}
  inline void operator()(uintE u_orig, uintE v_orig, sequence<parent>& Parents,
                         sequence<edge>& Edges) {
    uintE u = u_orig;
    uintE v = v_orig;
    while (u != v) {
      /* link high -> low */
      if (v > u) std::swap(u, v);
      if (Parents[u] == u && gbbs::atomic_compare_and_swap(&Parents[u], u, v)) {
        Edges[u] = std::make_pair(u_orig, v_orig);
        break;
      }
      parent z = Parents[u];
      parent w = Parents[z];
      gbbs::atomic_compare_and_swap(&Parents[u], z, w);
      u = w;
    }
  }
};

template <class Find>
struct UniteND {
  Find find;
  sequence<uintE> hooks;
  UniteND(size_t n, Find& find) : find(find) {
    hooks = sequence<uintE>(n, UINT_E_MAX);
  }

  inline void operator()(uintE u_orig, uintE v_orig, sequence<parent>& Parents,
                         sequence<edge>& Edges) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, Parents);
      v = find(v, Parents);
      if (u == v) break;
      /* link high -> low */
      if (u < v) std::swap(u, v);
      if (hooks[u] == UINT_E_MAX &&
          gbbs::atomic_compare_and_swap(&hooks[u], UINT_E_MAX, v)) {
        Edges[u] = std::make_pair(u_orig, v_orig);
        Parents[u] = v;
        break;
      }
    }
  }
};

}  // namespace unite_variants
}  // namespace gbbs
