#pragma once

#include "gbbs/bridge.h"
#include "pbbslib/seq.h"
#include "jayanti.h"
#include "benchmarks/Connectivity/common.h"

#include <mutex>

namespace gbbs {
namespace find_variants {
  inline uintE find_naive(uintE i, pbbs::sequence<parent>& parents) {
    uintE pathlen = 1;
    while(i != parents[i]) {
      i = parents[i];
      pathlen++;
    }
    report_pathlen(pathlen);
    return i;
  }

  inline uintE find_compress(uintE i, pbbs::sequence<parent>& parents) {
    uintE pathlen = 1;
    parent j = i;
    if (parents[j] == j) return j;
    do {
      j = parents[j];
      pathlen++;
    } while (parents[j] != j);
    parent tmp;
    while ((tmp=parents[i])>j) {
      parents[i] = j; i=tmp;
    }
    report_pathlen(pathlen);
    return j;
  }

  inline uintE find_atomic_split(uintE i, pbbs::sequence<parent>& parents) {
    uintE pathlen = 1;
    while(1) {
      parent v = parents[i];
      parent w = parents[v];
      if (v == w) {
        report_pathlen(pathlen);
        return v;
      }
      else {
        pbbs::atomic_compare_and_swap(&parents[i],v,w);
        // i = its parents
        i = v;
      }
      pathlen++;
    }
  }

  inline uintE find_atomic_halve(uintE i, pbbs::sequence<parent>& parents) {
    uintE pathlen = 1;
    while(1) {
      parent v = parents[i];
      parent w = parents[v];
      if(v == w) {
        report_pathlen(pathlen);
        return v;
      } else {
        pbbs::atomic_compare_and_swap(&parents[i],(parent)v,(parent)w);
        // i = its grandparent
        i = parents[i];
      }
      pathlen++;
    }
  }
} // namespace find_variants


namespace splice_variants {

  /* Used in Rem-CAS variants for splice */
  inline uintE split_atomic_one(uintE i, uintE x, pbbs::sequence<parent>& parents) {
    parent v = parents[i];
    parent w = parents[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parents[i],v,w);
      i = v;
      return i;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE halve_atomic_one(uintE i, uintE x, pbbs::sequence<parent>& parents) {
    parent v = parents[i];
    parent w = parents[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parents[i],v,w);
      i = w;
      return i;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE splice_atomic(uintE u, uintE v, pbbs::sequence<parent>& parents) {
    parent z = parents[u];
    pbbs::atomic_compare_and_swap(&parents[u], z, parents[v]);
    return z;
  }
} // namespace splice_variants


namespace unite_variants {

  template <class Find>
  struct Unite {
    Find& find;
    Unite(Find& find) : find(find) { }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& parents) {
      parent u = u_orig;
      parent v = v_orig;
      while(1) {
        u = find(u,parents);
        v = find(v,parents);
        if(u == v) break;
        else if (u > v && parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) {
          break;
        }
        else if (v > u && parents[v] == v && pbbs::atomic_compare_and_swap(&parents[v],v,u)) {
          break;
        }
      }
    }
  };

  template <class Splice, class Compress, FindOption find_option>
  struct UniteRemLock {
    uintE n;
    //pbbs::sequence<std::atomic<bool>> locks;
    pbbs::sequence<bool> locks;
    Compress& compress;
    Splice& splice;
    UniteRemLock(Compress& compress, Splice& splice, uintE n) : n(n), compress(compress), splice(splice) {
      locks = pbbs::sequence<bool>(n, false);
//      locks = pbbs::sequence<std::atomic<bool>>(n);
//      std::cout << "initializing locks, n = " << n << std::endl;
//      parallel_for(0, n, [&] (size_t i) {
//        locks[i].store(false);
//      });
//      parallel_for (0, n, [&] (size_t i) {
//        if (locks[i].load()) {
//          std::cout << "Bad memory at i = " << i << std::endl;
//        }
//      });
    }

    inline void fence() {
      std::atomic_thread_fence(std::memory_order_seq_cst);
    }

    bool acquire_lock(uintE u) {
      if (!locks[u] && pbbs::atomic_compare_and_swap(&locks[u], false, true)) {
        return true;
      }
      return false;
//      bool expected = false;
//      if (locks[u].compare_exchange_strong(expected, true)) {
//        return true;
//      }
//      return false;
    }

    void release_lock(uintE u) {
      locks[u] = false;
//      bool expected = true;
//      bool ret = locks[u].compare_exchange_strong(expected, false, std::memory_order_seq_cst);
//      assert(ret);
//      assert(expected);
      //locks[u] = false;
//      bool ret = pbbs::atomic_compare_and_swap<bool>(&locks[u], true, false);
//      assert(ret);
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& parents) {
      parent rx = u_orig;
      parent ry = v_orig;
      uintE pathlen = 1;
      uintE lock_attempts = 0;
      while (parents[rx] != parents[ry]) {
        /* link from high -> low */
        if (parents[rx] < parents[ry]) std::swap(rx,ry);
        if (rx == parents[rx]) {
          if (acquire_lock(rx)) {
            parent py = parents[ry];
            if (rx == parents[rx] && rx > py) {
              parents[rx] = py;
            }
            release_lock(rx);
          } else {
            lock_attempts++;
            if (lock_attempts % 10000000 == 0) {
              std::cout << "Spinning trying to acquire: " << rx << " lock says: " << locks[rx] << std::endl;
            }
          }
          /* link high -> low */
        } else {
          rx = splice(rx, ry, parents);
        }
        pathlen++;
      }
      if constexpr (find_option != find_naive) { /* aka find_none */
        compress(u_orig, parents);
        compress(v_orig, parents);
      }
      report_pathlen(pathlen);
      return;
    }
  };

  template <class Splice, class Compress, FindOption find_option>
  struct UniteRemCAS {
    Compress& compress;
    Splice& splice;
    UniteRemCAS(Compress& compress, Splice& splice) : compress(compress), splice(splice) { }

    inline void operator()(uintE x, uintE y, pbbs::sequence<parent>& parents) {
      uintE rx = x; uintE ry = y;
      uintE pathlen = 1;
      while (parents[rx] != parents[ry]) {
        /* link high -> low */
        parent p_ry = parents[ry];
        parent p_rx = parents[rx];
        if (p_rx < p_ry) {
          std::swap(rx, ry);
          std::swap(p_rx, p_ry);
        }
        if (rx == parents[rx] && pbbs::atomic_compare_and_swap(&parents[rx], rx, p_ry)) {
          // success
          if constexpr (find_option != find_naive) { /* aka find_none */
            compress(x, parents);
            compress(y, parents);
          }
          break;
        } else {
          // failure: locally compress by splicing and try again
          rx = splice(rx, ry, parents);
        }
        pathlen++;
      }
      report_pathlen(pathlen);
      return;
    }
  };

  template <class Find, FindOption find_option>
  struct UniteEarly {
    Find& find;
    UniteEarly(Find& find) : find(find) {}
    inline void operator()(uintE u, uintE v, pbbs::sequence<parent>& parents) {
      [[maybe_unused]] uintE u_orig = u, v_orig = v;
      int pathlen = 0;
      while(u != v) {
        /* link high -> low */
        if(v > u) std::swap(u,v);
        if (parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) {
          break;
        }
        parent z = parents[u];
        parent w = parents[z];
        pbbs::atomic_compare_and_swap(&parents[u],z,w);
        u = w;
        pathlen++;
      }
      if constexpr (find_option != find_naive) {
        u = find(u_orig, parents); /* force */
        v = find(v_orig, parents); /* force */
      }
      report_pathlen(pathlen);
    }
  };

  template <class Find>
  struct UniteND {
    Find find;
    pbbs::sequence<uintE> hooks;
    UniteND(size_t n, Find& find) : find(find) {
      hooks = pbbs::sequence<uintE>(n, UINT_E_MAX);
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& parents) {
      parent u = u_orig;
      parent v = v_orig;
      while(1) {
        u = find(u,parents);
        v = find(v,parents);
        if(u == v) break;
        /* link high -> low */
        if(u < v) std::swap(u,v);
        if (hooks[u] == UINT_E_MAX && pbbs::atomic_compare_and_swap(&hooks[u], UINT_E_MAX,v)) {
          parents[u] = v;
          break;
        }
      }
    }
  };

  /* Add unite-by-size (lock-based?) */

}  // namespace unite_variants

}  // namespace gbbs
