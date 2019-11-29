#pragma once

#include "ligra/bridge.h"
#include "pbbslib/seq.h"
#include "jayanti.h"
#include "benchmarks/Connectivity/common.h"

#include <mutex>

namespace find_variants {
  inline uintE find_naive(uintE i, pbbs::sequence<parent>& parents) {
#ifdef REPORT_PATH_LENGTHS
    uintE pathlen = 1;
#endif
    while(i != parents[i]) {
      i = parents[i];
#ifdef REPORT_PATH_LENGTHS
      pathlen++;
#endif
    }
#ifdef REPORT_PATH_LENGTHS
    max_pathlen.update_value(pathlen);
    total_pathlen.update_value(pathlen);
#endif
    return i;
  }

  inline uintE find_compress(uintE i, pbbs::sequence<parent>& parents) {
#ifdef REPORT_PATH_LENGTHS
    uintE pathlen = 1;
#endif
    parent j = i;
    if (parents[j] == j) return j;
    do {
      j = parents[j];
#ifdef REPORT_PATH_LENGTHS
      pathlen++;
#endif
    } while (parents[j] != j);
    parent tmp;
    while ((tmp=parents[i])>j) {
      parents[i] = j; i=tmp;
    }
#ifdef REPORT_PATH_LENGTHS
    max_pathlen.update_value(pathlen);
    total_pathlen.update_value(pathlen);
#endif
    return j;
  }

  inline uintE find_atomic_split(uintE i, pbbs::sequence<parent>& parents) {
#ifdef REPORT_PATH_LENGTHS
    uintE pathlen = 1;
#endif
    while(1) {
      parent v = parents[i];
      parent w = parents[v];
      if (v == w) {
#ifdef REPORT_PATH_LENGTHS
    max_pathlen.update_value(pathlen);
    total_pathlen.update_value(pathlen);
#endif
        return v;
      }
      else {
        pbbs::atomic_compare_and_swap(&parents[i],v,w);
        // i = its parents
        i = v;
      }
#ifdef REPORT_PATH_LENGTHS
      pathlen++;
#endif
    }
  }

  inline uintE find_atomic_halve(uintE i, pbbs::sequence<parent>& parents) {
#ifdef REPORT_PATH_LENGTHS
    uintE pathlen = 1;
#endif
    while(1) {
      parent v = parents[i];
      parent w = parents[v];
      if(v == w) {
#ifdef REPORT_PATH_LENGTHS
    max_pathlen.update_value(pathlen);
    total_pathlen.update_value(pathlen);
#endif
        return v;
      } else {
        pbbs::atomic_compare_and_swap(&parents[i],(parent)v,(parent)w);
        // i = its grandparent
        i = parents[i];
      }
#ifdef REPORT_PATH_LENGTHS
      pathlen++;
#endif
    }
  }

  inline uintE find_split(uintE i, pbbs::sequence<parent>& parents) {
#ifdef REPORT_PATH_LENGTHS
    uintE pathlen = 1;
#endif
    while(1) {
      parent v = parents[i];
      parent w = parents[v];
      if(v == w) {
#ifdef REPORT_PATH_LENGTHS
        max_pathlen.update_value(pathlen);
        total_pathlen.update_value(pathlen);
#endif
        return v;
      }
      else {
        parents[i] = w;
        i = v;
        std::atomic_thread_fence(std::memory_order_seq_cst);
      }
#ifdef REPORT_PATH_LENGTHS
      pathlen++;
#endif
    }
  }

  inline uintE find_halve(uintE i, pbbs::sequence<parent>& parents) {
#ifdef REPORT_PATH_LENGTHS
    uintE pathlen = 1;
#endif
    while(1) {
      parent v = parents[i];
      parent w = parents[v];
      if(v == w) {
#ifdef REPORT_PATH_LENGTHS
        max_pathlen.update_value(pathlen);
        total_pathlen.update_value(pathlen);
#endif
        return v;
      }
      else {
        parents[i] = w;
        //i = w;
        i = parents[i];
        std::atomic_thread_fence(std::memory_order_seq_cst);
      }
#ifdef REPORT_PATH_LENGTHS
      pathlen++;
#endif
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
  inline uintE splice(uintE u, uintE v, pbbs::sequence<parent>& parents) {
    parent z = parents[u];
    parents[u] = parents[v];
    return z;
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
      uintE tries = 0;
      while(1) {
        tries++;
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
#ifdef REPORT_MAX_TRIES
      max_uf_tries.update_value(tries);
      total_uf_tries.update_value(tries);
#endif
    }
  };

  struct UniteRemLock {
    uintE n;
    std::mutex* locks;
    UniteRemLock(uintE n) : n(n) {
      locks = pbbs::new_array<std::mutex>(n);
    }

    ~UniteRemLock() {
      pbbs::free_array(locks);
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& parents) {
      parent rx = u_orig;
      parent ry = v_orig;
      parent z;
#ifdef REPORT_PATH_LENGTHS
      uintE pathlen = 1;
#endif
      while (parents[rx] != parents[ry]) {
        /* link from high -> low */
        if (parents[rx] < parents[ry]) std::swap(rx,ry);
        if (rx == parents[rx]) {
          locks[rx].lock();
          parent py = parents[ry];
          /* link high -> low */
          if (rx == parents[rx] && rx > py) {
            parents[rx] = py;
          }
          locks[rx].unlock();
        } else {
          z = parents[rx];
          pbbs::atomic_compare_and_swap(&parents[rx], z, parents[ry]);
          rx = z;
        }
#ifdef REPORT_PATH_LENGTHS
        pathlen++;
#endif
      }
#ifdef REPORT_PATH_LENGTHS
      max_pathlen.update_value(pathlen);
      total_pathlen.update_value(pathlen);
#endif
      return;
    }
  };

  template <class Splice, class Compress, FindOption find_option>
  struct UniteRemCAS {
    Splice& splice;
    Compress& compress;
    UniteRemCAS(Splice& splice, Compress& compress) : splice(splice), compress(compress) { }

    inline void operator()(uintE x, uintE y, pbbs::sequence<parent>& parents) {
      uintE rx = x; uintE ry = y;
#ifdef REPORT_PATH_LENGTHS
      uintE pathlen = 1;
#endif
      while (parents[rx] != parents[ry]) {
        /* link high -> low */
        parent p_ry = parents[ry];
        if (parents[rx] < p_ry) {
          std::swap(rx, ry);
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
#ifdef REPORT_PATH_LENGTHS
        pathlen++;
#endif
      }
#ifdef REPORT_PATH_LENGTHS
      max_pathlen.update_value(pathlen);
      total_pathlen.update_value(pathlen);
#endif
      return;
    }
  };

  struct UniteEarly {
    UniteEarly() {}
    inline void operator()(uintE u, uintE v, pbbs::sequence<parent>& parents) {
      uintE tries = 0;
      while(u != v) {
        tries++;
        /* link high -> low */
        if(v > u) std::swap(u,v);
        if (parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) { return; }
        parent z = parents[u];
        parent w = parents[z];
        pbbs::atomic_compare_and_swap(&parents[u],z,w);
        u = w;
      }
#ifdef REPORT_MAX_TRIES
      max_uf_tries.update_value(tries);
      total_uf_tries.update_value(tries);
#endif
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
      uintE tries = 0;
      while(1) {
        tries++;
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
#ifdef REPORT_MAX_TRIES
      max_uf_tries.update_value(tries);
      total_uf_tries.update_value(tries);
#endif
    }
  };

  /* Add unite-by-size (lock-based?) */

} // namespace unite_variants
