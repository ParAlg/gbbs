#pragma once

#include "ligra/bridge.h"
#include "pbbslib/seq.h"
#include "jayanti.h"
#include "benchmarks/Connectivity/common.h"

#include <mutex>

namespace find_variants {
  inline uintE find_naive(uintE i, pbbs::sequence<parent>& parents) {
    while(i != parents[i])
      i = parents[i];
    return i;
  }

  inline uintE find_compress(uintE i, pbbs::sequence<parent>& parents) {
    uintE j = i;
    if (parents[j] == j) return j;
    do {
      j = parents[j];
    } while (parents[j] != j);
    uintE tmp;
    while ((tmp=parents[i])<j) {
      parents[i] = j; i=tmp;
    }
    return j;
  }

  inline uintE find_atomic_split(uintE i, pbbs::sequence<parent>& parents) {
    while(1) {
      uintE v = parents[i];
      uintE w = parents[v];
      if(v == w) return v;
      else {
        pbbs::atomic_compare_and_swap(&parents[i],v,w);
        // i = its parents
        i = v;
      }
    }
  }

  inline uintE find_atomic_halve(uintE i, pbbs::sequence<parent>& parents) {
    while(1) {
      uintE v = parents[i];
      uintE w = parents[v];
      if(v == w) return v;
      else {
        pbbs::atomic_compare_and_swap(&parents[i],(parent)v,(parent)w);
        // i = its grandparent
        i = parents[i];
      }
    }
  }

  inline uintE find_split(uintE i, pbbs::sequence<parent>& parents) {
    while(1) {
      uintE v = parents[i];
      uintE w = parents[v];
      if(v == w) return v;
      else {
        parents[i] = w;
        i = v;
        std::atomic_thread_fence(std::memory_order_seq_cst);
      }
    }
  }

  inline uintE find_halve(uintE i, pbbs::sequence<parent>& parents) {
    while(1) {
      uintE v = parents[i];
      uintE w = parents[v];
      if(v == w) return v;
      else {
        parents[i] = w;
        //i = w;
        i = parents[i];
        std::atomic_thread_fence(std::memory_order_seq_cst);
      }
    }
  }
} // namespace find_variants


namespace splice_variants {


  /* Used in Rem-CAS variants for splice */
  inline uintE split_atomic_one(uintE i, uintE x, pbbs::sequence<parent>& parents) {
    uintE v = parents[i];
    uintE w = parents[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parents[i],v,w);
      i = v;
      return i;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE halve_atomic_one(uintE i, uintE x, pbbs::sequence<parent>& parents) {
    uintE v = parents[i];
    uintE w = parents[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parents[i],v,w);
      //i = w;
      i = parents[i];
      return i;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE splice(uintE u, uintE v, pbbs::sequence<parent>& parents) {
    uintE z = parents[u];
    parents[u] = parents[v];
    return z;
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE splice_atomic(uintE u, uintE v, pbbs::sequence<parent>& parents) {
    uintE z = parents[u];
    pbbs::atomic_compare_and_swap(&parents[u], z, parents[v]);
    return z;
  }
} // namespace splice_variants


namespace unite_variants {

  template <class Find>
  struct Unite {
    Find& find;
    Unite(Find& find) : find(find) {
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& parents) {
      uintE u = u_orig;
      uintE v = v_orig;
      while(1) {
        u = find(u,parents);
        v = find(v,parents);
        if(u == v) break;
        else if (u < v && parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) {
          return;
        }
        else if (v < u && parents[v] == v && pbbs::atomic_compare_and_swap(&parents[v],v,u)) {
          return;
        }
      }
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
      uintE rx = u_orig;
      uintE ry = v_orig;
      uintE z;
      while (parents[rx] != parents[ry]) {
        if (parents[rx] > parents[ry]) std::swap(rx,ry);
        if (rx == parents[rx]) {
          locks[rx].lock();
          uintE py = parents[ry];
          if (rx == parents[rx] && rx < py) {
            parents[rx] = py;
          }
          locks[rx].unlock();
        } else {
          z = parents[rx];
          // parents[rx] = parents[ry];
          pbbs::atomic_compare_and_swap(&parents[rx], z, parents[ry]);
          rx = z;
        }
      }
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
      while (parents[rx] != parents[ry]) {
        if (parents[rx] > parents[ry]) {
          std::swap(rx, ry);
        }
        if (rx == parents[rx] && pbbs::atomic_compare_and_swap(&parents[rx], rx, parents[ry])) {
          // success
          if constexpr (find_option != find_naive) { /* aka find_none */
            compress(x, parents);
            compress(y, parents);
          }
          return;
        } else {
          // failure: locally compress by splicing and try again
          rx = splice(rx, ry, parents);
        }
      }
      return;
    }
  };

  struct UniteEarly {
    UniteEarly() {}
    inline void operator()(uintE u, uintE v, pbbs::sequence<parent>& parents) {
      while(1) {
        if(u == v) return;
        if(v > u) std::swap(u,v);
        if (parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) { return; }
        uintE z = parents[u];
        uintE w = parents[z];
        pbbs::atomic_compare_and_swap(&parents[u],z,w);
        u = w;
      }
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
      uintE u = u_orig;
      uintE v = v_orig;
      while(1) {
        u = find(u,parents);
        v = find(v,parents);
        if(u == v) break;
        if(u > v) std::swap(u,v);
        if (hooks[u] == UINT_E_MAX && pbbs::atomic_compare_and_swap(&hooks[u], UINT_E_MAX,v)) {
          parents[u] = v;
          break;
        }
      }
    }
  };

  /* Add unite-by-size (lock-based?) */

} // namespace unite_variants
