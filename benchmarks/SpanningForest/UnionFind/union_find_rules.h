#pragma once

#include "ligra/bridge.h"
#include "pbbslib/seq.h"
#include "jayanti.h"
#include "benchmarks/Connectivity/common.h"

#include <mutex>

namespace find_variants {
  inline uintE find_naive(uintE i, pbbs::sequence<parent>& Parents) {
    uintE pathlen = 1;
    while(i != Parents[i]) {
      i = Parents[i];
      pathlen++;
    }
    report_pathlen(pathlen);
    return i;
  }

  inline uintE find_compress(uintE i, pbbs::sequence<parent>& Parents) {
    uintE pathlen = 1;
    parent j = i;
    if (Parents[j] == j) return j;
    do {
      j = Parents[j];
      pathlen++;
    } while (Parents[j] != j);
    parent tmp;
    while ((tmp=Parents[i])>j) {
      Parents[i] = j; i=tmp;
    }
    report_pathlen(pathlen);
    return j;
  }

  inline uintE find_atomic_split(uintE i, pbbs::sequence<parent>& Parents) {
    uintE pathlen = 1;
    while(1) {
      parent v = Parents[i];
      parent w = Parents[v];
      if (v == w) {
        report_pathlen(pathlen);
        return v;
      }
      else {
        pbbs::atomic_compare_and_swap(&Parents[i],v,w);
        // i = its Parents
        i = v;
      }
      pathlen++;
    }
  }

  inline uintE find_atomic_halve(uintE i, pbbs::sequence<parent>& Parents) {
    uintE pathlen = 1;
    while(1) {
      parent v = Parents[i];
      parent w = Parents[v];
      if(v == w) {
        report_pathlen(pathlen);
        return v;
      } else {
        pbbs::atomic_compare_and_swap(&Parents[i],(parent)v,(parent)w);
        // i = its grandparent
        i = Parents[i];
      }
      pathlen++;
    }
  }
} // namespace find_variants


namespace splice_variants {

  /* Used in Rem-CAS variants for splice */
  inline uintE split_atomic_one(uintE i, uintE x, pbbs::sequence<parent>& Parents) {
    parent v = Parents[i];
    parent w = Parents[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&Parents[i],v,w);
      i = v;
      return i;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE halve_atomic_one(uintE i, uintE x, pbbs::sequence<parent>& Parents) {
    parent v = Parents[i];
    parent w = Parents[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&Parents[i],v,w);
      i = w;
      return i;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE splice_atomic(uintE u, uintE v, pbbs::sequence<parent>& Parents) {
    parent z = Parents[u];
    pbbs::atomic_compare_and_swap(&Parents[u], z, Parents[v]);
    return z;
  }
} // namespace splice_variants


namespace unite_variants {

  template <class Find>
  struct Unite {
    Find& find;
    Unite(Find& find) : find(find) { }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& Parents, pbbs::sequence<edge>& Edges) {
      parent u = u_orig;
      parent v = v_orig;
      uintE tries = 0;
      while(1) {
        tries++;
        u = find(u,Parents);
        v = find(v,Parents);
        if(u == v) break;
        else if (u > v && Parents[u] == u && pbbs::atomic_compare_and_swap(&Parents[u],u,v)) {
          Edges[u] = std::make_pair(u_orig, v_orig);
          break;
        }
        else if (v > u && Parents[v] == v && pbbs::atomic_compare_and_swap(&Parents[v],v,u)) {
          Edges[v] = std::make_pair(u_orig, v_orig);
          break;
        }
      }
      report_tries(tries);
    }
  };

  template <class Compress, class Splice>
  struct UniteRemLock {
    uintE n;
    std::mutex* locks;
    Compress& compress;
    Splice& splice;
    UniteRemLock(Compress& compress, Splice& splice, uintE n) : n(n), compress(compress), splice(splice) {
      locks = pbbs::new_array<std::mutex>(n);
    }

    ~UniteRemLock() {
      pbbs::free_array(locks);
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& Parents, pbbs::sequence<edge>& Edges) {
      parent rx = u_orig;
      parent ry = v_orig;
      parent z;
      uintE pathlen = 1;
      while (Parents[rx] != Parents[ry]) {
        /* link from high -> low */
        if (Parents[rx] < Parents[ry]) std::swap(rx,ry);
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
      return;
    }
  };

  template <class Splice, class Compress, FindOption find_option>
  struct UniteRemCAS {
    Compress& compress;
    Splice& splice;
    UniteRemCAS(Compress& compress, Splice& splice) : compress(compress), splice(splice) { }

    inline void operator()(uintE x, uintE y, pbbs::sequence<parent>& Parents, pbbs::sequence<edge>& Edges) {
      uintE rx = x; uintE ry = y;
      uintE pathlen = 1;
      while (Parents[rx] != Parents[ry]) {
        /* link high -> low */
        parent p_ry = Parents[ry];
        if (Parents[rx] < p_ry) {
          std::swap(rx, ry);
        }
        if (rx == Parents[rx] && pbbs::atomic_compare_and_swap(&Parents[rx], rx, p_ry)) {
          Edges[rx] = std::make_pair(u_orig, v_orig);
          // success
          if constexpr (find_option != find_naive) { /* aka find_none */
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
      return;
    }
  };

  struct UniteEarly {
    UniteEarly() {}
    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& Parents, pbbs::sequence<edge>& Edges) {
      uintE u = u_orig;
      uintE v = v_orig;
      uintE tries = 1;
      while(u != v) {
        tries++;
        /* link high -> low */
        if(v > u) std::swap(u,v);
        if (Parents[u] == u && pbbs::atomic_compare_and_swap(&Parents[u],u,v)) {
          Edges[u] = std::make_pair(u_orig, v_orig);
          return;
        }
        parent z = Parents[u];
        parent w = Parents[z];
        pbbs::atomic_compare_and_swap(&Parents[u],z,w);
        u = w;
      }
      report_tries(tries);
    }
  };

  template <class Find>
  struct UniteND {
    Find find;
    pbbs::sequence<uintE> hooks;
    UniteND(size_t n, Find& find) : find(find) {
      hooks = pbbs::sequence<uintE>(n, UINT_E_MAX);
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<parent>& Parents, pbbs::sequence<edge>& Edges) {
      parent u = u_orig;
      parent v = v_orig;
      uintE tries = 0;
      while(1) {
        tries++;
        u = find(u,Parents);
        v = find(v,Parents);
        if(u == v) break;
        /* link high -> low */
        if(u < v) std::swap(u,v);
        if (hooks[u] == UINT_E_MAX && pbbs::atomic_compare_and_swap(&hooks[u], UINT_E_MAX,v)) {
          Edges[u] = std::make_pair(u_orig, v_orig);
          Parents[u] = v;
          break;
        }
      }
      report_tries(tries);
    }
  };

} // namespace unite_variants
