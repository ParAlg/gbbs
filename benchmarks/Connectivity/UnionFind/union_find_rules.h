#pragma once

#include "ligra/bridge.h"
#include "pbbslib/seq.h"
#include "jayanti.h"

#include <mutex>

namespace find_variants {
  inline uintE find_naive(uintE i, pbbs::sequence<uintE>& parent) {
    while(i != parent[i])
      i = parent[i];
    return i;
  }

  inline uintE find_compress(uintE i, pbbs::sequence<uintE>& parent) {
    uintE j = i;
    if (parent[j] == j) return j;
    do {
      j = parent[j];
    } while (parent[j] != j);
    uintE tmp;
    while ((tmp=parent[i])<j) {
      parent[i]=j; i=tmp;
    }
    return j;
  }

  inline uintE find_atomic_split(uintE i, pbbs::sequence<uintE>& parent) {
    while(1) {
      uintE v = parent[i];
      uintE w = parent[v];
      if(v == w) return v;
      else {
        pbbs::atomic_compare_and_swap(&parent[i],v,w);
        // i = its parent
        i = v;
      }
    }
  }

  inline uintE find_atomic_halve(uintE i, pbbs::sequence<uintE>& parent) {
    while(1) {
      uintE v = parent[i];
      uintE w = parent[v];
      if(v == w) return v;
      else {
        pbbs::atomic_compare_and_swap(&parent[i],v,w);
        // i = its grandparent
        i = parent[i];
      }
    }
  }

  inline uintE find_split(uintE i, pbbs::sequence<uintE>& parent) {
    while(1) {
      uintE v = parent[i];
      uintE w = parent[v];
      if(v == w) return v;
      else {
        parent[i] = w;
        i = v;
      }
    }
  }

  inline uintE find_halve(uintE i, pbbs::sequence<uintE>& parent) {
    while(1) {
      uintE v = parent[i];
      uintE w = parent[v];
      if(v == w) return v;
      else {
        parent[i] = w;
        //i = w;
        i = parent[i];
      }
    }
  }



  /* Used in Rem-CAS variants */
  // compresses the path up to root
  inline void compress(uintE start, uintE root, pbbs::sequence<uintE>& parent) {
    uintE tmp;
    while ((tmp = parent[start]) < root) {
      parent[start] = root;
      start = tmp;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE split_atomic_one(uintE i, uintE x, pbbs::sequence<uintE>& parent) {
    uintE v = parent[i];
    uintE w = parent[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parent[i],v,w);
      i = v;
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE halve_atomic_one(uintE i, uintE x, pbbs::sequence<uintE>& parent) {
    uintE v = parent[i];
    uintE w = parent[v];
    if(v == w) return v;
    else {
      pbbs::atomic_compare_and_swap(&parent[i],v,w);
      //i = w;
      i = parent[i];
    }
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE splice(uintE u, uintE v, pbbs::sequence<uintE>& parent) {
    uintE z = parent[u];
    parent[u] = parent[v];
    return z;
  }

  /* Used in Rem-CAS variants for splice */
  inline uintE splice_atomic(uintE u, uintE v, pbbs::sequence<uintE>& parent) {
    uintE z = parent[u];
    pbbs::atomic_compare_and_swap(&parent[u], z, parent[v]);
    return z;
  }
} // namespace find_variants




namespace unite_variants {

  template <class Find>
  struct Unite {
    Find& find;
    Unite(Find& find) : find(find) {
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents) {
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
      std::cout << "created locks;" << std::endl;
    }

    ~UniteRemLock() {
      pbbs::free_array(locks);
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parent) {
      uintE rx = u_orig;
      uintE ry = v_orig;
      uintE z;
      while (parent[rx] != parent[ry]) {
        if (parent[rx] < parent[ry]) std::swap(rx,ry);
        if (rx == parent[rx]) {
          locks[rx].lock();
          if (rx == parent[rx]) {
    	parent[rx] = parent[ry];
          }
          locks[rx].unlock();
        } else {
          z = parent[rx];
          parent[rx] = parent[ry];
          rx = z;
        }
      }
      return;
    }
  };

  template <class Splice, class Compress>
  struct UniteRemCAS {
    Splice splice;
    Compress compress;
    UniteRemCAS() {
      splice = Splice();
      compress = Compress();
    }

    inline void operator()(uintE x, uintE y, pbbs::sequence<uintE>& parent) {
      uintE rx = x; uintE ry = y;
      while (parent[rx] != parent[ry]) {
        if (parent[rx] > parent[ry]) {
          std::swap(rx, ry);
        }
        if (rx == parent[rx] && pbbs::atomic_compare_and_swap(&parent[rx], rx, parent[ry])) {
          // success
          compress(x, parent);
          compress(y, parent);
          return;
        } else {
          // failure: locally compress by splicing and try again
          rx = splice(rx, ry, parent);
        }
      }
      return;
      exit(0);
    }
  };

  struct UniteEarly {
    UniteEarly() {}
    inline void operator()(uintE u, uintE v, pbbs::sequence<uintE>& parents) {
      while(1) {
        if(u == v) return;
        if(v > u) std::swap(u,v);
        if (parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) { return; }
        uintE z = parents[u];
        uintE w = parents[z];
        pbbs::atomic_compare_and_swap(&parents[u],z,w);
        u = z;
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

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents) {
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
