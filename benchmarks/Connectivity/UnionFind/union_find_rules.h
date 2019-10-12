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

  template <class Splice, class Compress>
  struct UniteRemCAS {
    Splice& splice;
    Compress& compress;
    UniteRemCAS(uintE n, Splice& splice, Compress& compress) : splice(splice), compress(compress) {}

    inline void unite(uintE rx, uintE ry, pbbs::sequence<uintE>& parent, pbbs::sequence<uintE>& hooks) {
      assert(false);
      exit(0);
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

//  inline void unite(uintT i, edge<uintT>* E, uintT* parent, uintT* hooks) {
//  uintT rx = E[i].u;
//  uintT ry = E[i].v;
//  uintT z;
//  while (parent[rx] != parent[ry]) {
//    if (parent[rx] > parent[ry]) {
//      swap(rx, ry);
//    }
//
//    if (rx == parent[rx]) {
//      if (__sync_bool_compare_and_swap(&parent[rx], rx, parent[ry])) {
//        hooks[rx] = i;
//#ifdef COMPRESS
//        compress(E[i].u, parent[ry], parent);
//        compress(E[i].v, parent[ry], parent);
//#endif
//
//#ifdef SPLIT
//        split(E[i].u, parent);
//        split(E[i].v, parent);
//#endif
//
//#ifdef HALVE
//        halve(E[i].u, parent);
//        halve(E[i].v, parent);
//#endif
//        return;
//      }
//    } else {
//#ifdef SPLICE
//      rx = splice(rx, ry, parent);
//#endif
//
//#ifdef SPLICE_CAS
//      rx = splice_CAS(rx, ry, parent);
//#endif
//
//#ifdef SPLIT_ONE
//      rx = split_one(rx, parent);
//#endif
//
//    }
//  }
//  return;
//  }


  template <class Find>
  struct Unite {
    Find& find;
    Unite(uintE n, Find& find) : find(find) {}

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents, pbbs::sequence<uintE>& hooks) {
      uintE u = u_orig;
      uintE v = v_orig;
      while(1) {
        u = find(u,parents);
        v = find(v,parents);
        if(u == v) break;
        else if (u < v && parents[u] == u && pbbs::atomic_compare_and_swap(&parents[u],u,v)) {
          hooks[u] = v;
          return;
        }
        else if (v < u && parents[v] == v && pbbs::atomic_compare_and_swap(&parents[v],v,u)) {
          hooks[v] = u;
          return;
        }
      }
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

  template <class Find>
  struct UniteRemLock {
    uintE n;
    std::mutex* locks;
    UniteRemLock(uintE n, Find& find) : n(n) {
      locks = pbbs::new_array<std::mutex>(n);
      std::cout << "created locks;" << std::endl;
    }

    ~UniteRemLock() {
      pbbs::free_array(locks);
    }

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parent, pbbs::sequence<uintE>& hooks) {
      uintE rx = u_orig;
      uintE ry = v_orig;
      uintE z;
      while (parent[rx] != parent[ry]) {
        if (parent[rx] < parent[ry]) std::swap(rx,ry);
        if (rx == parent[rx]) {
          locks[rx].lock();
          bool success = false;
          if (rx == parent[rx]) {
    	parent[rx] = parent[ry];
    	success = true;
          }
          locks[rx].unlock();
          if (success) {
    	hooks[rx] = u_orig; // todo
    	return;
          }
        } else {
          z = parent[rx];
          parent[rx] = parent[ry];
          rx = z;
        }
      }
      return;
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


  template <class Find>
  struct UniteEarly {
    UniteEarly(uintE n, Find& find) {}
    inline void operator()(uintE u, uintE v, pbbs::sequence<uintE>& parents, pbbs::sequence<uintE>& hooks) {
      while(1) {
        if(u == v) return;
        if(v > u) std::swap(u,v);
        if(pbbs::atomic_compare_and_swap(&parents[u],u,v)) { hooks[u] = v; return; }
        uintE z = parents[u];
        uintE w = parents[z];
        pbbs::atomic_compare_and_swap(&parents[u],z,w);
        u = z;
      }
    }
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
    Find& find;
    UniteND(uintE n, Find& find) : find(find) {}

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents, pbbs::sequence<uintE>& hooks) {
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

    inline void operator()(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parents) {
      assert(false);
      exit(0);
    }
  };

  /* Add unite-by-size (lock-based?) */

} // namespace unite_variants

namespace union_find {

  template <class Find, class Unite, class Graph>
  struct UFAlgorithm {
    Graph& GA;
    Unite& unite;
    Find& find;
    bool use_hooks;
    UFAlgorithm(Graph& GA, Unite& unite, Find& find, bool use_hooks=false) :
      GA(GA), unite(unite), find(find), use_hooks(use_hooks) {}

    void components(pbbs::sequence& parents, long frequent_comp=-1) {
      size_t n = GA.n;

      timer ut; ut.start();
      parallel_for(0, n, [&] (size_t u) {
	// Only process edges for vertices not linked to the main component
	// note that this is safe only for undirected graphs. For directed graphs,
	// the in-edges must also be explored for all vertices.
	if (parents[u] != frequent_comp) {
	  auto map_f = [&] (uintE u, uintE v, const W& wgh) {
	    if (u < v) {
	      if (use_hooks) {
		unite(u, v, parents, hooks);
	      } else {
		unite(u, v, parents);
	      }
	    }
	  };
	  GA.get_vertex(u).mapOutNgh(u, map_f); // in parallel
	}
      }, 1);
      ut.stop(); ut.reportTotal("union time");

      timer ft; ft.start();
      parallel_for(0, n, [&] (size_t i) {
	parents[i] = find(i,parents);
      });
      ft.stop(); ft.reportTotal("find time");

    }


  }

} // namespace union_find
