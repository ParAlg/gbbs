#pragma once
// Based on RemST/CC.C

// Usage:
// numactl -i all ./RemST -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include <algorithm>
#include "ligra.h"
#include "utils/stats.h"

#include <iostream>
#include <limits.h>
#include <vector>
#include <mutex>
#include <atomic>

namespace remsf {
struct notMax { bool operator() (uintE i) {return i < UINT_E_MAX;}};

inline void unite(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parent, pbbs::sequence<uintE>& hooks, pbbs::sequence<std::mutex>& locks) {
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
} //namespace remsf

namespace remcc {
struct notMax { bool operator() (uintE i) {return i < UINT_E_MAX;}};

inline void unite(uintE u_orig, uintE v_orig, pbbs::sequence<uintE>& parent, pbbs::sequence<std::mutex>& locks) {
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
    } else {
      z = parent[rx];
      parent[rx] = parent[ry];
      rx = z;
    }
  }
  return;
}
} //namespace remcc


// returns a component labeling
template <class G>
pbbs::sequence<uintE> RemSF(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  auto hooks = pbbs::sequence<uintE>(n, [&] (size_t i) { return UINT_E_MAX; });
  auto locks = pbbs::sequence<std::mutex>(n);

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        remsf::unite(u, v, parents, hooks, locks);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  });

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = async::find_compress(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  // filter UINT_E_MAX from hooks if spanning forest edges are desired.
  return parents;
}

// returns a component labeling
template <class G>
pbbs::sequence<uintE> RemCC(G& GA) {
  using W = typename G::weight_type;
  size_t n = GA.n;
  size_t m = GA.m;

  auto parents = pbbs::sequence<uintE>(n, [&] (size_t i) { return i; });
  auto locks = pbbs::sequence<std::mutex>(n);

  parallel_for(0, n, [&] (size_t i) {
    auto map_f = [&] (uintE u, uintE v, const W& wgh) {
      if (u < v) {
        remcc::unite(u, v, parents, locks);
      }
    };
    GA.get_vertex(i).mapOutNgh(i, map_f); // in parallel
  });

  timer ft; ft.start();
  parallel_for(0, n, [&] (size_t i) {
    parents[i] = async::find_compress(i,parents);
  });
  ft.stop(); debug(ft.reportTotal("find time"););
  return parents;
}
