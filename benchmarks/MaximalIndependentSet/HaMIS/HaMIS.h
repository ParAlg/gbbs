// HaMIS: asynchronous maximal independent set, in the spirit of HaMM.
//
// Implicitly orient each edge from the higher-id endpoint to the lower-id
// endpoint. A vertex v sequentially scans its lower-id neighbors (CSRs are
// sorted ascending, so we stop as soon as we hit a neighbor w >= v):
//   - if any lower-id neighbor is IN_MIS, v becomes COVERED;
//   - if a lower-id neighbor is still UNDECIDED, v helps resolve it;
//   - once all lower-id neighbors are COVERED, v CASes itself UNDECIDED -> IN_MIS
//     and then walks its full adjacency list, CASing each UNDECIDED neighbor
//     to COVERED.
//
// Recursion is replaced by an explicit per-thread work stack, so the maximum
// resolution depth (a long ascending-id chain) is bounded by heap, not the
// program stack. Each frame holds (v, idx-into-adj-list); when v's scan hits
// an UNDECIDED lower-id neighbor w, we save v's position and push (w, 0).
// When w finishes, we resume v at the saved idx and re-read state[w] (now
// IN_MIS or COVERED) without recursing.
//
// State encoding (one byte per vertex):
//   0 = UNDECIDED, 1 = COVERED, 2 = IN_MIS.
// Only a vertex writes itself to IN_MIS; only neighbors write a vertex to
// COVERED. The UNDECIDED -> IN_MIS CAS can fail only if a concurrent neighbor
// transitioned to IN_MIS first and marked us COVERED, so a failed CAS is
// equivalent to becoming COVERED.

#pragma once

#include "gbbs/gbbs.h"

#include <vector>

namespace gbbs {
namespace hamis {

constexpr uint8_t kUndecided = 0;
constexpr uint8_t kCovered = 1;
constexpr uint8_t kInMis = 2;

struct Frame {
  uintE v;
  uintE idx;  // next neighbor index to examine
};

template <template <class W> class vertex, class W>
inline void Resolve(uint8_t* state, symmetric_graph<vertex, W>& G,
                    uintE start_v, std::vector<Frame>& stack) {
  if (state[start_v] != kUndecided) return;
  stack.push_back({start_v, 0});

  while (!stack.empty()) {
    const size_t top = stack.size() - 1;
    const uintE v = stack[top].v;

    if (state[v] != kUndecided) {
      stack.pop_back();
      continue;
    }

    auto nghs = G.get_vertex(v).out_neighbors();
    const uintE d = G.get_vertex(v).out_degree();

    uintE idx = stack[top].idx;
    bool forced_covered = false;
    bool pushed = false;
    uintE push_target = 0;

    while (idx < d) {
      const uintE w = nghs.get_neighbor(idx);
      // CSR adjacency lists are sorted ascending; once w >= v, all
      // remaining neighbors are also >= v, so we stop.
      if (w >= v) break;
      const uint8_t s = state[w];
      if (s == kUndecided) {
        pushed = true;
        push_target = w;
        break;
      }
      if (s == kInMis) {
        forced_covered = true;
        break;
      }
      // s == kCovered: keep scanning.
      ++idx;
    }

    if (pushed) {
      // Save position; when we return to v, state[push_target] will be
      // resolved and we'll re-read it at this idx without recursing.
      stack[top].idx = idx;
      stack.push_back({push_target, 0});
      continue;
    }

    if (forced_covered) {
      gbbs::atomic_compare_and_swap(&state[v], kUndecided, kCovered);
    } else if (gbbs::atomic_compare_and_swap(&state[v], kUndecided, kInMis)) {
      // Mark every UNDECIDED neighbor COVERED. Lower-id neighbors are
      // already COVERED by invariant; CAS is a cheap no-op for them.
      for (uintE i = 0; i < d; ++i) {
        const uintE w = nghs.get_neighbor(i);
        if (state[w] == kUndecided) {
          gbbs::atomic_compare_and_swap(&state[w], kUndecided, kCovered);
        }
      }
    }
    // Otherwise: a concurrent helper finished v, or an IN_MIS neighbor
    // marked us COVERED. Either way, v is now resolved.
    stack.pop_back();
  }
}

template <template <class W> class vertex, class W>
inline sequence<uint8_t> MaximalIndependentSet(
    symmetric_graph<vertex, W>& G) {
  size_t n = G.n;
  auto state = sequence<uint8_t>(n, kUndecided);

  timer mis_t;
  mis_t.start();
  parlay::parallel_for(0, n, [&](uintE v) {
    thread_local std::vector<Frame> stack;
    stack.clear();
    Resolve(state.begin(), G, v, stack);
  });
  double mis_time = mis_t.stop();
  std::cout << "### MIS Time: " << mis_time << std::endl;

  return state;
}

template <template <class W> class vertex, class W, class Seq>
inline void verify_mis(symmetric_graph<vertex, W>& G, const Seq& state) {
  size_t n = G.n;
  auto ok = parlay::sequence<bool>::from_function(n, [](size_t) { return true; });

  parlay::parallel_for(0, n, [&](size_t i) {
    auto count_mis = [&](const uintE& src, const uintE& ngh, const W&) {
      return src != ngh && state[ngh] == kInMis;
    };
    size_t ct = G.get_vertex(i).out_neighbors().count(count_mis);
    if (state[i] == kInMis) {
      // No neighbor may also be in the MIS.
      ok[i] = (ct == 0);
    } else if (state[i] == kCovered) {
      // At least one neighbor must be in the MIS.
      ok[i] = (ct > 0);
    } else {
      // No vertex should be left UNDECIDED.
      ok[i] = false;
    }
  });

  size_t n_ok = parlay::reduce(parlay::delayed_seq<size_t>(
      n, [&](size_t i) -> size_t { return ok[i] ? 1 : 0; }));
  size_t mis_size = parlay::reduce(parlay::delayed_seq<size_t>(
      n, [&](size_t i) -> size_t { return state[i] == kInMis ? 1 : 0; }));
  if (n_ok == n) {
    std::cout << "MIS OK! independent set size is: " << mis_size << "\n";
  } else {
    std::cout << "MIS invalid---" << (n - n_ok)
              << " vertices saw bad neighborhoods.\n";
  }
}

}  // namespace hamis
}  // namespace gbbs
