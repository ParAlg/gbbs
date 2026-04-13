// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/helpers/speculative_for.h"

namespace gbbs {

class PerWorkerStats {
 private:
  // Define the cache line size as requested.
  // 128 bytes is a common L2 cache line size and a safe bet.
  static constexpr size_t CACHE_LINE_SIZE = 128;

  struct alignas(CACHE_LINE_SIZE) PaddedCounter {
    size_t count = 0;
  };

  size_t num_workers;
  std::unique_ptr<PaddedCounter[]> per_worker_storage;

 public:
  PerWorkerStats() : num_workers(parlay::num_workers()) {
    per_worker_storage = std::make_unique<PaddedCounter[]>(num_workers);
  }

  void increment() {
    size_t id = parlay::worker_id();
    per_worker_storage[id].count++;
  }

  size_t finalize() {
    size_t total_count = 0;
    for (size_t i = 0; i < num_workers; ++i) {
      total_count += per_worker_storage[i].count;
    }
    return total_count;
  }
};

template <typename T>
class PerWorkerBag {
private:
    // Define the cache line size as requested.
    // 128 bytes is a common L2 cache line size and a safe bet.
    static constexpr size_t CACHE_LINE_SIZE = 128;

    /**
     * @brief A wrapper struct that aligns its contents to a cache line.
     *
     * `alignas(CACHE_LINE_SIZE)` ensures that any object of this type
     * starts on a 128-byte boundary. It also forces the total
     * `sizeof(PaddedStorage)` to be a multiple of 128.
     *
     * In an array of `PaddedStorage`, this guarantees that
     * `array[i]` and `array[i+1]` are in different cache lines.
     */
    struct alignas(CACHE_LINE_SIZE) PaddedStorage {
        std::vector<T> vec;
    };

    // We need to store the number of workers, since the unique_ptr
    // to an array doesn't store its own size.
    size_t num_workers;

    // A smart pointer to a C-style array of PaddedStorage objects.
    // `std::make_unique<T[]>(N)` correctly allocates an array of N
    // objects, respecting the `alignas` specifier.
    std::unique_ptr<PaddedStorage[]> per_worker_storage;

public:
    /**
     * @brief Constructs the PerWorkerBag.
     *
     * Allocates an array of padded, aligned vectors, one for each worker.
     */
    PerWorkerBag() : num_workers(parlay::num_workers()) {
      // `std::make_unique` for arrays (C++14) will correctly
      // allocate and initialize an array of aligned PaddedStorage objects.
      per_worker_storage = std::make_unique<PaddedStorage[]>(num_workers);
    }

    // --- Thread-Safe Methods (callable in parallel) ---

    /**
     * @brief Adds an element to the bag by copying it.
     *
     * Accesses the correct worker's private (and padded) vector.
     *
     * @param element The element to add (lvalue).
     */
    void add(const T& element) {
        size_t id = parlay::worker_id();
        per_worker_storage[id].vec.push_back(element);
    }

    /**
     * @brief Adds an element to the bag by moving it.
     *
     * Accesses the correct worker's private (and padded) vector.
     *
     * @param element The element to add (rvalue).
     */
    void add(T&& element) {
        size_t id = parlay::worker_id();
        per_worker_storage[id].vec.push_back(std::move(element));
    }


    // --- Serial-Only Method ---

    /**
     * @brief Combines all elements from all worker bags into a single vector.
     *
     * @warning This method is **NOT** thread-safe with respect to `add`.
     * It should only be called *after* all parallel threads have
     * finished adding elements.
     *
     * @return A std::vector<T> containing all elements.
     */
    parlay::sequence<T> finalize() {
      // 1. Calculate the total size needed.
      size_t total_size = 0;
      for (size_t i = 0; i < num_workers; ++i) {
          total_size += per_worker_storage[i].vec.size();
      }

      // 2. Reserve space.
      parlay::sequence<T> result;
      result.reserve(total_size);

      // 3. Move elements from each worker bag into the result vector.
      for (size_t i = 0; i < num_workers; ++i) {
          auto& bag = per_worker_storage[i].vec;
          std::move(bag.begin(), bag.end(), std::back_inserter(result));
          bag.clear();
      }
      return result;
    }

    /**
     * @brief Clears all per-worker bags for reuse.
     *
     * Like finalize(), this should only be called serially.
     */
    void clear() {
        for (size_t i = 0; i < num_workers; ++i) {
            per_worker_storage[i].vec.clear();
        }
    }

    // --- Rule of 5/0 ---
    // Disable copying, as copying a parallel-use bag is undefined.
    PerWorkerBag(const PerWorkerBag&) = delete;
    PerWorkerBag& operator=(const PerWorkerBag&) = delete;
    // Default move operations are fine, as std::unique_ptr is movable.
    PerWorkerBag(PerWorkerBag&&) = default;
    PerWorkerBag& operator=(PerWorkerBag&&) = default;
    // Default destructor is fine, std::unique_ptr will handle cleanup.
    ~PerWorkerBag() = default;
};


template <class W, class Bag>
uintE HelpMatch(Bag& bag, PerWorkerStats& stats, uintE* P, uintE v, uintE c, parlay::sequence<std::atomic<bool>>& matched) {
  while (true) {
    auto p_v = P[v];
    if (p_v > v) return p_v;
    if (P[p_v] == p_v && gbbs::atomic_compare_and_swap(&P[p_v], p_v, v)) {
      matched[v] = true;
      matched[p_v] = true;
      bag.add(std::make_tuple(v, p_v, W{}));
      return p_v;
    }
    stats.increment();
    auto r = HelpMatch<W, Bag>(bag, stats, P, p_v, v, matched);
    if (r == v) return p_v;
    if (P[v] == p_v && gbbs::atomic_compare_and_swap(&P[v], p_v, c)) {
      matched[v] = true;
      matched[c] = true;
      bag.add(std::make_tuple(v, c, W{}));
      return c;
    }
    stats.increment();
  }
};


template <template <class W> class vertex, class W>
inline sequence<std::tuple<uintE, uintE, W>> MaximalMatching(symmetric_graph<vertex, W>& G) {
  using edge = std::tuple<uintE, uintE, W>;
  parlay::sequence<std::tuple<uintE, uintE, W>> output;

  auto P = sequence<uintE>::from_function(G.n, [&](uintE i) { return i; });

  PerWorkerBag<edge> bag;
  PerWorkerStats stats;
  PerWorkerStats help_matches;
  PerWorkerStats attempts;

  parlay::sequence<std::atomic<bool>> matched(G.n);

  parlay::parallel_for(0, G.n, [&] (uintE v) {
    uintE prev = v;
    auto map_f = [&](const uintE& v, const uintE& x, const W& wgh) {
      if (v > x) { // lower neighbor; try to match
        // if we fail the next CAS, we must be matched (and the edge
        // is already added).
        if (P[v] != prev) {
          return false;
        }

        // Early terminate.
        if (matched[x]) { return true; }

        if (!gbbs::atomic_compare_and_swap(&P[v], prev, x)) {
          stats.increment();
          return false;
        }
        attempts.increment();

        // If we succeed for the next CAS, we matched x to us. Add the
        // edge to the bag.
        if (P[x] == x && gbbs::atomic_compare_and_swap(&P[x], x, v)) {
          bag.add(std::make_tuple(x, v, W{}));
          matched[x] = true;
          matched[v] = true;
          return false;  // don't continue looping
        }
        stats.increment();



        // Otherwise, x is trying to match with someone else. Go help
        // it.
        help_matches.increment();
        if (HelpMatch<W>(bag, stats, P.begin(), x, v, matched) == v) {
          return false;
        }
        prev = x;
      }
      // Continue looping over the edges.
      return true;
    };
    G.get_vertex(v).out_neighbors().map_seq(map_f);
  });

  std::cout << "Failed CASes: " << stats.finalize() << "\n";
  std::cout << "Num help matches: " << help_matches.finalize() << "\n";
  std::cout << "Num vertex attempts: " << attempts.finalize() << std::endl;

  return bag.finalize();
}

template <template <class W> class vertex, class W, class Seq>
inline void verify_matching(symmetric_graph<vertex, W>& G, Seq& matching) {
  size_t n = G.n;
  auto ok = parlay::sequence<bool>::from_function(n, [](size_t i) { return 1; });
  auto matched = parlay::sequence<uintE>::from_function(n, [](size_t i) { return 0; });

  // Check that this is a valid matching
  parlay::parallel_for(0, matching.size(), [&](size_t i) {
    const auto& edge = matching[i];
    gbbs::write_add(&matched[std::get<0>(edge)], 1);
    gbbs::write_add(&matched[std::get<1>(edge)], 1);
  });

  bool valid = true;
  parlay::parallel_for(0, n, [&](size_t i) {
    if (matched[i] > 1) valid = false;
  });
  assert(valid == true);

  // Check maximality of the matching (skip self-loops: can't match a vertex
  // with itself)
  auto map2_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    if (src != ngh && !matched[src] && !matched[ngh]) {
      // could have added this edge, increasing the size of the matching
      ok[src] = 0;
      ok[ngh] = 0;
    }
  };
  parlay::parallel_for(0, n,
               [&](size_t i) { G.get_vertex(i).out_neighbors().map(map2_f); }, 1);

  auto ok_f = [&](size_t i) { return ok[i]; };
  auto ok_im = parlay::delayed_seq<size_t>(n, ok_f);
  size_t n_ok = parlay::reduce(ok_im);
  if (n == n_ok) {
    std::cout << "Matching OK! matching size is: " << matching.size() << "\n";
  } else {
    std::cout << "Matching invalid---" << (n - n_ok)
              << " vertices saw bad neighborhoods."
              << "\n";
  }
}

}  // namespace gbbs
