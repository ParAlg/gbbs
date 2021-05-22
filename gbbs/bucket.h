// This code is part of the project "Julienne: A Framework for Parallel Graph
// Algorithms using Work-efficient Bucketing", presented at Symposium on
// Parallelism in Algorithms and Architectures, 2017.
// Copyright (c) 2017 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Implementation of a bucketing structure as described in "Julienne: A
// Framework for Parallel Graph Algorithms using Work-efficient Bucketing"
// (SPAA'17)
//
// Note that in practice very few buckets are opened (while vertices may
// start off in a massive number of buckets). To optimize for this situation we
// only explicitly represent "total_buckets" buckets, and filter the remaining
// vertices that have not yet been processed every time the current range runs
// out of buckets. Experimenting with different values of total_buckets (the
// -nb parameter in applications) is necessary to obtain the best performance.
// This also means that the current code could be optimized to run much faster
// in a case where many buckets will be processed; please contact us if you have
// such a use-case.
#pragma once

#include <limits>
#include <tuple>
#include <cassert>
#include <optional>

#include "vertex_subset.h"
#include "bridge.h"


#define CACHE_LINE_S 64

namespace gbbs {

enum bucket_order { decreasing, increasing };

// Maintains a dynamic mapping from a set of ident_t's to a set of buckets with
// integer type bucket_t.
template <class D, class ident_t, class bucket_t>
struct buckets {
 public:
   using bucket_id = bucket_t; // for readability
  // node that bucket_dests are just bucket_t's (bucket_ids)

  struct bucket {
    size_t id;
    size_t num_filtered;
    sequence<ident_t> identifiers;
    bucket(size_t _id, sequence<ident_t>&& _identifiers)
        : id(_id), identifiers(std::move(_identifiers)) {
    }
  };

  using id_dyn_arr = sequence<ident_t>;

  const bucket_id null_bkt = std::numeric_limits<bucket_id>::max();

  // Create a bucketing structure.
  //   n : the number of identifiers
  //   d : map from identifier -> bucket
  //   order : the order to iterate over the buckets
  //   total_buckets: the total buckets to materialize
  //
  //   For an identifier i:
  //   d[i] is the bucket currently containing i
  //   d[i] = std::numeric_limits<bucket_id>::max() if i is not in any bucket
  buckets(size_t _n, D& _d, bucket_order _order, size_t _total_buckets)
      : n(_n),
        d(_d),
        order(_order),
        open_buckets(_total_buckets - 1),
        total_buckets(_total_buckets),
        cur_bkt(0),
        max_bkt(_total_buckets),
        num_elms(0),
        allocated(true) {
    // Initialize array consisting of the materialized buckets.
    bkts = sequence<id_dyn_arr>(total_buckets);
    prev_bkt_sizes = sequence<size_t>::uninitialized(total_buckets);
    bkt_pointers = sequence<ident_t*>::uninitialized(total_buckets);

    // Set the current range being processed based on the order.
    if (order == increasing) {
      auto imap_f = [&](size_t i) { return d[i]; };
      auto imap = parlay::delayed_seq<bucket_id>(n, imap_f);
      size_t min_b = parlay::reduce(imap, parlay::minm<bucket_id>());
      cur_range = min_b / open_buckets;
    } else if (order == decreasing) {
      auto imap_f = [&](size_t i) { return (d[i] == null_bkt) ? 0 : d[i]; };
      auto imap = parlay::delayed_seq<bucket_id>(n, imap_f);
      size_t max_b = parlay::reduce(imap, parlay::maxm<bucket_id>());
      cur_range = (max_b + open_buckets) / open_buckets;
    } else {
      std::cout << "Unknown order: " << order
                << ". Must be one of {increasing, decreasing}"
                << "\n";
      abort();
    }

    // Update buckets with all (id, bucket) pairs. Identifiers with bkt =
    // null_bkt are ignored by update_buckets.
    auto get_id_and_bkt = [&](ident_t i) -> std::optional<std::tuple<ident_t, bucket_id> > {
      bucket_id bkt = _d[i];
      if (bkt != null_bkt) {
        bkt = to_range(bkt);
      }
      return std::optional<std::tuple<ident_t, bucket_id> >(std::make_tuple(i, bkt));
    };
    update_buckets(get_id_and_bkt, n);
  }

  void del() {
    if (allocated) {
      for (size_t i = 0; i < total_buckets; i++) {
        bkts[i].clear();
      }
      bkts.clear();
      allocated = false;
    }
  }

  void report() {
    t1.reportTotal("unpack time");
    t4.reportTotal("nextbucket time");
    t2.reportTotal("parallel bucket time");
    t3.reportTotal("sequential bucket time");
  }

  // Returns the next non-empty bucket from the bucket structure. The return
  // value's bkt_id is null_bkt when no further buckets remain.
  inline bucket next_bucket() {
    while (!curBucketNonEmpty() && num_elms > 0) {
      _next_bucket();
    }
    if (num_elms == 0) {
      size_t bkt_num = null_bkt;  // no buckets remain
      return bucket(bkt_num, sequence<ident_t>());
    }
    return get_cur_bucket();
  }

  // Computes a bucket_dest for an identifier moving from bucket_id prev to
  // bucket_id next.
  inline bucket_id get_bucket(const bucket_id& prev,
                              const bucket_id& next) const {
    bucket_id pb = to_range(prev);
    bucket_id nb = to_range(next);
    if ((nb != null_bkt) &&
        ((prev == null_bkt) || (pb != nb) || (nb == cur_bkt))) {
      return nb;
    }
    return null_bkt;
  }

  // Computes a bucket_dest for an identifier moving to bucket_id next.
  inline bucket_id get_bucket(const bucket_id& next) const {
    bucket_t nb = to_range(next);
    // Note that the interface currently only implements strictly_decreasing
    // priority, which is why the code below does not check pri_order.
    if (order == increasing) {
      // case for strictly_decreasing priorities, assuming elements start out
      // in the structure.
      if (nb != null_bkt && nb != open_buckets) {
        return nb;
      } // case for strictly_increasing elided
    } else { // bkt_order == decreasing
      if (nb != null_bkt) {
      // strictly_decreasing priorities, assuming elements start out in the structure.
        return nb;
      }
    }
    return null_bkt;
  }

  template <class F>
  inline size_t update_buckets_seq(F& f, size_t k) {
    t3.start();
    size_t ne_before = num_elms;
    for (size_t i = 0; i < k; i++) {
      auto m = f(i);
      bucket_id bkt = std::get<1>(*m);
      if (m.has_value() && bkt != null_bkt) {
        bkts[bkt].emplace_back(std::get<0>(*m));
        num_elms++;
      }
    }
    t3.stop();
    return num_elms - ne_before;
  }

  // Updates k identifiers in the bucket structure. The i'th identifier and
  // its bucket_dest are given by F(i).
  template <class F>
  inline size_t update_buckets(F f, size_t k) {
    size_t num_blocks = k / 4096;
    int num_threads = num_workers();
    if (k < 4096 || num_threads == 1) {
      return update_buckets_seq(f, k);
    }
    t2.start();

    size_t ne_before = num_elms;

    size_t block_bits = parlay::log2_up(num_blocks);
    num_blocks = 1 << block_bits;
    size_t block_size = (k + num_blocks - 1) / num_blocks;

    auto hists_seq = sequence<bucket_id>::uninitialized((num_blocks + 1) * total_buckets * CACHE_LINE_S);
    auto hists = hists_seq.data();

    // 1. Compute per-block histograms
    parallel_for(0, num_blocks, [&] (size_t block_id) {
      size_t s = block_id * block_size;
      size_t e = std::min(s + block_size, k);
      bucket_id* hist = &(hists[block_id * total_buckets]);

      for (size_t j = 0; j < total_buckets; j++) {
        hist[j] = 0;
      }
      for (size_t j = s; j < e; j++) {
        auto m = f(j);
        bucket_id b = std::get<1>(*m);
        if (m.has_value() && b != null_bkt) {
          hist[b]++;
        }
      }
    }, 1);

    // 2. Aggregate histograms into a single histogram.
    auto get = [&](size_t i) {
      size_t col = i % num_blocks;
      size_t row = i / num_blocks;
      return hists[col * total_buckets + row];
    };

    size_t last_ind = (num_blocks * total_buckets);
    auto outs_seq = sequence<bucket_id>::uninitialized(last_ind + 1);
    auto outs = outs_seq.data();
    parallel_for(0, last_ind, [&] (size_t i) {
      outs[i] = get(i);
    });
    outs[last_ind] = 0;

    parlay::scan_inplace(parlay::make_slice(outs_seq), parlay::addm<bucket_id>());

    size_t* const prev_sizes = prev_bkt_sizes.data();
    auto bkt_ptrs = bkt_pointers.data();

    // 3. Resize buckets based on the summed histogram.
    for (size_t i = 0; i < total_buckets; i++) {
      size_t num_inc = outs[(i + 1) * num_blocks] - outs[i * num_blocks];
      size_t cur_size = bkts[i].size();
      prev_sizes[i] = cur_size;
      bkts[i].resize(cur_size + num_inc);
      bkt_ptrs[i] = bkts[i].data();
      num_elms += num_inc;
    }

    // 4. Compute the starting offsets for each block.
    parallel_for(0, total_buckets, [&] (size_t bucket_id) {
      size_t start = outs[bucket_id * num_blocks];
      for (size_t block_id = 0; block_id < num_blocks; block_id++) {
        hists[(bucket_id * num_blocks + block_id) * CACHE_LINE_S] =
            outs[bucket_id * num_blocks + block_id] - start;
      }
    }, 1);


    // 5. Iterate over blocks again. Insert (id, bkt) into bkt[hists[bkt]]
    // and increment hists[bkt].
    parallel_for(0, num_blocks, [&] (size_t block_id) {
      size_t s = block_id * block_size;
      size_t e = std::min(s + block_size, k);
      // our buckets are now spread out, across outs
      for (size_t j = s; j < e; j++) {
        auto m = f(j);
        bucket_id bucket_id = std::get<1>(*m);
        if (m.has_value() && bucket_id != null_bkt) {
          ident_t v = std::get<0>(*m);
          size_t ind = prev_sizes[bucket_id] + hists[(bucket_id * num_blocks + block_id) * CACHE_LINE_S];
          bkt_ptrs[bucket_id][ind] = v;
          hists[(bucket_id * num_blocks + block_id) * CACHE_LINE_S]++;
        }
      }
    }, 1);

//    // 6. Finally, update the size of each bucket.
//    for (size_t i = 0; i < total_buckets; i++) {
//      size_t num_inc = outs[(i + 1) * num_blocks] - outs[i * num_blocks];
//      size_t& m = bkts[i].size();
//      m += num_inc;
//    }

    t2.stop();
    return num_elms - ne_before;
  }


 private:
  size_t n;  // total number of identifiers in the system
  D& d;
  const bucket_order order;
  size_t open_buckets;
  size_t total_buckets;
  size_t cur_bkt;
  size_t max_bkt;
  size_t num_elms;
  bool allocated;

  size_t cur_range;
  sequence<id_dyn_arr> bkts;
  sequence<size_t> prev_bkt_sizes;
  sequence<ident_t*> bkt_pointers;

  timer t1;
  timer t2;
  timer t3;
  timer t4;

  inline bool curBucketNonEmpty() { return bkts[cur_bkt].size() > 0; }

  inline void unpack() {
    t1.start();
    size_t m = bkts[open_buckets].size();
    auto tmp_seq = std::move(bkts[open_buckets]);
    auto tmp = tmp_seq.data();
    bkts[open_buckets] = id_dyn_arr();

    if (order == increasing) {
      cur_range++;  // increment range
    } else {
      cur_range--;
    }

    auto g = [&](ident_t i) -> std::optional<std::tuple<ident_t, bucket_id> > {
      ident_t v = tmp[i];
      bucket_id bkt = to_range(d[v]);
      return std::optional<std::tuple<ident_t, bucket_id> >(std::make_tuple(v, bkt));
    };

    if (m != num_elms) {
      std::cout << "m = " << m << " num_elms = " << num_elms << "\n";
      cur_bkt = 0;
      std::cout << "curBkt = " << get_cur_bucket_num() << "\n";
      std::cout << "mismatch"
                << "\n";
      assert(m == num_elms);  // corrruption in bucket structure.
    }
    size_t updated = update_buckets(g, m);
    size_t num_in_range = updated - bkts[open_buckets].size();
    //none in range
    if(num_in_range == 0 && bkts[open_buckets].size() > 0) {
      auto open_bkt = bkts[open_buckets].begin();
      auto imap = parlay::delayed_seq<bucket_t>(bkts[open_buckets].size(), [&] (size_t j) {
        return (size_t)d[open_bkt[j]];
      });
      if(order == increasing) {
        size_t minBkt = parlay::reduce(imap, parlay::minm<size_t>());
        cur_range = minBkt/open_buckets-1; //will be incremented in next unpack() call
      }
      else if(order == decreasing) {
        size_t minBkt = parlay::reduce(imap, parlay::maxm<size_t>());
        cur_range = (open_buckets+minBkt)/open_buckets+1; //will be decremented in next unpack() call
      }
    }
    num_elms -= m;
    t1.stop();
  }

  inline void _next_bucket() {
    cur_bkt++;
    if (cur_bkt == open_buckets) {
      unpack();
      cur_bkt = 0;
    }
  }

  // increasing: [cur_range*open_buckets, (cur_range+1)*open_buckets)
  // decreasing: [(cur_range-1)*open_buckets, cur_range*open_buckets)
  inline bucket_id to_range(bucket_id bkt) const {
    if (order == increasing) {
      if (bkt <
          cur_range *
              open_buckets) {  // this can happen because of the lazy bucketing
        return null_bkt;
      }
      return (bkt < (cur_range + 1) * open_buckets) ? (bkt % open_buckets)
                                                    : open_buckets;
    } else {
      if (bkt >= (cur_range)*open_buckets) {
        return null_bkt;
      }
      return (bkt >= (cur_range - 1) * open_buckets)
                 ? ((open_buckets - (bkt % open_buckets)) - 1)
                 : open_buckets;
    }
  }

  size_t get_cur_bucket_num() const {
    if (order == increasing) {
      return cur_range * open_buckets + cur_bkt;
    } else {
      return (cur_range) * (open_buckets)-cur_bkt - 1;
    }
  }

  inline bucket get_cur_bucket() {
    t4.start();
    auto& bkt = bkts[cur_bkt];
    size_t size = bkt.size();
    num_elms -= size;
    size_t cur_bkt_num = get_cur_bucket_num();
    auto p = [&](size_t i) { return d[i] == cur_bkt_num; };
    auto out = parlay::filter(parlay::make_slice(bkt), p);
    bkts[cur_bkt].clear();  // reset size
    if (out.size() == 0) {
      return next_bucket();
    }
    auto ret = bucket(cur_bkt_num, std::move(out));
    ret.num_filtered = size;
    t4.stop();
    return ret;
  }

};

inline const std::optional<std::tuple<uintE, uintE> > wrap(const uintE& l,
                                                           const uintE& r) {
  if ((l != UINT_E_MAX) && (r != UINT_E_MAX)) {
    return {std::make_tuple(l, r)};
  } else {
    return std::nullopt;
  }
}

template <class ident_t, class bucket_t, class D>
inline buckets<D, ident_t, bucket_t> make_buckets(size_t n, D d, bucket_order order,
                               size_t total_buckets = 128) {
  return buckets<D, ident_t, bucket_t>(n, d, order, total_buckets);
}

// ident_t := uintE, bucket_t := uintE
template <class D>
inline buckets<D, uintE, uintE> make_vertex_buckets(size_t n, D& d, bucket_order
      order, size_t total_buckets = 128) {
  return buckets<D, uintE, uintE>(n, d, order, total_buckets);
}

// ident_t := uintE, bucket_t := bucket_t
template <class bucket_t, class D>
inline buckets<D, uintE, bucket_t> make_vertex_custom_buckets(size_t n, D& d, bucket_order
      order, size_t total_buckets = 128) {
  return buckets<D, uintE, bucket_t>(n, d, order, total_buckets);
}

}  // namespace gbbs
