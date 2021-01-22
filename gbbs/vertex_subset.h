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

#include "bridge.h"
#include "flags.h"
#include "macros.h"

#include <functional>
#include <limits>
#include <optional>

namespace gbbs {

template <class data>
struct vertexSubsetData {
  using S = std::tuple<uintE, data>;
  using D = std::tuple<bool, data>;

  // An empty vertex set.
  vertexSubsetData(size_t _n) : n(_n), m(0), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  // A vertexSubset from array of vertex indices.
  vertexSubsetData(size_t _n, size_t _m, parlay::sequence<S>&& indices)
      : n(_n), m(_m), s(std::move(indices)), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

//  // A vertexSubset from a sequence.
//  vertexSubsetData(size_t _n, sequence<S>& seq, bool transfer = true)
//      : n(_n), m(seq.size()), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {
//    if (transfer) {
//      s = seq.to_array();
//    } else {
//      s = seq.begin();
//    }
//  }

  // A vertexSubset from boolean array giving number of true values.
  vertexSubsetData(size_t _n, size_t _m, parlay::sequence<D>&& _d)
      : n(_n), m(_m), d(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  // A vertexSubset from boolean array giving number of true values. Calculate
  // number of nonzeros and store in m.
  vertexSubsetData(size_t _n, parlay::sequence<D>&& _d) : n(_n), d(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    auto d_map = parlay::delayed_seq<size_t>(n, [&](size_t i) { return (size_t)std::get<0>(_d[i]); });
    m = pbbslib::reduce_add(d_map);
  }

  vertexSubsetData() : n(0), m(0), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  bool out_degrees_set() {
    return (sum_out_degrees != std::numeric_limits<size_t>::max());
  }
  size_t get_out_degrees() {
    return sum_out_degrees;
  }
  void set_out_degrees(size_t _sum_out_degrees) {
    sum_out_degrees = _sum_out_degrees;
  }

  // Sparse
  inline uintE& vtx(const uintE& i) const { return std::get<0>(s[i]); }
  inline data& vtxData(const uintE& i) const { return std::get<1>(s[i]); }
  inline std::tuple<uintE, data> vtxAndData(const uintE& i) const {
    return s[i];
  }

  // Dense
  __attribute__((always_inline)) inline bool isIn(const uintE& v) const {
    return std::get<0>(d.begin()[v]);
  }
  inline data& ithData(const uintE& v) const { return std::get<1>(d[v]); }

  // Returns (uintE) -> std::optional<std::tuple<vertex, vertex-data>>.
  auto get_fn_repr() const
      -> std::function<std::optional<std::tuple<uintE, data>>(uintE)> {
    std::function<std::optional<std::tuple<uintE, data>>(const uintE&)> fn;
    if (isDense) {
      fn = [&](const uintE& v) -> std::optional<std::tuple<uintE, data>> {
        const auto& dv = d[v];
        if (std::get<0>(dv)) {
          return std::optional<std::tuple<uintE, data>>(std::make_tuple(v, std::get<1>(d[v])));
        } else {
          return std::nullopt;
        }
      };
    } else {
      fn = [&](const uintE& i) -> std::optional<std::tuple<uintE, data>> {
        return std::optional<std::tuple<uintE, data>>(s[i]);
      };
    }
    return fn;
  }

  size_t size() const { return m; }
  size_t numVertices() const { return n; }

  size_t numRows() const { return n; }
  size_t numNonzeros() const { return m; }

  bool isEmpty() const { return m == 0; }
  bool dense() const { return isDense; }

  void toSparse() {
    if (s.size() == 0 && m > 0) {
      auto f_seq = parlay::delayed_seq<D>(n, [&](size_t i) -> std::tuple<bool, data> { return d[i]; });
      auto out = pbbslib::pack_index_and_data<uintE, data>(f_seq, n);
      if (out.size() != m) {
        std::cout << "# m is " << m << " but out.size says" << out.size() << std::endl;
        std::cout << "# bad stored value of m"
                  << "\n";
        abort();
      }
      s = out.to_array();
    }
    isDense = false;
  }

  // Convert to dense but keep sparse representation if it exists.
  void toDense() {
    if (d.size() == 0) {
      d = parlay::sequence<D>(n);
      par_for(0, n, [&](size_t i) { std::get<0>(d[i]) = false; });
      par_for(0, m, [&](size_t i) {
        d[std::get<0>(s[i])] = std::make_tuple(true, std::get<1>(s[i]));
      });
    }
    isDense = true;
  }

  size_t n, m;
  parlay::sequence<S> s;
  parlay::sequence<D> d;
  bool isDense;
  size_t sum_out_degrees;
};

// Specialized version where data = gbbs::empty.
template <>
struct vertexSubsetData<gbbs::empty> {
  using S = uintE;
  using D = bool;

  // An empty vertex set.
  vertexSubsetData<gbbs::empty>(size_t _n)
      : n(_n), m(0), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  // A vertexSubset with a single vertex.
  vertexSubsetData<gbbs::empty>(size_t _n, uintE v)
      : n(_n), m(1), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    s = parlay::sequence<uintE>(1);
    s[0] = v;
  }

  // A vertexSubset from array of vertex indices.
  vertexSubsetData<gbbs::empty>(size_t _n, size_t _m, parlay::sequence<S>&& indices)
      : n(_n), m(_m), s(std::move(indices)), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

//  // A vertexSubset from a sequence.
//  vertexSubsetData<gbbs::empty>(size_t _n, sequence<S>& seq, bool transfer =
//  true)
//      : n(_n), m(seq.size()), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {
//    if (transfer) {
//      s = seq.to_array();
//    } else {
//      s = seq.begin();
//    }
//  }

  // A vertexSubset from boolean array giving number of true values.
  vertexSubsetData<gbbs::empty>(size_t _n, size_t _m, parlay::sequence<bool>&& _d)
      : n(_n), m(_m), d(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  // A vertexSubset from boolean array giving number of true values. Calculate
  // number of nonzeros and store in m.
  vertexSubsetData<gbbs::empty>(size_t _n, parlay::sequence<bool>&& _d)
      : n(_n), d(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    auto d_map = parlay::delayed_seq<size_t>(n, [&](size_t i) { return d[i]; });
    m = pbbslib::reduce_add(d_map);
  }

  bool out_degrees_set() {
    return (sum_out_degrees != std::numeric_limits<size_t>::max());
  }
  size_t get_out_degrees() {
    return sum_out_degrees;
  }
  void set_out_degrees(size_t _sum_out_degrees) {
    sum_out_degrees = _sum_out_degrees;
  }

  // Sparse
  inline uintE vtx(const uintE& i) const { return s[i]; }

  inline gbbs::empty vtxData(const uintE& i) const {
    return gbbs::empty();
  }
  inline std::tuple<uintE, gbbs::empty> vtxAndData(const uintE& i) const {
    return std::make_tuple(s[i], gbbs::empty());
  }

  // Dense
  __attribute__((always_inline)) inline bool isIn(const uintE& v) const {
    return d.begin()[v];
  }
  inline gbbs::empty ithData(const uintE& v) const {
    return gbbs::empty();
  }

  // Returns (uintE) -> std::optional<std::tuple<vertex, vertex-data>>.
  auto get_fn_repr() const
      -> std::function<std::optional<std::tuple<uintE, gbbs::empty>>(uintE)> {
    std::function<std::optional<std::tuple<uintE, gbbs::empty>>(const uintE&)> fn;
    if (isDense) {
      fn = [&](const uintE& v) -> std::optional<std::tuple<uintE, gbbs::empty>> {
        if (d[v]) {
          return std::optional<std::tuple<uintE, gbbs::empty>>(std::make_tuple(v, gbbs::empty()));
        } else {
          return std::nullopt;
        }
      };
    } else {
      fn = [&](const uintE& i) -> std::optional<std::tuple<uintE, gbbs::empty>> {
        return std::optional<std::tuple<uintE, gbbs::empty>>(
            std::make_tuple(s[i], gbbs::empty()));
      };
    }
    return fn;
  }

  size_t size() const { return m; }
  size_t numVertices() const { return n; }

  size_t numRows() const { return n; }
  size_t numNonzeros() const { return m; }

  bool isEmpty() const { return m == 0; }
  bool dense() const { return isDense; }

  void toSparse() {
    if (s.size() == 0 && m > 0) {
      auto _d = d;
      auto f_in = parlay::delayed_seq<bool>(n, [&](size_t i) { return _d[i]; });
      s = parlay::pack_index<uintE>(f_in);
      if (s.size() != m) {
        std::cout << "# m is " << m << " but s.size says" << s.size() << std::endl;
        std::cout << "# bad stored value of m"
                  << "\n";
        std::cout << "# s.size = " << s.size() << " m = " << m << " n = " << n
                  << "\n";
        abort();
      }
    }
    isDense = false;
  }

  // Converts to dense but keeps sparse representation if it exists.
  void toDense() {
    if (d.size() == 0) {
      d = parlay::sequence<bool>(n);
      par_for(0, n, [&](size_t i) { d[i] = 0; });
      par_for(0, m, [&](size_t i) { d[s[i]] = 1; });
    }
    isDense = true;
  }

  size_t n, m;
  parlay::sequence<S> s;
  parlay::sequence<bool> d;
  bool isDense;
  size_t sum_out_degrees;
};

using vertexSubset = vertexSubsetData<gbbs::empty>;

/* ======================== Functions on VertexSubsets ====================== */

// Takes a vertexSubsetData (with some non-trivial Data) and applies a map
// function f : (uintE x Data) -> void over each vertex in the vertexSubset, in
// parallel.
template <class F, class VS,
          typename std::enable_if<!std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f, size_t granularity=gbbs::kSequentialForThreshold) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    parallel_for(0, n, [&] (size_t i) {
      if (V.isIn(i)) {
        f(i, V.ithData(i));
      }
    }, granularity);
  } else {
    parallel_for(0, m, [&] (size_t i) { f(V.vtx(i), V.vtxData(i)); }, granularity);
  }
}

// Takes a vertexSubset (with no extra data per-vertex) and applies a map
// function f : uintE -> void over each vertex in the vertexSubset, in
// parallel.
template <class VS, class F,
          typename std::enable_if<std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f, size_t granularity=gbbs::kSequentialForThreshold) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    parallel_for(0, n, [&] (size_t i) {
      if (V.isIn(i)) {
        f(i);
      }
    }, granularity);
  } else {
    parallel_for(0, m, [&] (size_t i)
                    { f(V.vtx(i)); }, granularity);
  }
}

template <class F, class Data>
inline vertexSubset vertexFilter_dense(vertexSubsetData<Data>& V, F filter, size_t granularity=gbbs::kSequentialForThreshold) {
  size_t n = V.numRows();
  V.toDense();
  auto d_out = parlay::sequence<bool>(n);
  parallel_for(0, n, [&] (size_t i) { d_out[i] = 0; }, granularity);
  parallel_for(0, n, [&] (size_t i) {
    if constexpr (std::is_same<Data, gbbs::empty>::value) {
      if (V.isIn(i)) d_out[i] = filter(i);
    } else {
      if (V.isIn(i)) d_out[i] = filter(i, V.ithData(i));
    }
  }, granularity);
  return vertexSubset(n, std::move(d_out));
}

template <class F, class Data>
inline vertexSubset vertexFilter_sparse(vertexSubsetData<Data>& V, F filter, size_t granularity=gbbs::kSequentialForThreshold) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  auto bits = parlay::sequence<bool>(m);
  V.toSparse();
  parallel_for(0, m, [&] (size_t i) {
    uintE v = V.vtx(i);
    if constexpr (std::is_same<Data, gbbs::empty>::value) {
      bits[i] = filter(v);
    } else {
      bits[i] = filter(v, V.vtxData(i));
    }
  }, granularity);
  auto v_imap = parlay::delayed_seq<uintE>(m, [&](size_t i) { return V.vtx(i); });
  auto bits_m = parlay::delayed_seq<bool>(m, [&](size_t i) { return bits[i]; });
  auto out = parlay::pack(v_imap, bits_m);
  size_t out_size = out.size();
  return vertexSubset(n, std::move(out));
}

// Note that this currently strips vertices of their associated data (which is
// the intended use-case in all current uses). Should refactor at some point to
// make keeping/removing the data a choice.
template <class F, class VS>
inline vertexSubset vertexFilter(VS& vs, F filter, flags fl = 0) {
  if (fl == dense_only) {
    return vertexFilter_dense(vs, filter);
  } else if (fl == no_dense) {
    return vertexFilter_sparse(vs, filter);
  }
  // TODO: can measure selectivity and call sparse/dense based on a sample.
  if (vs.dense()) {
    return vertexFilter_dense(vs, filter);
  }
  return vertexFilter_sparse(vs, filter);
}

// TODO(laxman): update interface
//void add_to_vsubset(vertexSubset& vs, uintE* new_verts, uintE num_new_verts);

}  // namespace gbbs
