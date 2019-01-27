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

#include <functional>
#include <limits>

#include "lib/sequence_ops.h"
#include "lib/extra_sequence_ops.h"
#include "maybe.h"

template <class data>
struct vertexSubsetData {
  using S = std::tuple<uintE, data>;
  using D = std::tuple<bool, data>;

  // An empty vertex set.
  vertexSubsetData(size_t _n) : n(_n), m(0), d(NULL), s(NULL), isDense(0) {}

  // A vertexSubset from array of vertex indices.
  vertexSubsetData(long _n, long _m, S* indices)
      : n(_n), m(_m), s(indices), d(NULL), isDense(0) {}

  // A vertexSubset from boolean array giving number of true values.
  vertexSubsetData(long _n, long _m, D* _d)
      : n(_n), m(_m), s(NULL), d(_d), isDense(1) {}

  // A vertexSubset from boolean array giving number of true values. Calculate
  // number of nonzeros and store in m.
  vertexSubsetData(long _n, D* _d) : n(_n), s(NULL), d(_d), isDense(1) {
    auto df = [&](size_t i) { return (size_t)std::get<0>(_d[i]); };
    auto d_map = make_sequence<size_t>(
        n, df);
    m = pbbs::reduce_add(d_map);
  }

  vertexSubsetData() : n(0), m(0), s(NULL), d(NULL), isDense(0) {}

  void del() {
    if (d != NULL) pbbs::free_array(d);
    if (s != NULL) pbbs::free_array(s);
    d = NULL;
    s = NULL;
  }

  // Sparse
  inline uintE& vtx(const uintE& i) const { return std::get<0>(s[i]); }
  inline data& vtxData(const uintE& i) const { return std::get<1>(s[i]); }
  inline std::tuple<uintE, data> vtxAndData(const uintE& i) const {
    return s[i];
  }

  // Dense
  inline bool isIn(const uintE& v) const { return std::get<0>(d[v]); }
  inline data& ithData(const uintE& v) const { return std::get<1>(d[v]); }

  // Returns (uintE) -> Maybe<std::tuple<vertex, vertex-data>>.
  auto get_fn_repr() const
      -> std::function<Maybe<std::tuple<uintE, data>>(uintE)> {
    std::function<Maybe<std::tuple<uintE, data>>(const uintE&)> fn;
    if (isDense) {
      fn = [&](const uintE& v) -> Maybe<std::tuple<uintE, data>> {
        auto ret = Maybe<std::tuple<uintE, data>>(
            std::make_tuple(v, std::get<1>(d[v])));
        ret.exists = std::get<0>(d[v]);
        return ret;
      };
    } else {
      fn = [&](const uintE& i) -> Maybe<std::tuple<uintE, data>> {
        return Maybe<std::tuple<uintE, data>>(s[i]);
      };
    }
    return fn;
  }

  long size() { return m; }
  long numVertices() { return n; }

  long numRows() { return n; }
  long numNonzeros() { return m; }

  bool isEmpty() { return m == 0; }
  bool dense() { return isDense; }

  void toSparse() {
    if (s == NULL && m > 0) {
      auto f = [&](size_t i) -> std::tuple<bool, data> { return d[i]; };
      auto f_seq = make_sequence<D>(n, f);
      auto out = pbbs::pack_index_and_data<uintE, data>(f_seq, n);
      s = out.get_array();
      if (out.size() != m) {
        std::cout << "bad stored value of m"
                  << "\n";
        abort();
      }
    }
    isDense = false;
  }

  // Convert to dense but keep sparse representation if it exists.
  void toDense() {
    if (d == NULL) {
      d = pbbs::new_array_no_init<D>(n);
      par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                      { std::get<0>(d[i]) = false; });
      par_for(0, m, pbbs::kSequentialForThreshold, [&] (size_t i) {
        d[std::get<0>(s[i])] = std::make_tuple(true, std::get<1>(s[i]));
      });
    }
    isDense = true;
  }

  size_t n, m;
  S* s;
  D* d;
  bool isDense;
};

// Specialized version where data = pbbs::empty.
template <>
struct vertexSubsetData<pbbs::empty> {
  using S = uintE;

  // An empty vertex set.
  vertexSubsetData<pbbs::empty>(size_t _n)
      : n(_n), m(0), s(NULL), d(NULL), isDense(0) {}

  // A vertexSubset with a single vertex.
  vertexSubsetData<pbbs::empty>(long _n, uintE v)
      : n(_n), m(1), d(NULL), isDense(0) {
    s = pbbs::new_array_no_init<uintE>(1);
    s[0] = v;
  }

  // A vertexSubset from array of vertex indices.
  vertexSubsetData<pbbs::empty>(long _n, long _m, S* indices)
      : n(_n), m(_m), s(indices), d(NULL), isDense(0) {}

  // A vertexSubset from array of vertex indices.
  vertexSubsetData<pbbs::empty>(long _n, long _m,
                                std::tuple<uintE, pbbs::empty>* indices)
      : n(_n), m(_m), s((uintE*)indices), d(NULL), isDense(0) {}

  // A vertexSubset from boolean array giving number of true values.
  vertexSubsetData<pbbs::empty>(long _n, long _m, bool* _d)
      : n(_n), m(_m), s(NULL), d(_d), isDense(1) {}

  // A vertexSubset from boolean array giving number of true values. Calculate
  // number of nonzeros and store in m.
  vertexSubsetData<pbbs::empty>(long _n, bool* _d)
      : n(_n), s(NULL), d(_d), isDense(1) {
    auto d_f = [&](size_t i) { return _d[i]; };
    auto d_map = make_sequence<size_t>(n, d_f);
    auto f = [&](size_t i, size_t j) { return i + j; };
    m = pbbs::reduce(d_map, f);
  }

  // A vertexSubset from boolean array giving number of true values. Calculate
  // number of nonzeros and store in m.
  vertexSubsetData<pbbs::empty>(long _n, std::tuple<bool, pbbs::empty>* _d)
      : n(_n), s(NULL), d((bool*)_d), isDense(1) {
    auto d_f = [&](size_t i) { return std::get<0>(_d[i]); };
    auto d_map = make_sequence<size_t>(n, d_f);
    auto f = [&](size_t i, size_t j) { return i + j; };
    m = pbbs::reduce(d_map, f);
  }

  void del() {
    if (d != NULL) {
      pbbs::free_array(d);
    }
    if (s != NULL) {
      pbbs::free_array(s);
    }
    d = NULL;
    s = NULL;
  }

  // Sparse
  inline uintE& vtx(const uintE& i) const { return s[i]; }
  inline pbbs::empty vtxData(const uintE& i) const { return pbbs::empty(); }
  inline std::tuple<uintE, pbbs::empty> vtxAndData(const uintE& i) const {
    return std::make_tuple(s[i], pbbs::empty());
  }

  // Dense
  inline bool isIn(const uintE& v) const { return d[v]; }
  inline pbbs::empty ithData(const uintE& v) const { return pbbs::empty(); }

  // Returns (uintE) -> Maybe<std::tuple<vertex, vertex-data>>.
  auto get_fn_repr() const
      -> std::function<Maybe<std::tuple<uintE, pbbs::empty>>(uintE)> {
    std::function<Maybe<std::tuple<uintE, pbbs::empty>>(const uintE&)> fn;
    if (isDense) {
      fn = [&](const uintE& v) -> Maybe<std::tuple<uintE, pbbs::empty>> {
        auto ret = Maybe<std::tuple<uintE, pbbs::empty>>(
            std::make_tuple(v, pbbs::empty()));
        ret.exists = d[v];
        return ret;
      };
    } else {
      fn = [&](const uintE& i) -> Maybe<std::tuple<uintE, pbbs::empty>> {
        return Maybe<std::tuple<uintE, pbbs::empty>>(
            std::make_tuple(s[i], pbbs::empty()));
      };
    }
    return fn;
  }

  long size() { return m; }
  long numVertices() { return n; }

  long numRows() { return n; }
  long numNonzeros() { return m; }

  bool isEmpty() { return m == 0; }
  bool dense() { return isDense; }

  void toSparse() {
    if (s == NULL && m > 0) {
      auto _d = d;
      auto f = [&](size_t i) { return _d[i]; };
      auto f_in = make_sequence<bool>(n, f);
      auto out = pbbs::pack_index<uintE>(f_in);
      s = out.get_array();
      if (out.size() != m) {
        std::cout << "bad stored value of m"
                  << "\n";
        std::cout << "out.size = " << out.size() << " m = " << m << " n = " << n
                  << "\n";
        abort();
      }
    }
    isDense = false;
  }

  // Converts to dense but keeps sparse representation if it exists.
  void toDense() {
    if (d == NULL) {
      d = pbbs::new_array_no_init<bool>(n);
      par_for(0, n, pbbs::kSequentialForThreshold, [&] (size_t i)
                      { d[i] = 0; });
      par_for(0, m, pbbs::kSequentialForThreshold, [&] (size_t i)
                      { d[s[i]] = 1; });
    }
    isDense = true;
  }

  size_t n, m;
  S* s;
  bool* d;
  bool isDense;
};

using vertexSubset = vertexSubsetData<pbbs::empty>;
