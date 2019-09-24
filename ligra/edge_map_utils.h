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

// Contains helper functions and special cases of edgeMap. Most of these
// functions are specialized for the type of data written per vertex using tools
// from type_traits.

#pragma once

#include <functional>
#include <tuple>
#include <type_traits>

#include "pbbslib/binary_search.h"
#include "pbbslib/utilities.h"
#include "maybe.h"

// Standard version of edgeMapDense.
template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_gen(std::tuple<bool, data>* next) {
  return [next](uintE ngh, bool m = false) __attribute__((always_inline)) {
    if (m) next[ngh] = std::make_tuple(1, pbbslib::empty());
  };
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_gen(std::tuple<bool, data>* next) {
  return [next](uintE ngh, Maybe<data> m = Maybe<data>()) __attribute__((always_inline)) {
    if (m.exists) next[ngh] = std::make_tuple(1, m.t);
  };
}

// Standard version of edgeMapDenseForward.
template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_forward_gen(std::tuple<bool, data>* next) {
  return [next](uintE ngh, bool m = false) __attribute__((always_inline)) {
    if (m) next[ngh] = std::make_tuple(1, pbbslib::empty());
  };
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_forward_gen(std::tuple<bool, data>* next) {
  return [next](uintE ngh, Maybe<data> m = Maybe<data>()) __attribute__((always_inline)) {
    if (m.exists) next[ngh] = std::make_tuple(1, m.t);
  };
}

// Standard version of edgeMapSparse.
template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emsparse_gen_full(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset, bool m) __attribute__((always_inline)) {
    if (m) {
      outEdges[offset] = std::make_tuple(ngh, pbbslib::empty());
    } else {
      outEdges[offset] = std::make_tuple(UINT_E_MAX, pbbslib::empty());
    }
  };
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emsparse_gen_full(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset, Maybe<data> m = Maybe<data>()) __attribute__((always_inline)) {
    if (m.exists) {
      outEdges[offset] = std::make_tuple(ngh, m.t);
    } else {
      std::get<0>(outEdges[offset]) = UINT_E_MAX;
    }
  };
}

template <typename data>
inline auto get_emsparse_gen_empty(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset) __attribute__((always_inline)) {
    std::get<0>(outEdges[offset]) = UINT_E_MAX;
  };
}

// edgeMapSparse_no_filter
// Version of edgeMapSparse that binary-searches and packs out blocks of the
// next frontier.
template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emsparse_blocked_gen(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset, bool m = false) __attribute__((always_inline)) {
    if (m) {
      outEdges[offset] = std::make_tuple(ngh, pbbslib::empty());
      return true;
    }
    return false;
  };
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emsparse_blocked_gen(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset, Maybe<data> m = Maybe<data>()) __attribute__((always_inline)) {
    if (m.exists) {
      outEdges[offset] = std::make_tuple(ngh, m.t);
      return true;
    }
    return false;
  };
}

template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emblock_gen(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset, bool m = false)
      __attribute__((always_inline)) {
    if (m) {
      outEdges[offset] = std::make_tuple(ngh, pbbslib::empty());
      return true;
    }
    return false;
  };
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emblock_gen(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset, Maybe<data> m = Maybe<data>())
      __attribute__((always_inline)) {
    if (m.exists) {
      outEdges[offset] = std::make_tuple(ngh, m.t);
      return true;
    }
    return false;
  };
}

// Gen-functions that produce no output
template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emsparse_nooutput_gen() {
  return [&](uintE ngh, uintE offset, bool m = false) {};
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emsparse_nooutput_gen() {
  return [&](uintE ngh, uintE offset, Maybe<data> m = Maybe<data>()) {};
}

template <typename data>
inline auto get_emsparse_nooutput_gen_empty() {
  return [&](uintE ngh, uintE offset) {};
}

template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_nooutput_gen() {
  return [&](uintE ngh, bool m = false) {};
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_nooutput_gen() {
  return [&](uintE ngh, Maybe<data> m = Maybe<data>()) {};
}

template <typename data,
          typename std::enable_if<std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_forward_nooutput_gen() {
  return [&](uintE ngh, bool m = false) {};
}

template <typename data,
          typename std::enable_if<!std::is_same<data, pbbslib::empty>::value,
                                  int>::type = 0>
inline auto get_emdense_forward_nooutput_gen() {
  return [&](uintE ngh, Maybe<data> m = Maybe<data>()) {};
}

template <class W, class F>
struct Wrap_F {
  F f;
  Wrap_F(F _f) : f(_f) {}
  inline bool update(const uintE& s, const uintE& d, const W& e) {
    return f.update(s, d);
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& e) {
    return f.updateAtomic(s, d);
  }
  inline bool cond(const uintE& d) { return f.cond(d); }
};

template <class W, class D, class F>
struct Wrap_Default_F {
  F f;
  D def;
  Wrap_Default_F(F _f, D _def) : f(_f), def(_def) {}
  inline bool update(const uintE& s, const uintE& d, const W& e) {
    return f.update(s, d, def);
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& e) {
    return f.updateAtomic(s, d, def);
  }
  inline bool cond(const uintE& d) { return f.cond(d); }
};

template <class W, class F>
inline auto wrap_em_f(F f) -> Wrap_F<W, F> {
  return Wrap_F<W, F>(f);
}

template <class W, class D, class F,
          typename std::enable_if<!std::is_same<W, D>::value, int>::type = 0>
inline auto wrap_with_default(F f, D def) -> Wrap_Default_F<W, D, F> {
  return Wrap_Default_F<W, D, F>(f, def);
}

template <class W, class D, class F,
          typename std::enable_if<std::is_same<W, D>::value, int>::type = 0>
inline auto wrap_with_default(F f, D def) -> decltype(f) {
  return f;
}
