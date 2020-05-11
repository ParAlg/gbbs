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

#include <tuple>
#include "macros.h"

template <class T>
struct Maybe {
  // pays an extra byte for NONE/SOME
  // avoid Maybe<T>*---this is basically free otherwise
  bool exists;
  T t;
  Maybe(const T& _t) : exists(true), t(_t) {}
  Maybe() : exists(false) {}
};

inline const Maybe<std::tuple<uintE, uintE> > wrap(const uintE& l,
                                                   const uintE& r) {
  auto t = Maybe<std::tuple<uintE, uintE> >(std::make_tuple(l, r));
  t.exists = (l != UINT_E_MAX) && (r != UINT_E_MAX);
  return t;
}

// integer types L and R
template <class L, class R>
inline const Maybe<std::tuple<L, R> > wrap(const L& l, const R& r) {
  auto t = Maybe<std::tuple<L, R> >(std::make_tuple(l, r));
  t.exists = (l != std::numeric_limits<L>::max()) && (r != std::numeric_limits<R>::max());
  return t;
}


template <class L, class R>
inline const Maybe<std::tuple<L, R> > wrap(const L& l, const Maybe<R>& r) {
  auto t = Maybe<std::tuple<L, R> >(std::make_tuple(l, getT(r)));
  t.exists = r.exists;
  return t;
}

template <class L, class R>
inline Maybe<std::tuple<L, R> > wrap(const Maybe<L>& l, const R& r) {
  auto t = Maybe<std::tuple<L, R> >(std::make_tuple(getT(l), r));
  t.exists = l.exists;
  return t;
}

template <class L, class R>
inline Maybe<std::tuple<L, R> > wrap(const Maybe<L>& l, const Maybe<R>& r) {
  auto t = Maybe<std::tuple<L, R> >(std::make_tuple(getT(l), getT(r)));
  t.exists = l.exists && r.exists;
  return t;
}

template <class T>
inline bool isSome(const Maybe<T>& m) {
  return m.exists;
}

template <class T>
inline T getT(const Maybe<T>& m) {
  return m.t;
}
