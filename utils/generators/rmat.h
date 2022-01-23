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

#include "gbbs/bridge.h"
#include "gbbs/macros.h"

namespace gbbs {
namespace rmat {
double hashDouble(uintE i) {
  return ((double)(parlay::hash32(i)) / ((double)UINT_E_MAX));
}

struct rMat_gen {
  double a, ab, abc;
  uintE n;
  size_t h;
  rMat_gen(uintE _n, uintE _seed, double _a, double _b, double _c) {
    n = _n;
    a = _a;
    ab = _a + _b;
    abc = _a + _b + _c;
    h = parlay::hash64(_seed);
    if (abc > 1) {
      std::cout << "in rMat: a + b + c add to more than 1\n";
      abort();
    }
    if ((1U << parlay::log2_up(n)) != n) {
      std::cout << "in rMat: n not a power of 2";
      exit(0);
    }
  }

  std::tuple<uintE, uintE> rMatRec(uintE nn, uintE randStart,
                                   uintE randStride) {
    if (nn == 1)
      return std::make_tuple((uintE)0, (uintE)0);
    else {
      auto x = rMatRec(nn / 2, randStart + randStride, randStride);
      double r = hashDouble(randStart);
      if (r < a)
        return x;
      else if (r < ab)
        return std::make_tuple(std::get<0>(x), std::get<1>(x) + nn / 2);
      else if (r < abc)
        return std::make_tuple(std::get<0>(x) + nn / 2, std::get<1>(x));
      else
        return std::make_tuple(std::get<0>(x) + nn / 2,
                               std::get<1>(x) + nn / 2);
    }
  }

  std::tuple<uintE, uintE> operator()(uintE i) {
    uintE randStart = (uintE)parlay::hash64((uintE)(2 * i) * h);
    uintE randStride = (uintE)parlay::hash64((uintE)(2 * i + 1) * h);
    return rMatRec(n, randStart, randStride);
  }
};

parlay::sequence<std::tuple<uintE, uintE>> generate_updates(uintE n, size_t m,
                                                            uintE seed,
                                                            float a = 0.5,
                                                            float b = 0.1,
                                                            float c = 0.1) {
  uintE nn = (1 << parlay::log2_up(n));
  rMat_gen g(nn, seed, a, b, c);
  auto E = parlay::sequence<std::tuple<uintE, uintE>>(m);
  parallel_for(0, m, [&](size_t i) { E[i] = g(i); });
  return E;
}
}  // namespace rmat
}  // namespace gbbs
