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

// Usage:
// numactl -i all ./Triangle -rounds 2 -s -c -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm

#include "Triangle.h"

template <class G>
double Triangle_runner(G& GA, commandLine P) {
  std::cout << "### Application: Triangle Counting" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: n/a" << std::endl;
  std::cout << "### ------------------------------------" << endl;
  assert(P.getOption("-s"));
  size_t count = 0;
  auto f = [&] (uintE u, uintE v, uintE w) { };
  timer t; t.start();
  count = Triangle(GA, f);
  double tt = t.stop();
  if (P.getOption("-stats")) {
    auto wedge_im_f = [&](size_t i) {
      size_t deg = GA.V[i].getOutDegree();
      return (deg * deg - 1) / 2;
    };
    auto wedge_im = pbbslib::make_sequence<size_t>(GA.n, wedge_im_f);
    size_t n_wedges = pbbslib::reduce_add(wedge_im);
    std::cout << "### n_wedges = " << n_wedges << "\n";
    std::cout << "### triangle density = " << ((3.0 * count) / n_wedges) << "\n";
  }

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

generate_main(Triangle_runner, false);
