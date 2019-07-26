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
// numactl -i all ./MaximalMatching -s -c -m clueweb_sym.bytepda
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "MaximalMatching.h"

#include "bridge.h"
#include "ligra.h"


#include <fstream>
#include <iostream>

template <class G>
double MaximalMatching_runner(G& GA, commandLine P) {
  std::cout << "### Application: CC (Connectivity)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: (n/a)" << std::endl;
  std::cout << "### ------------------------------------" << endl;

  assert(P.getOption("-s"));  // input graph must be symmetric
//  auto in_f = P.getOptionValue("-if");
//  if (in_f) {
//    auto S = readStringFromFile(in_f);
//    auto W = pbbslib::tokenize(S, [] (const char c) { return pbbs::is_space(c); });
//    size_t ms = atol(W[0]);
//    using edge = std::tuple<uintE, uintE>;
//    auto matching = sequence<edge>(ms);
//    par_for(0, ms, pbbslib::kSequentialForThreshold, [&] (size_t i) {
//      matching[i] =
//          std::make_tuple(atol(W[1 + 2 * i]), atol(W[2 * (i + 1)]));
//    });
//    verify_matching(GA, matching);
//    exit(0);
//  }
  timer t; t.start();
  auto matching = MaximalMatching(GA);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;

  // Note that as we mutate the graph by packing out edges, we can't verify that
  // the returned set of edges is a valid maximal matching on the graph
  // currently in-memory. Instead, we write the matching to disk and read it
  // back in to verify correctness.
  auto of = P.getOptionValue("-of");
  if (of) {
    std::cout << "outfile is = " << of << "\n";
    std::ofstream out(of, std::ofstream::out);
    out << matching.size() << "\n";
    for (size_t i = 0; i < matching.size(); i++) {
      auto e = matching[i];
      out << (std::get<0>(e) & mm::VAL_MASK) << " " << std::get<1>(e) << "\n";
    }
    out.close();
  }
  // Maximal-matching mutates the underlying graph (unless it is copied, which
  // we don't do to prevent memory issues), so we make sure the algorithm is run
  // exactly once.
  exit(0);
}

generate_main(MaximalMatching_runner, true);
