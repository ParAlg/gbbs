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
// numactl -i all ./MIS -stats -rounds 4 -s com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc
//     -specfor : run the speculative_for based algorithm from pbbs

#include "MIS.h"
#include "ligra.h"

template <class vertex>
double MIS_runner(graph<vertex>& GA, commandLine P) {
  bool spec_for = P.getOption("-specfor");
  std::cout << "### Application: MIS" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: -specfor (deterministic reservations) = " << spec_for << std::endl;

  assert(P.getOption("-s"));
  double tt = 0.0;

  // Code below looks duplicated; this is because the return types of specfor
  // and rootset are different
  if (spec_for) {
    timer t; t.start();
    auto MIS = MIS_spec_for::MIS(GA);
    // in spec_for, MIS[i] == 1 indicates that i was chosen
    tt = t.stop();
    auto size_f = [&](size_t i) { return (MIS[i] == 1); };
    auto size_imap =
        pbbslib::make_sequence<size_t>(GA.n, size_f);
    if (P.getOptionValue("-stats")) {
      std::cout << "MIS size: " << pbbslib::reduce_add(size_imap) << "\n";
    }
    if (P.getOptionValue("-verify")) {
      verify_MIS(GA, size_imap);
    }
  } else {
    timer t; t.start();
    auto MIS = MIS_rootset::MIS(GA);
    tt = t.stop();
    auto size_f = [&](size_t i) { return MIS[i]; };
    auto size_imap =
        pbbslib::make_sequence<size_t>(GA.n, size_f);
    if (P.getOptionValue("-stats")) {
      std::cout << "MIS size: " << pbbslib::reduce_add(size_imap) << "\n";
    }
    if (P.getOptionValue("-verify")) {
      verify_MIS(GA, size_imap);
    }
  }

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

generate_main(MIS_runner, false);
