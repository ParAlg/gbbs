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
// numactl -i all ./CC -rounds 3 -s -m twitter_SJ
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -stats : print the #ccs, and the #vertices in the largest cc

#include "Connectivity.h"
#include "union_find_rules.h"

template <class Seq>
inline size_t num_cc(Seq& labels) {
  size_t n = labels.size();
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    if (!flags[labels[i]]) {
      flags[labels[i]] = 1;
    }
  });
  pbbslib::scan_add_inplace(flags);
  std::cout << "n_cc = " << flags[n] << "\n";
  return flags[n];
}

template <class Seq>
inline size_t largest_cc(Seq& labels) {
  size_t n = labels.size();
  // could histogram to do this in parallel.
  auto flags = sequence<uintE>(n + 1, [&](size_t i) { return 0; });
  for (size_t i = 0; i < n; i++) {
    flags[labels[i]] += 1;
  }
  size_t sz = pbbslib::reduce_max(flags);
  std::cout << "largest connected component has size: " << sz << "\n";
  return sz;
}

template <class W>
double CC_runner(edge_array<W>& G, commandLine P) {
  std::cout << "### Application: CC (Connectivity)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.num_rows << std::endl;
  std::cout << "### m: " << G.non_zeros << std::endl;
  std::cout << "### ------------------------------------" << endl;

  std::string test = P.getOptionValue("-testname", "async");

  timer t;
  t.start();
  pbbs::sequence<uintE> components;
  if (test == "async_hook") {
    auto find = find_variants::find_compress;
    auto unite = unite_variants::Unite<decltype(find)>(find);
    components = union_find::UnionFindHookTemplate_coo(G, unite, find);
  } else if (test == "async") {
    auto find = find_variants::find_compress;
    auto unite = unite_variants::Unite<decltype(find)>(find);
    components = union_find::UnionFindTemplate_coo(G, unite, find);
  } else if (test == "rem") {
    auto find = find_variants::find_compress;
    auto unite = unite_variants::UniteRem<decltype(find)>(find, G.num_rows);
    components = union_find::UnionFindTemplate_coo(G, unite, find);
  } else {
    cout << "Unknown test named: " << test << endl;
  }

  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  if (P.getOption("-stats")) {
    auto cc_f = [&](size_t i) { return components[i]; };
    auto cc_im =
        pbbslib::make_sequence<uintE>(G.num_rows, cc_f);
    num_cc(cc_im);
    largest_cc(cc_im);
  }
  return tt;
}

generate_coo_main(CC_runner, false);
