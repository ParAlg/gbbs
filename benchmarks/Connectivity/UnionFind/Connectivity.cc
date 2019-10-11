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
#include "jayanti.h"
#include "union_find_rules.h"
#include "benchmarks/Connectivity/WorkEfficientSDB14/Connectivity.h"

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




template <class Seq>
inline size_t RelabelDet(Seq& ids) {
  using T = typename Seq::value_type;
  size_t n = ids.size();
  auto component_map = pbbs::sequence<T>(n + 1, (T)0);
  T cur_comp = 1;
  for (size_t i=0; i<n; i++) {
    T comp = ids[i];
    T new_comp = cur_comp++;
    if (component_map[comp] == 0) {
      component_map[comp] = new_comp;
    }
    ids[i] = new_comp;
  }
  return cur_comp;
//  pbbslib::scan_add_inplace(inverse_map);
//
//  size_t new_n = inverse_map[n];
//  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
//                  { ids[i] = inverse_map[ids[i]]; });
//  return new_n;
}

template <class S1, class S2>
inline void cc_check(S1& correct, S2& check) {
  RelabelDet(check);

  bool is_correct = true;
  uintE max_cor = 0;
  uintE max_chk = 0;
//  parallel_for(0, correct.size(), [&] (size_t i) {
  for (size_t i=0; i<correct.size(); i++) {
    assert(correct[i] == check[i]);
    if ((correct[i] != check[i])) {
      is_correct = false;
      cout << "at i = " << i << " cor = " << correct[i] << " got: " << check[i] << endl;
    }
    if (correct[i] > max_cor) {
      pbbs::write_max(&max_cor, correct[i], std::less<uintE>());
    }
    if (check[i] > max_chk) {
      pbbs::write_max(&max_chk, check[i], std::less<uintE>());
    }
  }//);
  cout << "correctness check: " << is_correct << endl;
  cout << "max_cor = " << max_cor << " max_chk = " << max_chk << endl;
}

template <class Graph>
double CC_runner(Graph& G, commandLine P) {
  std::cout << "### Application: CC (Connectivity)" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### ------------------------------------" << endl;

  std::string unite_arg = P.getOptionValue("-unite", "unite");
  std::string find_arg = P.getOptionValue("-find", "find_compress");

  std::string sampling_arg = P.getOptionValue("-sample", "none"); /* oneof { kout, bfs, none } */

  std::cout << "Params: -unite = " << unite_arg << " -find = " << find_arg << std::endl;


  timer t;
  t.start();
  pbbs::sequence<uintE> components;
  int sampling_num_rounds = P.getOptionLongValue("-sample_rounds", /*default rounds=*/2);
  if (sampling_arg == "kout") {
    if (unite_arg == "unite") {
      components = union_find::select_framework_algorithm<
        unite_variants::Unite,
        union_find::UnionFindSampleTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_early") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteEarly,
        union_find::UnionFindSampleTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_nd") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteND,
        union_find::UnionFindSampleTemplate,
        Graph>(G, find_arg, sampling_num_rounds, /* use_hooks = */true);
    } else if (unite_arg == "unite_rem") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteRemLock,
        union_find::UnionFindSampleTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else {
      std::cout << "Unknown unite variant: " << unite_arg << std::endl;
    }
  } else if (sampling_arg == "bfs") {
    if (unite_arg == "unite") {
      components = union_find::select_framework_algorithm<
        unite_variants::Unite,
        union_find::UnionFindSampledBFSTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_early") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteEarly,
        union_find::UnionFindSampledBFSTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_nd") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteND,
        union_find::UnionFindSampledBFSTemplate,
        Graph>(G, find_arg, sampling_num_rounds, /* use_hooks = */true);
    } else if (unite_arg == "unite_rem") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteRemLock,
        union_find::UnionFindSampledBFSTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else {
      std::cout << "Unknown unite variant: " << unite_arg << std::endl;
    }
  } else if (sampling_arg == "ldd") {
    if (unite_arg == "unite") {
      components = union_find::select_framework_algorithm<
        unite_variants::Unite,
        union_find::UnionFindLDDTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_early") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteEarly,
        union_find::UnionFindLDDTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_nd") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteND,
        union_find::UnionFindLDDTemplate,
        Graph>(G, find_arg, sampling_num_rounds, /* use_hooks = */true);
    } else if (unite_arg == "unite_rem") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteRemLock,
        union_find::UnionFindLDDTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else {
      std::cout << "Unknown unite variant: " << unite_arg << std::endl;
    }
  } else if (sampling_arg == "none") {
    if (unite_arg == "unite") {
      components = union_find::select_framework_algorithm<
        unite_variants::Unite,
        union_find::UnionFindTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_early") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteEarly,
        union_find::UnionFindTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_nd") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteND,
        union_find::UnionFindTemplate,
        Graph>(G, find_arg, sampling_num_rounds, /* use_hooks = */true);
    } else if (unite_arg == "unite_rem") {
      components = union_find::select_framework_algorithm<
        unite_variants::UniteRemLock,
        union_find::UnionFindTemplate,
        Graph>(G, find_arg, sampling_num_rounds, false);
    } else if (unite_arg == "unite_rank") {
      // Note that tht there are no options for alternate find implementations
      // in this algorithm, as it uses its own specialized find implementation.
      auto q = jayanti_rank::JayantiTBUnite<Graph>(G);
      components = q.components();
    } else {
      std::cout << "Unknown unite variant: " << unite_arg << std::endl;
    }
  } else {
    std::cout << "Unknown argument for sampling parameter: " << sampling_arg << std::endl;
  }
  double tt = t.stop();
  std::cout << "### Running Time: " << tt << std::endl;

  if (P.getOptionValue("-check")) {
    auto correct = workefficient_cc::CC(G, 0.2, false, true);
    RelabelDet(correct);
    cc_check(correct, components);
  }

  if (P.getOption("-stats")) {
    auto cc_f = [&](size_t i) { return components[i]; };
    auto cc_im =
        pbbslib::make_sequence<uintE>(G.n, cc_f);
    num_cc(cc_im);
    largest_cc(cc_im);
  }
  return tt;
}

generate_main(CC_runner, false);
