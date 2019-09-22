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
// numactl -i all ./Biconnectivity -s -m com-orkut.ungraph.txt_SJ
// flags:
//   required:
//     -s : indicate that the graph is symmetric
//   optional:
//     -c : indicate that the graph is compressed
//     -m : indicate that the graph should be mmap'd
//     -of : the output file to write the biconnectivity labels
//     -if : the file storing bicc labels, specify if you want to compute the
//           number of biconnected components.
//
// ex computing #biccs:
// > numactl -i all ./Biconnectivity -of orkut_bcs.out -s -m
// com-orkut.ungraph.txt_SJ > numactl -i all ./Biconnectivity -if orkut_bcs.out
// -s -m com-orkut.ungraph.txt_SJ

#include "Biconnectivity.h"
#include "ligra/pbbslib/sparse_additive_map.h"
#include "ligra/ligra.h"

template <template <typename W> class vertex, class W>
void BiconnectivityStats(graph<vertex<W>>& GA, char* s,
                         uintE component_id = UINT_E_MAX) {
  size_t n = GA.n;
  auto S = pbbslib::char_seq_from_file(s);
  auto tokens = pbbslib::tokenize(S, [] (const char c) { return pbbs::is_space(c); });
  auto labels = sequence<std::tuple<uintE, uintE>>(n);
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i) {
    labels[i] =
        std::make_tuple(atol(tokens[2 * i]), atol(tokens[2 * i + 1]));
  });

  auto bits = sequence<uintE>(n, (uintE)0);
  auto flags = sequence<bool>(n, false);

  auto bicc_label = [&](const uintE& u, const uintE& v) {
    auto lab_u = labels[u];
    auto lab_v = labels[v];
    uintE p_u = std::get<0>(lab_u);
    uintE p_v = std::get<0>(lab_v);
    if (v == p_u) {
      return std::get<1>(lab_u);
    } else if (u == p_v) {
      return std::get<1>(lab_v);
    } else {
      assert(std::get<1>(lab_u) == std::get<1>(lab_v));
      return std::get<1>(lab_u);
    }
    exit(0);
    return (uintE)1;
  };

  size_t mask = (1 << 12) - 1;
  auto empty = std::make_tuple(UINT_E_MAX, UINT_E_MAX);
  auto ST = sparse_additive_map<uintE, uintE>(n, empty);

  auto map_bc_label = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    auto label = bicc_label(src, ngh);
    if (!bits[label]) {
      bits[label] = 1;
    }
    size_t key = (static_cast<size_t>(src) << 32) + static_cast<size_t>(ngh);
    if ((pbbslib::hash64(key) & mask) == 0) {
      ST.insert(std::make_tuple(label, 1));
    }
    if (component_id != UINT_E_MAX) {
      if (label == component_id && !flags[src]) {
        flags[src] = 1;
      }
    }
  };
  par_for(0, n, pbbslib::kSequentialForThreshold, [&] (size_t i)
                  { GA.V[i].mapOutNgh(i, map_bc_label); });

  if (component_id == UINT_E_MAX) {
    auto ET = ST.entries();
    auto cmp_snd = [&](const std::tuple<uintE, uintE>& l,
                       const std::tuple<uintE, uintE>& r) {
      return std::get<1>(l) > std::get<1>(r);
    };
    pbbslib::sample_sort_inplace(ET.slice(), cmp_snd, true);
    for (size_t i = 0; i < std::min((size_t)10, ET.size()); i++) {
      std::cout << std::get<0>(ET[i]) << " " << std::get<1>(ET[i]) << "\n";
    }
  } else {
    // reduce flags
    auto flags_f = [&](size_t i) { return (size_t)flags[i]; };
    auto flags_imap = pbbslib::make_sequence<size_t>(n, flags_f);
    std::cout << "Largest component size = " << pbbslib::reduce_add(flags_imap)
              << "\n";
  }

  // The size of a biconnected component as the number of vertices that
  // participate in it. #edges is another option, but #vertices matches the
  // similar metrics for SCC and CC.
  // Can use sketching to find the top components, and then do an exact
  // component on the top-k using hash-maps, space usage would be O(nk).

  // Note that this is the number of biconnected components excluding isolated
  // vertices (the definition maps edges -> components, so isolated vertices
  // don't contribute to any meaningful components).
  uintE total_biccs = pbbslib::scan_add_inplace(bits);
  std::cout << "num biconnected components = " << total_biccs << "\n";
}

template <class vertex>
double Biconnectivity_runner(graph<vertex>& GA, commandLine P) {
  std::cout << "### Application: Biconnectivity" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << GA.n << std::endl;
  std::cout << "### m: " << GA.m << std::endl;
  std::cout << "### Params: n/a" << std::endl;
  std::cout << "### ------------------------------------" << endl;

  auto in_f = P.getOptionValue("-if");
  auto out_f = P.getOptionValue("-of");
  assert(P.getOptionValue("-s"));
  if (in_f) {
    BiconnectivityStats(GA, in_f);
  } else {
    timer t; t.start();
    Biconnectivity(GA, out_f);
    double tt = t.stop();
    std::cout << "### Running Time: " << tt << std::endl;
  }
  // Note that Biconnectivity mutates the graph, so we only run the algorithm
  // once.
  exit(0);
}

generate_main(Biconnectivity_runner, true);
